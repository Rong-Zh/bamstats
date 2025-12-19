
#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include "ArgumentParser.h"
#include "stats.h"

int main(int argc, char *argv[])
{
    // Create argument parser
    ArgumentParser parser("BAM file statistics tool");
    parser.set_program_name("bamstats");
    parser.set_version("1.1.0");
    // Add arguments
    parser.add_required<File>('b', "bam", "Input BAM file path");
    parser.add_required<File>('s', "stats", "Output statistics file path");
    parser.add_optional<int>('t', "thread", "Number of threads", 1);
    parser.add_optional<int>('q', "mapping-quality", "Mapping quality cutoff for unique mapping", 20);
    parser.add_optional<std::string>('L', "target", "Target regions: BED file path OR region string ");
    parser.add_optional<std::string>('d', "depth-thresholds", "Comma-separated depth thresholds for coverage calculation", "1,20,30,50,100,500,1000,3000");
    parser.add_help();

    parser.parse_args(argc, argv);

    // Get parsed values
    std::string bam_file = parser.get<File>("bam");
    std::string output_file = parser.get<File>("stats");
    int n_threads = parser.get<int>("thread");
    int mapq_cut = parser.get<int>("mapping-quality");
    std::string depth_thresholds_str = parser.get<std::string>("depth-thresholds");
    std::string target_spec = parser.get<std::string>("target");
    // Parse depth thresholds
    std::vector<uint32_t> depth_thresholds = parse_depth_thresholds(depth_thresholds_str);

    // Open BAM file first to get header (needed for region parsing)
    samFile *in = sam_open(bam_file.c_str(), "r");
    if (in == NULL)
    {
        std::cerr << "Error: Failed to open BAM file: " << bam_file << std::endl;
        return 1;
    }

    bam_hdr_t *header = sam_hdr_read(in);
    if (header == NULL)
    {
        std::cerr << "Error: Failed to read BAM header" << std::endl;
        sam_close(in);
        return 1;
    }

    // Check if index exists
    hts_idx_t *idx = sam_index_load(in, bam_file.c_str());
    if (idx == NULL)
    {
        std::cerr << "Error: Failed to load BAM index. Please index the BAM file first." << std::endl;
        bam_hdr_destroy(header);
        sam_close(in);
        return 1;
    }
    hts_idx_destroy(idx);

    // Load target regions (BED file or region string)
    std::unique_ptr<Bed> bed_regions;
    bool is_region_string = false;
    if (!target_spec.empty())
    {
        bed_regions = parse_region(target_spec, header);
        if (!bed_regions || bed_regions->is_empty())
        {
            std::cerr << "Error: Failed to parse target specification or no valid regions" << std::endl;
            bam_hdr_destroy(header);
            sam_close(in);
            return 1;
        }
        // Sort and merge overlapping intervals
        bed_regions->sort(n_threads);
        bed_regions->merge(n_threads);
        std::cout << "Target regions length: " << bed_regions->get_length() << " bp\n";
    }

    std::cout << "BAM Statistics Tool\n";
    std::cout << "===================\n";
    std::cout << "Input BAM: " << bam_file << "\n";
    std::cout << "Output file: " << output_file << "\n";
    std::cout << "Threads: " << n_threads << "\n";
    std::cout << "Mapping quality cutoff: " << mapq_cut << "\n";
    if (bed_regions)
    {
        if (is_region_string)
        {
            std::cout << "Target regions: " << target_spec << "\n";
        }
        else
        {
            std::cout << "BED file: " << target_spec << "\n";
        }
        std::cout << "Depth thresholds: ";
        for (size_t i = 0; i < depth_thresholds.size(); i++)
        {
            if (i > 0)
                std::cout << ", ";
            std::cout << depth_thresholds[i] << "x";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    // Get regions (chromosomes) for parallel processing
    std::vector<std::string> regions;
    if (bed_regions)
    {
        const auto &intervals = bed_regions->get_intervals();
        regions.reserve(intervals.size());
        for (const auto &[chrom, interval_list] : intervals)
        {
            (void)interval_list;
            regions.emplace_back(chrom);
        }
    }

    bam_hdr_destroy(header);
    sam_close(in);

    // Create result queue for threads to send back results
    // Larger queue reduces thread blocking
    constexpr uint32_t QUEUE_SIZE = 8192;
    ProducerConsumerQueue<ThreadStats> result_queue(QUEUE_SIZE);
    // Create tasks for each thread
    std::vector<RegionTask> tasks(n_threads);
    for (int i = 0; i < n_threads; ++i)
    {
        tasks[i].bam_path = bam_file;
        tasks[i].mapq_cut = mapq_cut;
        tasks[i].process_unmapped = (i == 0);
        tasks[i].bed_regions = bed_regions.get();
    }
    // Distribute chromosomes to threads in round-robin fashion
    for (size_t i = 0; i < regions.size(); ++i)
    {
        int thread_idx = i % n_threads;
        tasks[thread_idx].regions.push_back(regions[i]);
    }
    // Launch worker threads
    std::vector<std::thread> worker_threads;
    for (int i = 0; i < n_threads; ++i)
    {
        worker_threads.emplace_back(
            region_worker_thread,
            std::cref(tasks[i]),
            i,
            std::ref(result_queue));
    }
    // Wait for all worker threads to complete
    for (auto &thread : worker_threads)
    {
        thread.join();
    }

    std::cout << "All threads completed. Collecting results...\n";
    // Collect results from result queue
    BamStats global_stats;
    ThreadStats thread_stats;
    int results_collected = 0;
    while (result_queue.read(thread_stats))
    {
        global_stats.merge(thread_stats);
        results_collected++;
    }

    std::cout << "Collected results from " << results_collected << " threads.\n";
    std::cout << "Processing complete!\n";

    // Calculate depth statistics if BED file is provided
    DepthStats depth_stats;
    if (bed_regions)
    {
        std::cout << "\n=== Calculating Depth Statistics ===\n";
        std::cout << "BED region total length: " << bed_regions->get_length() << " bp\n";
        auto start_time = std::chrono::high_resolution_clock::now();

        depth_stats = calculate_depth_stats(
            bam_file, bed_regions.get(), mapq_cut, n_threads);

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        std::cout << "Depth calculation time: " << duration.count() << " seconds\n";
        std::cout << "Total bases counted: " << depth_stats.total_bases << " bp\n";
        std::cout << "Covered bases (>0x): " << depth_stats.covered_bases << " bp\n";
        std::cout << "Total rawdepth: " << depth_stats.total_rawdepth << "\n";
        std::cout << "Total rmdupdepth: " << depth_stats.total_rmdupdepth << "\n";
        std::cout << "Total covdepth: " << depth_stats.total_covdepth << "\n";
    }

    // Write statistics to file
    std::ofstream outfile(output_file);
    if (!outfile)
    {
        std::cerr << "Error: Cannot open output file: " << output_file << std::endl;
        return 1;
    }
    write_stats(outfile, global_stats, mapq_cut, bed_regions.get(),
                bed_regions ? &depth_stats : nullptr,
                bed_regions ? &depth_thresholds : nullptr);
    outfile.close();

    // Also write to stdout
    write_stats(std::cout, global_stats, mapq_cut, bed_regions.get(),
                bed_regions ? &depth_stats : nullptr,
                bed_regions ? &depth_thresholds : nullptr);
    std::cout << "\nStatistics written to: " << output_file << "\n";

    return 0;
}
