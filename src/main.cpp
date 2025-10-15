#include <iostream>
#include <fstream>
#include <string>
#include <atomic>
#include <vector>
#include <thread>
#include <mutex>
#include <memory>
#include <map>
#include <unordered_set>
#include <chrono>
#include <cstring>
#include <functional>
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "ProducerConsumerQueue.h"
#include "../ArgumentsParser/ArgumentParser.h"
#include "bed.h"

// 1: paired
// 2: properly paired
// 4: unmapped (BAM_FUNMAP)
// 8: mate unmapped
// 16: reverse strand
// 32: mate reverse strand
// 64: first in pair (BAM_FREAD1)
// 128: second in pair (BAM_FREAD2)
// 256: secondary alignment (BAM_FSECONDARY)
// 512: QC failed (BAM_FQCFAIL) ← 这个
// 1024: PCR/optical duplicate (BAM_FDUP)
// 2048: supplementary alignment (BAM_FSUPPLEMENTARY)

// Get alignment start and end positions (excluding soft clips)
// This matches bamdst's behavior
void get_alignment_span(bam1_t *b, int32_t& aln_start, int32_t& aln_end) {
    aln_start = b->core.pos;  // 0-based leftmost position
    
    uint32_t *cigar = bam_get_cigar(b);
    int n_cigar = b->core.n_cigar;
    
    // Skip leading soft clips to get true alignment start
    int32_t pos = aln_start;
    int i = 0;
    while (i < n_cigar && bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
        i++;
    }
    
    // Calculate alignment end by walking through CIGAR
    for (; i < n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        
        if (bam_cigar_type(op) & 2) {  // Consumes reference
            pos += len;
        }
    }
    
    aln_end = pos;  // 0-based, exclusive (position after last aligned base)
}

// Structure to hold statistics for each thread
struct ThreadStats {
    uint64_t total_records = 0;
    uint64_t total_reads = 0;
    uint64_t qc_failed = 0;            // QC failed(flag 512)
    uint64_t duplicates = 0;           // PCR/optical duplicate (flag 1024)
    uint64_t non_primary = 0;
    uint64_t secondary = 0;            // Secondary alignments (flag 256)
    uint64_t supplementary = 0;        // Supplementary alignments (flag 2048)
    uint64_t unmapped = 0;             // Unmapped reads (flag 4)
    uint64_t mapped = 0;               // Mapped reads (primary alignments that are not unmapped)
    uint64_t mapq_below_cut = 0;
    uint64_t mapq_above_cut = 0;
    uint64_t read1 = 0;
    uint64_t read2 = 0;
    uint64_t plus_strand = 0;
    uint64_t minus_strand = 0;
    uint64_t non_splice = 0;
    uint64_t splice = 0;
    uint64_t proper_pairs = 0;
    uint64_t reads_mate_diff_chrom = 0;  // All reads where mate maps to different chromosome
    uint64_t on_target = 0;  // Reads overlapping BED regions (only counted when BED file is provided)
};

// Structure to hold final statistics
struct BamStats {
    std::atomic<uint64_t> total_records{0};
    std::atomic<uint64_t> total_reads{0};
    std::atomic<uint64_t> qc_failed{0};
    std::atomic<uint64_t> duplicates{0};
    std::atomic<uint64_t> non_primary{0};
    std::atomic<uint64_t> secondary{0};
    std::atomic<uint64_t> supplementary{0};
    std::atomic<uint64_t> unmapped{0};
    std::atomic<uint64_t> mapped{0};
    std::atomic<uint64_t> mapq_below_cut{0};
    std::atomic<uint64_t> mapq_above_cut{0};
    std::atomic<uint64_t> read1{0};
    std::atomic<uint64_t> read2{0};
    std::atomic<uint64_t> plus_strand{0};
    std::atomic<uint64_t> minus_strand{0};
    std::atomic<uint64_t> non_splice{0};
    std::atomic<uint64_t> splice{0};
    std::atomic<uint64_t> proper_pairs{0};
    std::atomic<uint64_t> proper_pairs_diff_chrom{0};
    std::atomic<uint64_t> reads_mate_diff_chrom{0};  // All reads where mate maps to different chromosome
    std::atomic<uint64_t> on_target{0};  // Reads overlapping BED regions
    
    // Merge thread statistics into global statistics
    void merge(const ThreadStats& stats) {
        total_records += stats.total_records;
        total_reads += stats.total_reads;
        qc_failed += stats.qc_failed;
        duplicates += stats.duplicates;
        non_primary += stats.non_primary;
        secondary += stats.secondary;
        supplementary += stats.supplementary;
        unmapped += stats.unmapped;
        mapped += stats.mapped;
        mapq_below_cut += stats.mapq_below_cut;
        mapq_above_cut += stats.mapq_above_cut;
        read1 += stats.read1;
        read2 += stats.read2;
        plus_strand += stats.plus_strand;
        minus_strand += stats.minus_strand;
        non_splice += stats.non_splice;
        splice += stats.splice;
        proper_pairs += stats.proper_pairs;
        reads_mate_diff_chrom += stats.reads_mate_diff_chrom;
        on_target += stats.on_target;
    }
};

// Task structure for processing (region-based approach)
struct RegionTask {
    std::string bam_path;
    int mapq_cut;
    std::vector<std::string> regions;  // List of chromosome names or regions (e.g., "chr1", "chr2")
    bool process_unmapped;  // Whether this thread should process unmapped reads
    const Bed* bed_regions;  // Pointer to BED regions (shared across threads), nullptr if no BED file
};

// Function to check if read is spliced (has 'N' in CIGAR)
bool is_spliced(bam1_t *b) {
    uint32_t *cigar = bam_get_cigar(b);
    int n_cigar = b->core.n_cigar;
    
    for (int i = 0; i < n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        if (op == BAM_CREF_SKIP) {  // 'N' operation
            return true;
        }
    }
    return false;
}

// Process statistics for a single read
void process_single_read(bam1_t *b, int mapq_cut, const bam_hdr_t* header, 
                         const Bed* bed_regions, ThreadStats& stats) {
    stats.total_records++;
    
    // Non-primary alignment - skip for read-level statistics
    bool is_secondary = (b->core.flag & BAM_FSECONDARY);
    bool is_supplementary = (b->core.flag & BAM_FSUPPLEMENTARY);
    if (is_secondary || is_supplementary) {
        stats.non_primary++;
        stats.secondary += is_secondary;
        stats.supplementary += is_supplementary;
        return;  // Skip non-primary for read-level stats
    }
    
    // From here on, we only count primary alignments + unmapped reads
    // This is read-level statistics (each read counted once)
    stats.total_reads++;
    
    // QC failed (read-level)
    if (b->core.flag & BAM_FQCFAIL) {
        stats.qc_failed++;
    }
    
    // Optical/PCR duplicate (read-level)
    if (b->core.flag & BAM_FDUP) {
        stats.duplicates++;
    }
    
    // Unmapped reads
    if (b->core.flag & BAM_FUNMAP) {
        stats.unmapped++;
        return;  // Skip unmapped for other stats
    }
    
    // This is a mapped primary read
    stats.mapped++;
    
    // Check if read and mate map to different chromosomes (for all paired reads, not just proper pairs)
    // tid: chromosome ID of current read, mtid: chromosome ID of mate
    // mtid >= 0 means mate is mapped
    if ((b->core.flag & BAM_FPAIRED) && b->core.mtid >= 0 && b->core.tid != b->core.mtid) {
        stats.reads_mate_diff_chrom++;
    }
    
    // Check if read is on target (BED regions)
    // Use alignment span (excluding soft clips) to match bamdst behavior
    if (bed_regions != nullptr && !bed_regions->is_empty()) {
        const char* chrom = header->target_name[b->core.tid];
        int32_t aln_start, aln_end;
        get_alignment_span(b, aln_start, aln_end);
        
        if (bed_regions->is_ontarget(chrom, aln_start, aln_end)) {
            stats.on_target++;
        }
    }
    
    // Proper pairs (count all mapped primary reads, regardless of MAPQ)
    if (b->core.flag & BAM_FPROPER_PAIR) {
        stats.proper_pairs++;
    }
    
    // MAPQ filtering
    if (b->core.qual < mapq_cut) {
        stats.mapq_below_cut++;
    } else {
        stats.mapq_above_cut++;
        
        // Read 1 or Read 2
        if (b->core.flag & BAM_FREAD1) {
            stats.read1++;
        }
        if (b->core.flag & BAM_FREAD2) {
            stats.read2++;
        }
        
        // Strand
        if (b->core.flag & BAM_FREVERSE) {
            stats.minus_strand++;
        } else {
            stats.plus_strand++;
        }
        
        // Splice reads
        if (is_spliced(b)) {
            stats.splice++;
        } else {
            stats.non_splice++;
        }
    }
}

// Consumer thread function: processes specific chromosome regions in parallel
void region_worker_thread(
    const RegionTask& task,
    int thread_id,
    ProducerConsumerQueue<ThreadStats>& result_queue)
{
    // Open BAM file for this thread
    samFile *in = sam_open(task.bam_path.c_str(), "r");
    if (in == NULL) {
        std::cerr << "Thread " << thread_id << ": Failed to open BAM file" << std::endl;
        return;
    }
    
    // Enable htslib threading for decompression
    hts_set_threads(in, 1);
    
    bam_hdr_t *header = sam_hdr_read(in);
    if (header == NULL) {
        std::cerr << "Thread " << thread_id << ": Failed to read header" << std::endl;
        sam_close(in);
        return;
    }
    
    // Load BAM index
    hts_idx_t *idx = sam_index_load(in, task.bam_path.c_str());
    if (idx == NULL) {
        std::cerr << "Thread " << thread_id << ": Failed to load BAM index. Make sure .bai file exists." << std::endl;
        bam_hdr_destroy(header);
        sam_close(in);
        return;
    }
    
    ThreadStats local_stats;
    bam1_t *b = bam_init1();
    
    // Process each region assigned to this thread
    for (const std::string& region : task.regions) {
        hts_itr_t *iter = sam_itr_querys(idx, header, region.c_str());
        if (iter == NULL) {
            std::cerr << "Thread " << thread_id << ": Failed to parse region: " << region << std::endl;
            continue;
        }
        
        // Iterate through reads in this region
        while (sam_itr_next(in, iter, b) >= 0) {
            process_single_read(b, task.mapq_cut, header, task.bed_regions, local_stats);
        }
        
        hts_itr_destroy(iter);
    }
    
    // Process unmapped reads if assigned to this thread
    if (task.process_unmapped) {
        // Query for unmapped reads using special region "*"
        hts_itr_t *iter = sam_itr_querys(idx, header, "*");
        if (iter != NULL) {
            while (sam_itr_next(in, iter, b) >= 0) {
                process_single_read(b, task.mapq_cut, header, task.bed_regions, local_stats);
            }
            hts_itr_destroy(iter);
        }
    }
    
    // Send results back
    bool sent = false;
    while (!sent) {
        sent = result_queue.write(local_stats);
        if (!sent) {
            std::this_thread::yield();
        }
    }
    
    // Cleanup
    bam_destroy1(b);
    hts_idx_destroy(idx);
    bam_hdr_destroy(header);
    sam_close(in);
}

void write_stats(std::ostream& out, const BamStats& stats, int mapq_cut, const Bed* bed) {
    out << "Total records: " << stats.total_records << "\n";
    out << "Total reads: " << stats.total_reads << "\n";
    out << "QC failed: " << stats.qc_failed << "\n";
    out << "Optical/PCR duplicate: " << stats.duplicates << "\n";
    out << "Non primary hits: " << stats.non_primary << "\n";
    out << "Secondary alignments: " << stats.secondary << "\n";
    out << "Supplementary alignments: " << stats.supplementary << "\n";
    out << "Unmapped reads: " << stats.unmapped << "\n";
    out << "Mapped reads: " << stats.mapped << "\n";
    out << "Mapping quality < cutoff (non-unique): " << stats.mapq_below_cut << "\n";
    out << "Mapping quality >= cutoff (unique): " << stats.mapq_above_cut << "\n";
    out << "Read-1: " << stats.read1 << "\n";
    out << "Read-2: " << stats.read2 << "\n";
    out << "Reads map to '+': " << stats.plus_strand << "\n";
    out << "Reads map to '-': " << stats.minus_strand << "\n";
    out << "Non-splice reads: " << stats.non_splice << "\n";
    out << "Splice reads: " << stats.splice << "\n";
    out << "Proper pairs: " << stats.proper_pairs << "\n";
    out << "Reads and mate map to different chrom: " << stats.reads_mate_diff_chrom << "\n";
    if (bed) {
        out << "Target region size(bp): " << bed->get_length() << "\n";
        out << "Ontarget reads:" << stats.on_target << "\n";
    }
}

int main(int argc, char *argv[]) {
    // Create argument parser
    ArgumentParser parser("BAM file statistics tool");
    parser.set_program_name("bamstats");
    parser.set_version("1.0.0");
    
    // Add arguments
    parser.add_required<File>('b', "bam", "Input BAM file path");
    parser.add_required<File>('s', "stats", "Output statistics file path");
    parser.add_optional<int>('t', "thread", "Number of threads", 1);
    parser.add_optional<int>('q', "mapping-quality", "Mapping quality cutoff for unique mapping", 30);
    parser.add_optional<File>('L', "target", "BED file with target regions");
    parser.add_flag('h', "help", "Show this help message");
    
    // If no arguments provided, print help and exit
    if (argc == 1) {
        parser.print_help();
        return 0;
    }
    // Parse arguments
    try {
        parser.parse_args(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "Error parsing arguments: " << e.what() << std::endl;
        parser.print_help();
        return 1;
    }
    // Get parsed values
    std::string bam_file = parser.get<File>("bam");
    std::string output_file = parser.get<File>("stats");
    int n_threads = parser.get<int>("thread");
    int mapq_cut = parser.get<int>("mapping-quality");
    bool show_help = parser.get<bool>("help");
    if (show_help) {
        parser.print_help();
        return 0;
    }
    
    // Load BED file if provided
    std::unique_ptr<Bed> bed_regions;
    std::string bed_file = parser.get<File>("target");
    if (!bed_file.empty()) {
        std::cout << "Loading BED file: " << bed_file << "\n";
        bed_regions = std::make_unique<Bed>(bed_file);
        if (bed_regions->is_empty()) {
            std::cerr << "Error: Failed to load BED file or BED file is empty" << std::endl;
            return 1;
        }
        // Sort and merge overlapping intervals to remove duplicates
        bed_regions->sort(n_threads);
        bed_regions->merge(n_threads);
        std::cout << "BED regions length: " << bed_regions->get_length() << " bp\n";
    }
    
    std::cout << "BAM Statistics Tool\n";
    std::cout << "===================\n";
    std::cout << "Input BAM: " << bam_file << "\n";
    std::cout << "Output file: " << output_file << "\n";
    std::cout << "Threads: " << n_threads << "\n";
    std::cout << "Mapping quality cutoff: " << mapq_cut << "\n";
    if (bed_regions) {
        std::cout << "BED file: " << bed_file << "\n";
    }
    std::cout << "\n";
    // Open BAM file to get header and regions
    samFile *in = sam_open(bam_file.c_str(), "r");
    if (in == NULL) {
        std::cerr << "Error: Failed to open BAM file: " << bam_file << std::endl;
        return 1;
    }
    
    bam_hdr_t *header = sam_hdr_read(in);
    if (header == NULL) {
        std::cerr << "Error: Failed to read BAM header" << std::endl;
        sam_close(in);
        return 1;
    }
    
    // Check if index exists
    hts_idx_t *idx = sam_index_load(in, bam_file.c_str());
    if (idx == NULL) {
        std::cerr << "Error: Failed to load BAM index. Please index the BAM file first." << std::endl;
        bam_hdr_destroy(header);
        sam_close(in);
        return 1;
    }
    hts_idx_destroy(idx);
    
    // Get regions (chromosomes) for parallel processing
    std::vector<std::string> regions;
    for (int i = 0; i < header->n_targets; i++) {
        regions.emplace_back(std::string(header->target_name[i]));
    }
    bam_hdr_destroy(header);
    sam_close(in);
    
    // Create result queue for threads to send back results
    constexpr uint32_t QUEUE_SIZE = 128;
    ProducerConsumerQueue<ThreadStats> result_queue(QUEUE_SIZE);
    // Create tasks for each thread
    std::vector<RegionTask> tasks(n_threads);
    for (int i = 0; i < n_threads; ++i) {
        tasks[i].bam_path = bam_file;
        tasks[i].mapq_cut = mapq_cut;
        tasks[i].process_unmapped = (i == 0);  // Only thread 0 processes unmapped reads
        tasks[i].bed_regions = bed_regions.get();  // Pass pointer to BED regions (nullptr if not provided)
    }
    // Distribute chromosomes to threads in round-robin fashion
    for (size_t i = 0; i < regions.size(); ++i) {
        int thread_idx = i % n_threads;
        tasks[thread_idx].regions.push_back(regions[i]);
    }
    // Launch worker threads
    std::vector<std::thread> worker_threads;
    for (int i = 0; i < n_threads; ++i) {
        worker_threads.emplace_back(
            region_worker_thread,
            std::cref(tasks[i]),
            i,
            std::ref(result_queue)
        );
    }
    // Wait for all worker threads to complete
    for (auto& thread : worker_threads) {
        thread.join();
    } 
    std::cout << "All threads completed. Collecting results...\n";
    // Collect results from result queue
    BamStats global_stats;
    ThreadStats thread_stats;
    int results_collected = 0;
    while (result_queue.read(thread_stats)) {
        global_stats.merge(thread_stats);
        results_collected++;
    }
    std::cout << "Collected results from " << results_collected << " threads.\n";
    std::cout << "Processing complete!\n";
    
    // Write statistics to file
    std::ofstream outfile(output_file);
    if (!outfile) {
        std::cerr << "Error: Cannot open output file: " << output_file << std::endl;
        return 1;
    }
    write_stats(outfile, global_stats, mapq_cut, bed_regions.get());
    outfile.close();
    
    // Also write to stdout
    write_stats(std::cout, global_stats, mapq_cut, bed_regions.get());
    std::cout << "\nStatistics written to: " << output_file << "\n";
    return 0;
}
