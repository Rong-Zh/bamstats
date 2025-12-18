#ifndef STATS_H
#define STATS_H

#include <algorithm>
#include <atomic>
#include <charconv>
#include <chrono>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <string_view>
#include <system_error>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "htslib/hts.h"
#include "htslib/sam.h"

#include "ProducerConsumerQueue.h"
#include "bed.h"

inline void trim_inplace(std::string &s)
{
    if (s.empty())
        return;
    
    size_t left = 0;
    size_t right = s.size() - 1;
    
    while (left <= right && (s[left] == ' ' || s[left] == '\t' || s[left] == '\r' || s[left] == '\n'))
        ++left;
    
    if (left > right)
    {
        s.clear();
        return;
    }
    
    while (right > left && (s[right] == ' ' || s[right] == '\t' || s[right] == '\r' || s[right] == '\n'))
        --right;
    
    if (left > 0 || right < s.size() - 1)
    {
        s = s.substr(left, right - left + 1);
    }
}

inline bool to_u32(std::string s, uint32_t &out)
{
    auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), out);
    return ec == std::errc() && ptr == s.data() + s.size();
}

inline bool to_u64(std::string s, uint64_t &out)
{
    auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), out);
    return ec == std::errc() && ptr == s.data() + s.size();
}

inline bool file_exists(const std::string &filename)
{
    return std::filesystem::exists(filename);
}

inline std::unique_ptr<Bed> parse_region(const std::string &target_spec, bam_hdr_t *header)
{
    if (target_spec.empty())
        return nullptr;
    if (file_exists(target_spec)) 
        return std::make_unique<Bed>(target_spec, false);
    
        // It's a region string - parse it
    if (!header)
        return nullptr;
    std::unordered_map<std::string, std::vector<Interval>> intervals;
    std::vector<std::string> parts;
    parts.reserve(16);
    split(target_spec, ',', parts);
    for (std::string &raw : parts)
    {
        trim_inplace(raw);
        if (raw.empty())
            continue;
        const size_t colon = raw.find(':');
        if (colon == std::string::npos)
        {
            // Chromosome name only - get full length from BAM header
            const std::string chrom = raw;
            int tid = bam_name2id(header, chrom.c_str());
            if (tid < 0)
            {
                std::cerr << "Warning: unknown chromosome in -L: '" << chrom << "' (ignored)\n";
                continue;
            }
            uint64_t len = static_cast<uint64_t>(header->target_len[tid]);
            if (len == 0)
                continue;
            // 0-based half-open: [0, len)
            intervals[chrom].push_back(Interval{0, len});
            continue;
        }
        // chr:start-end format (1-based inclusive)
        const std::string chrom = raw.substr(0, colon);
        std::string coord = raw.substr(colon + 1);
        trim_inplace(coord);
        const size_t dash = coord.find('-');
        if (dash == std::string::npos)
        {
            std::cerr << "Warning: invalid interval (expected chr:start-end): '" << raw << "' (ignored)\n";
            continue;
        }

        std::string start_s = coord.substr(0, dash);
        std::string end_s = coord.substr(dash + 1);
        trim_inplace(start_s);
        trim_inplace(end_s);
        uint64_t start_1based = 0;
        uint64_t end_1based = 0;
        if (!to_u64(start_s, start_1based) || 
            !to_u64(end_s, end_1based) || 
            start_1based == 0 || start_1based > end_1based)
        {
            std::cerr << "Warning: invalid interval: '" << raw << "' (ignored)\n";
            continue;
        }
        // Convert 1-based inclusive to 0-based half-open: [start-1, end)
        intervals[chrom].push_back(Interval{start_1based - 1, end_1based});
    }
    return std::make_unique<Bed>(std::move(intervals), false);
}

// Get alignment start and end positions (excluding soft clips)
// This matches bamdst's behavior
inline void get_alignment_span(bam1_t *b, int32_t &aln_start, int32_t &aln_end)
{
    aln_start = b->core.pos; // 0-based leftmost position

    uint32_t *cigar = bam_get_cigar(b);
    int n_cigar = b->core.n_cigar;

    // Skip leading soft clips to get true alignment start
    int32_t pos = aln_start;
    int i = 0;
    while (i < n_cigar && bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP)
    {
        i++;
    }

    // Calculate alignment end by walking through CIGAR
    for (; i < n_cigar; i++)
    {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);

        if (bam_cigar_type(op) & 2)
        { // Consumes reference
            pos += len;
        }
    }

    aln_end = pos; // 0-based, exclusive (position after last aligned base)
}

// Structure to hold depth statistics for BED regions
// Follows bamdst's three-depth model
struct DepthStats
{
    uint64_t total_bases = 0;                               // Total bases in BED regions
    uint64_t covered_bases = 0;                             // Bases with covdepth > 0
    uint64_t total_rawdepth = 0;                            // Sum of raw depths (with duplicates)
    uint64_t total_rmdupdepth = 0;                          // Sum of rmdup depths (without duplicates)
    uint64_t total_covdepth = 0;                            // Sum of coverage depths (with deletions)
    std::unordered_map<uint32_t, uint64_t> depth_histogram; // covdepth -> count

    double get_coverage() const
    {
        return total_bases > 0 ? (double)covered_bases / total_bases * 100.0 : 0.0;
    }

    // bamdst "Average depth" uses coverage depth: covdep
    double get_mean_depth() const
    {
        return total_bases > 0 ? (double)total_covdepth / total_bases : 0.0;
    }

    // Match bamdst "Average depth(rmdup)"
    double get_mean_depth_rmdup() const
    {
        return total_bases > 0 ? (double)total_rmdupdepth / total_bases : 0.0;
    }

    // Match bamdst behavior: avg is integer (floor), then threshold is floor(fraction * avg).
    double get_uniformity(double fraction = 0.2) const
    {
        if (total_bases == 0)
            return 0.0;

        uint64_t avg_int = total_covdepth / total_bases;
        if (avg_int == 0)
            return 0.0;

        uint32_t threshold = (uint32_t)(avg_int * fraction);
        uint64_t bases_above = 0;

        for (const auto &[depth, count] : depth_histogram)
        {
            if (depth >= threshold)
            {
                bases_above += count;
            }
        }

        return (double)bases_above / total_bases * 100.0;
    }

    double get_coverage_at_depth(uint32_t threshold) const
    {
        if (total_bases == 0)
            return 0.0;

        uint64_t bases_above = 0;
        for (const auto &[depth, count] : depth_histogram)
        {
            if (depth >= threshold)
            {
                bases_above += count;
            }
        }

        return (double)bases_above / total_bases * 100.0;
    }

    void merge(const DepthStats &other)
    {
        total_bases += other.total_bases;
        covered_bases += other.covered_bases;
        total_rawdepth += other.total_rawdepth;
        total_rmdupdepth += other.total_rmdupdepth;
        total_covdepth += other.total_covdepth;

        if (depth_histogram.size() + other.depth_histogram.size() > depth_histogram.bucket_count())
        {
            depth_histogram.reserve(depth_histogram.size() + other.depth_histogram.size());
        }

        for (const auto &[depth, count] : other.depth_histogram)
        {
            depth_histogram[depth] += count;
        }
    }
};

struct ThreadStats
{
    uint64_t total_records = 0;
    uint64_t total_reads = 0;
    uint64_t qc_failed = 0;  // QC failed(flag 512)
    uint64_t duplicates = 0; // PCR/optical duplicate (flag 1024)
    uint64_t non_primary = 0;
    uint64_t secondary = 0;     // Secondary alignments (flag 256)
    uint64_t supplementary = 0; // Supplementary alignments (flag 2048)
    uint64_t unmapped = 0;      // Unmapped reads (flag 4)
    uint64_t mapped = 0;        // Mapped reads (primary alignments that are not unmapped)
    uint64_t mapq_below_cut = 0;
    uint64_t mapq_above_cut = 0;
    uint64_t read1 = 0;
    uint64_t read2 = 0;
    uint64_t plus_strand = 0;
    uint64_t minus_strand = 0;
    uint64_t non_splice = 0;
    uint64_t splice = 0;
    uint64_t proper_pairs = 0;
    uint64_t reads_mate_diff_chrom = 0; // All reads where mate maps to different chromosome
    uint64_t on_target = 0;             // Reads overlapping BED regions (only counted when BED file is provided)
};

struct BamStats
{
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
    std::atomic<uint64_t> reads_mate_diff_chrom{0};
    std::atomic<uint64_t> on_target{0};

    void merge(const ThreadStats &stats)
    {
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

struct RegionTask
{
    std::string bam_path;
    int mapq_cut;
    std::vector<std::string> regions;
    bool process_unmapped;
    const Bed *bed_regions;
};

inline bool is_spliced(bam1_t *b)
{
    uint32_t *cigar = bam_get_cigar(b);
    int n_cigar = b->core.n_cigar;

    for (int i = 0; i < n_cigar; i++)
    {
        int op = bam_cigar_op(cigar[i]);
        if (op == BAM_CREF_SKIP)
        {
            return true;
        }
    }
    return false;
}

inline void process_single_read(bam1_t *b, int mapq_cut, const bam_hdr_t *header,
                                const Bed *bed_regions, ThreadStats &stats)
{
    stats.total_records++;

    bool is_secondary = (b->core.flag & BAM_FSECONDARY);
    bool is_supplementary = (b->core.flag & BAM_FSUPPLEMENTARY);
    if (is_secondary || is_supplementary)
    {
        stats.non_primary++;
        if (is_secondary)
            stats.secondary++;
        if (is_supplementary)
            stats.supplementary++;
        return;
    }

    stats.total_reads++;

    if (b->core.flag & BAM_FQCFAIL)
    {
        stats.qc_failed++;
    }

    if (b->core.flag & BAM_FDUP)
    {
        stats.duplicates++;
    }

    if (b->core.flag & BAM_FUNMAP)
    {
        stats.unmapped++;
        return;
    }

    stats.mapped++;

    if ((b->core.flag & BAM_FPAIRED) && b->core.mtid >= 0 && b->core.tid != b->core.mtid)
    {
        stats.reads_mate_diff_chrom++;
    }

    if (bed_regions != nullptr && !bed_regions->is_empty())
    {
        int32_t aln_start, aln_end;
        get_alignment_span(b, aln_start, aln_end);

        const char *chrom_name = header->target_name[b->core.tid];
        if (bed_regions->is_ontarget(chrom_name, aln_start, aln_end))
        {
            stats.on_target++;
        }
    }

    if (b->core.flag & BAM_FPROPER_PAIR)
    {
        stats.proper_pairs++;
    }

    if (b->core.qual < mapq_cut)
    {
        stats.mapq_below_cut++;
    }
    else
    {
        stats.mapq_above_cut++;
    }

    if (b->core.flag & BAM_FREAD1)
        stats.read1++;
    if (b->core.flag & BAM_FREAD2)
        stats.read2++;

    if (b->core.flag & BAM_FREVERSE)
    {
        stats.minus_strand++;
    }
    else
    {
        stats.plus_strand++;
    }

    if (is_spliced(b))
    {
        stats.splice++;
    }
    else
    {
        stats.non_splice++;
    }
}

inline DepthStats calculate_region_depth(
    samFile *in,
    bam_hdr_t *header,
    hts_idx_t *idx,
    const std::string &chrom,
    int64_t start,
    int64_t end,
    int mapq_cut)
{
    DepthStats stats;

    int tid = bam_name2id(header, chrom.c_str());
    if (tid < 0)
        return stats;

    int64_t region_len = end - start;
    std::vector<uint32_t> depth_array(region_len * 3, 0);

    hts_itr_t *iter = sam_itr_queryi(idx, tid, start, end);
    if (!iter)
        return stats;

    bam1_t *b = bam_init1();

    while (sam_itr_next(in, iter, b) >= 0)
    {
        if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY))
            continue;

        bool is_qcfail = (b->core.flag & BAM_FQCFAIL);
        bool is_duplicate = (!is_qcfail) && (b->core.flag & BAM_FDUP);
        bool pass_mapq = (b->core.qual >= mapq_cut);

        uint32_t *cigar = bam_get_cigar(b);
        int n_cigar = b->core.n_cigar;
        int64_t ref_pos = b->core.pos;

        for (int i = 0; i < n_cigar; i++)
        {
            int op = bam_cigar_op(cigar[i]);
            int len = bam_cigar_oplen(cigar[i]);

            if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF)
            {
                int64_t seg_start = std::max(ref_pos, start);
                int64_t seg_end = std::min(ref_pos + len, end);

                if (seg_start < seg_end)
                {
                    for (int64_t pos = seg_start; pos < seg_end; pos++)
                    {
                        int64_t offset = pos - start;
                        if (offset >= 0 && offset < region_len)
                        {
                            depth_array[offset * 3 + 0]++;
                            if (!is_duplicate && pass_mapq)
                                depth_array[offset * 3 + 1]++;
                            depth_array[offset * 3 + 2]++;
                        }
                    }
                }
                ref_pos += len;
            }
            else if (op == BAM_CDEL)
            {
                int64_t seg_start = std::max(ref_pos, start);
                int64_t seg_end = std::min(ref_pos + len, end);

                if (seg_start < seg_end)
                {
                    for (int64_t pos = seg_start; pos < seg_end; pos++)
                    {
                        int64_t offset = pos - start;
                        if (offset >= 0 && offset < region_len)
                        {
                            depth_array[offset * 3 + 2]++;
                        }
                    }
                }
                ref_pos += len;
            }
            else if (bam_cigar_type(op) & 2)
            {
                ref_pos += len;
            }
        }
    }

    bam_destroy1(b);
    hts_itr_destroy(iter);

    stats.total_bases = region_len;

    for (int64_t i = 0; i < region_len; i++)
    {
        uint32_t rawdepth = depth_array[i * 3 + 0];
        uint32_t rmdupdepth = depth_array[i * 3 + 1];
        uint32_t covdepth = depth_array[i * 3 + 2];

        if (covdepth > 0)
        {
            stats.covered_bases++;
        }

        stats.total_rawdepth += rawdepth;
        stats.total_rmdupdepth += rmdupdepth;
        stats.total_covdepth += covdepth;
        stats.depth_histogram[covdepth]++;
    }

    return stats;
}

inline void depth_worker_thread(
    const std::string &bam_path,
    const std::vector<std::tuple<std::string, int64_t, int64_t>> &regions,
    int mapq_cut,
    int thread_id,
    ProducerConsumerQueue<DepthStats> &result_queue)
{
    samFile *in = sam_open(bam_path.c_str(), "r");
    if (!in)
    {
        std::cerr << "Thread " << thread_id << ": Failed to open BAM file for depth calculation\n";
        return;
    }

    hts_set_threads(in, 1);
    bam_hdr_t *header = sam_hdr_read(in);
    if (!header)
    {
        sam_close(in);
        return;
    }

    hts_idx_t *idx = sam_index_load(in, bam_path.c_str());
    if (!idx)
    {
        bam_hdr_destroy(header);
        sam_close(in);
        return;
    }

    DepthStats local_stats;
    for (const auto &[chrom, start, end] : regions)
    {
        DepthStats region_stats = calculate_region_depth(
            in, header, idx, chrom, start, end, mapq_cut);
        local_stats.merge(region_stats);
    }

    bool sent = false;
    while (!sent)
    {
        sent = result_queue.write(local_stats);
        if (!sent)
        {
            std::this_thread::yield();
        }
    }

    hts_idx_destroy(idx);
    bam_hdr_destroy(header);
    sam_close(in);
}

inline DepthStats calculate_depth_stats(
    const std::string &bam_path,
    const Bed *bed_regions,
    int mapq_cut,
    int n_threads)
{
    DepthStats final_stats;

    if (!bed_regions || bed_regions->is_empty())
    {
        return final_stats;
    }

    const auto &intervals = bed_regions->get_intervals();

    size_t total_regions = 0;
    for (const auto &[chrom, interval_list] : intervals)
    {
        total_regions += interval_list.size();
    }

    std::vector<std::tuple<std::string, int64_t, int64_t>> all_regions;
    all_regions.reserve(total_regions);

    for (const auto &[chrom, interval_list] : intervals)
    {
        for (const auto &interval : interval_list)
        {
            all_regions.emplace_back(chrom, interval.start, interval.end);
        }
    }

    if (all_regions.empty())
    {
        return final_stats;
    }

    std::sort(all_regions.begin(), all_regions.end(),
              [](const auto &a, const auto &b)
              {
                  int64_t len_a = std::get<2>(a) - std::get<1>(a);
                  int64_t len_b = std::get<2>(b) - std::get<1>(b);
                  return len_a > len_b;
              });

    std::vector<std::vector<std::tuple<std::string, int64_t, int64_t>>> thread_regions(n_threads);
    std::vector<int64_t> thread_load(n_threads, 0);

    for (const auto &region : all_regions)
    {
        int64_t region_size = std::get<2>(region) - std::get<1>(region);
        int min_thread = (int)(std::min_element(thread_load.begin(), thread_load.end()) - thread_load.begin());
        thread_regions[min_thread].push_back(region);
        thread_load[min_thread] += region_size;
    }

    ProducerConsumerQueue<DepthStats> result_queue(8192);

    std::vector<std::thread> workers;
    for (int i = 0; i < n_threads; i++)
    {
        if (!thread_regions[i].empty())
        {
            workers.emplace_back(
                depth_worker_thread,
                std::cref(bam_path),
                std::cref(thread_regions[i]),
                mapq_cut,
                i,
                std::ref(result_queue));
        }
    }

    for (auto &thread : workers)
    {
        thread.join();
    }

    DepthStats thread_stats;
    while (result_queue.read(thread_stats))
    {
        final_stats.merge(thread_stats);
    }

    return final_stats;
}

inline std::vector<uint32_t> parse_depth_thresholds(const std::string &str)
{
    std::vector<uint32_t> thresholds;
    std::vector<std::string> items;
    items.reserve(16);
    split(str, ',', items);

    for (std::string &item : items)
    {
        trim_inplace(item);
        uint32_t val = 0;
        if (!to_u32(item, val))
        {
            std::cerr << "Warning: Invalid depth threshold '" << item << "', ignoring\n";
            continue;
        }
        thresholds.push_back(val);
    }
    std::sort(thresholds.begin(), thresholds.end());
    thresholds.erase(std::unique(thresholds.begin(), thresholds.end()), thresholds.end());
    return thresholds;
}

inline void region_worker_thread(
    const RegionTask &task,
    int thread_id,
    ProducerConsumerQueue<ThreadStats> &result_queue)
{
    samFile *in = sam_open(task.bam_path.c_str(), "r");
    if (in == NULL)
    {
        std::cerr << "Thread " << thread_id << ": Failed to open BAM file" << std::endl;
        return;
    }

    hts_set_threads(in, 1);

    bam_hdr_t *header = sam_hdr_read(in);
    if (header == NULL)
    {
        std::cerr << "Thread " << thread_id << ": Failed to read BAM header" << std::endl;
        sam_close(in);
        return;
    }

    hts_idx_t *idx = sam_index_load(in, task.bam_path.c_str());
    if (idx == NULL)
    {
        std::cerr << "Thread " << thread_id << ": Failed to load BAM index. Make sure .bai file exists." << std::endl;
        bam_hdr_destroy(header);
        sam_close(in);
        return;
    }

    ThreadStats local_stats;
    bam1_t *b = bam_init1();

    for (const std::string &region : task.regions)
    {
        hts_itr_t *iter = sam_itr_querys(idx, header, region.c_str());
        if (iter == NULL)
        {
            continue;
        }

        while (sam_itr_next(in, iter, b) >= 0)
        {
            process_single_read(b, task.mapq_cut, header, task.bed_regions, local_stats);
        }

        hts_itr_destroy(iter);
    }

    if (task.process_unmapped)
    {
        hts_itr_t *iter = sam_itr_querys(idx, header, "*");
        if (iter != NULL)
        {
            while (sam_itr_next(in, iter, b) >= 0)
            {
                process_single_read(b, task.mapq_cut, header, task.bed_regions, local_stats);
            }
            hts_itr_destroy(iter);
        }
    }

    bool sent = false;
    while (!sent)
    {
        sent = result_queue.write(local_stats);
        if (!sent)
        {
            std::this_thread::yield();
        }
    }

    bam_destroy1(b);
    hts_idx_destroy(idx);
    bam_hdr_destroy(header);
    sam_close(in);
}

inline void write_stats(std::ostream &out, const BamStats &stats, int mapq_cut, const Bed *bed,
                        const DepthStats *depth_stats = nullptr,
                        const std::vector<uint32_t> *depth_thresholds = nullptr)
{
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

    if (bed)
    {
        out << "Target region size(bp): " << bed->get_length() << "\n";
        out << "Ontarget reads: " << stats.on_target << "\n";

        if (depth_stats && depth_thresholds)
        {
            out << std::fixed << std::setprecision(2);
            out << "Coverage: " << depth_stats->get_coverage() << "%\n";
            out << "Mean depth: " << depth_stats->get_mean_depth() << "x\n";
            out << "Mean depth (rmdup): " << depth_stats->get_mean_depth_rmdup() << "x\n";
            out << "Uniformity (bases >= 0.2 * mean): " << depth_stats->get_uniformity(0.2) << "%\n";
            for (uint32_t threshold : *depth_thresholds)
            {
                double pct = depth_stats->get_coverage_at_depth(threshold);
                out << "Coverage (>= " << threshold << "x): " << pct << "%\n";
            }
        }
    }
}


#endif
