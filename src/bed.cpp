
#include "bed.h"
#include <iostream>
#include <thread>
#include <vector>
#include <atomic>
#include <string>
#include <fstream>
#include <unordered_map>
#include <algorithm>

void split(const std::string &str, char delimiter, std::vector<std::string> &result)
{
    if (str.empty())
        return;
    const char *data = str.data();
    const char *start = data;
    const char *end = data + str.size();

    for (const char *p = data; p < end; ++p)
    {
        if (*p == delimiter)
        {
            if (p > start)
            {
                result.emplace_back(start, p - start);
            }
            start = p + 1;
        }
    }
    if (start < end)
    {
        result.emplace_back(start, end - start);
    }
    return;
}

bool startswith(const std::string &str, const std::string &prefix)
{
    if (prefix.length() > str.length())
    {
        return false;
    }
    return str.compare(0, prefix.length(), prefix) == 0;
}

Bed::Bed(const std::string &bed_file, bool sorted) : sorted_(sorted), length_(0)
{
        std::ifstream file(bed_file);
        if (!file.is_open())
        {
            std::cerr << "Failed to open file: " << bed_file << std::endl;
            return;
        }
        std::string line;
        std::vector<std::string> fields;
        fields.reserve(12);
        std::string header;
        
        // Read and skip header lines
        size_t nline = 0;
        while (std::getline(file, line))
        {
            if (startswith(line, "#") || startswith(line, "browser") || startswith(line, "track"))
            {
                header += line + "\n";
                nline++;
            }
            else
            {
                break;
            }
        }
        // Process the first data line and all remaining lines
        // Use do-while to process the line already read by previous while loop
        do
        {
            nline++;         
            // Skip empty lines or additional header lines in the middle
            if (line.empty() || startswith(line, "#") || startswith(line, "browser") || startswith(line, "track"))
            {
                header += line + "\n";
                continue;
            }
            
            split(line, '\t', fields);
            if (fields.size() < 3)
            {
                std::cerr << "Invalid BED line (less than 3 fields) at line " << nline << ": " << line << std::endl;
                std::exit(EXIT_FAILURE);
            }
            
            const std::string &chrom = fields[0];
            uint64_t start = std::stoul(fields[1]);
            uint64_t end = std::stoul(fields[2]);
            
            if (start > end)
            {
                std::cerr << "Invalid interval in BED file at line " << nline << ": " << line << std::endl;
                std::exit(EXIT_FAILURE);
            }
            
            intervals_[chrom].emplace_back(Interval{start, end});
            fields.clear();
            
        } while (std::getline(file, line));
    }

Bed::~Bed() 
{
}

void Bed::sort(size_t max_threads)
    {
        if (sorted_)
            return;
        std::vector<std::string> chroms;
        for (const auto &kv : intervals_)
            chroms.emplace_back(kv.first);

        size_t num = chroms.size();
        std::atomic<size_t> idx(0);
        size_t thread_count = std::min(max_threads, num);
        std::vector<std::thread> threads;

        auto worker = [&]()
        {
            while (true)
            {
                size_t i = idx.fetch_add(1);
                if (i >= num)
                    break;
                std::sort(intervals_.at(chroms[i]).begin(), intervals_.at(chroms[i]).end());
            }
        };

        for (size_t t = 0; t < thread_count; ++t)
            threads.emplace_back(worker);

        for (auto &th : threads)
            th.join();

        sorted_ = true;
    }

std::unordered_map<std::string, std::vector<Interval>> Bed::merge(size_t max_threads)
    {
        length_ = 0;
        std::vector<std::string> chroms;
        for (const auto &kv : intervals_)
            chroms.emplace_back(kv.first);

        size_t num = chroms.size();
        std::vector<std::vector<Interval>> merged_results(num);
        std::vector<uint64_t> lengths(num, 0);

        std::atomic<size_t> idx(0);
        size_t thread_count = std::min(max_threads, num);
        std::vector<std::thread> threads;

        auto worker = [&]()
        {
            while (true)
            {
                size_t i = idx.fetch_add(1);
                if (i >= num)
                    break;
                const auto &range = intervals_.at(chroms[i]);
                std::vector<Interval> sorted = range;
                if (!sorted_)
                    std::sort(sorted.begin(), sorted.end());
                std::vector<Interval> merged;
                uint64_t length = 0;
                merged.reserve(sorted.size());
                auto current = sorted.front();
                for (size_t j = 1; j < sorted.size(); ++j)
                {
                    if (sorted[j].start <= current.end)
                        current.end = std::max(current.end, sorted[j].end);
                    else
                    {
                        merged.emplace_back(current);
                        length += current.length();
                        current = sorted[j];
                    }
                }
                merged.emplace_back(current);
                length += current.length();
                merged_results[i] = std::move(merged);
                lengths[i] = length;  // Fixed: should be assignment, not addition
            }
        };

        for (size_t t = 0; t < thread_count; ++t)
            threads.emplace_back(worker);

        for (auto &th : threads)
            th.join();

        std::unordered_map<std::string, std::vector<Interval>> merged_map;
        for (size_t i = 0; i < num; ++i)
        {
            merged_map[chroms[i]] = std::move(merged_results[i]);
            length_ += lengths[i];
        }

        // Update internal intervals with merged results
        intervals_ = std::move(merged_map);

        return intervals_;
    }

uint64_t Bed::get_length() const
{
    return length_;
}

bool Bed::is_ontarget(const char* chrom, uint64_t read_start, uint64_t read_end) const
{
    if (intervals_.empty()) {
        return true;  // No BED file provided, all reads are "on target"
    }
    
    auto it = intervals_.find(chrom);
    if (it == intervals_.end()) {
        return false;  // Chromosome not in BED file
    }
    
    // Check if read overlaps with any interval
    for (const auto& interval : it->second) {
        if (interval.overlaps(read_start, read_end)) {
            return true;
        }
    }
    
    return false;
}
