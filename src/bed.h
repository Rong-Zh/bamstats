#ifndef BED_H
#define BED_H

#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>

struct Interval
{
    uint64_t start;
    uint64_t end;

    bool operator<(const Interval &other) const
    {
        return start < other.start;
    }

    uint64_t length() const { return end - start; }

    bool overlaps(const Interval &other) const
    {
        return start < other.end && other.start < end;
    }
    
    // Check if a position range overlaps with this interval
    bool overlaps(uint64_t pos_start, uint64_t pos_end) const
    {
        return pos_start < end && pos_end > start;
    }
};

class Bed
{
public:
    Bed(const std::string &bed_file, bool sorted = false);
    ~Bed();
    
    void sort(size_t max_threads = 4);
    std::unordered_map<std::string, std::vector<Interval>> merge(size_t max_threads = 4);
    uint64_t get_length() const;
    
    // Check if a read position overlaps with any BED interval on the given chromosome
    bool is_ontarget(const char* chrom, uint64_t read_start, uint64_t read_end) const;
    
    // Get the intervals map (const reference)
    const std::unordered_map<std::string, std::vector<Interval>>& get_intervals() const {
        return intervals_;
    }
    
    bool is_empty() const {
        return intervals_.empty();
    }

private:
    bool sorted_;
    uint64_t length_;
    std::unordered_map<std::string, std::vector<Interval>> intervals_;
};

#endif // BED_H
