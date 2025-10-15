#!/bin/bash
# Example usage script for BAM Statistics Tool

# Configuration
BAM_FILE="/disk1/Project/KCWES/devel/20250922-SJZP/mapping/SJZP-202515/SJZP-202515.sort.bam"
OUTPUT_DIR="./SJZP-202515"
OUTPUT_FILE="$OUTPUT_DIR/SJZP-202515.bam.stats.tsv"
THREADS=4
MAPQ=30

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check if BAM file exists
if [ ! -f "$BAM_FILE" ]; then
    echo "Error: BAM file not found: $BAM_FILE"
    exit 1
fi

# Check if BAM index exists
if [ ! -f "$BAM_FILE.bai" ]; then
    echo "Warning: BAM index not found. Creating index..."
    samtools index "$BAM_FILE"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to create BAM index"
        exit 1
    fi
fi

# Run the statistics tool
echo "Running BAM statistics analysis..."
echo "BAM file: $BAM_FILE"
echo "Output: $OUTPUT_FILE"
echo "Threads: $THREADS"
echo "MAPQ cutoff: $MAPQ"
echo ""

./bam_stats -bam "$BAM_FILE" -stats "$OUTPUT_FILE" -nt $THREADS -mapq $MAPQ

if [ $? -eq 0 ]; then
    echo ""
    echo "Analysis complete! Results saved to: $OUTPUT_FILE"
else
    echo ""
    echo "Error: Analysis failed"
    exit 1
fi
