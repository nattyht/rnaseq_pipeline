#!/bin/bash

# Initialize variables with default values
output_dir="./fastqcresults"
fastq_dir="."

# Function to display script usage
usage() {
    echo "Usage: $0 [-o <output_dir>] [-f <fastq_dir>]"
    echo "  -o <output_dir>: Specify the directory for FastQC results to be stored in (default: ./fastqcresults)"
    echo "  -f <fastq_dir>: Specify the directory containing *.fastq.gz files (default: current directory)"
    exit 1
}

# Parse command-line options
while getopts "o:f:" opt; do
    case "$opt" in
        o) output_dir="$OPTARG";;
        f) fastq_dir="$OPTARG";;
        \?) usage;;
    esac
done

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Run FastQC with 6 threads on all *.fastq.gz files in the specified directory
fastqc -t 6 "$fastq_dir"/*.fastq.gz -o "$output_dir"

# Run multiqc in the current directory
multiqc .
