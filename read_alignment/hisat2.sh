#!/bin/bash

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
    -i|--input)
        input_dir="$2"
        shift # past argument
        shift # past value
        ;;
    -o|--output)
        hisat2_output_dir="$2"
        shift # past argument
        shift # past value
        ;;
    -g|--genome_version)
        genome_version="$2"
        shift # past argument
        shift # past value
        ;;
    -u|--genome_url)
        genome_url="$2"
        shift # past argument
        shift # past value
        ;;
    *)
        echo "Unknown option: $1"
        exit 1
        ;;
    esac
done

# Download the reference genome if not already downloaded
if [ ! -e "${genome_version}.genome.fa" ]; then
    wget -O "${genome_version}.genome.fa.gz" "${genome_url}"
    gunzip "${genome_version}.genome.fa.gz"
fi

# Create HISAT2 index
hisat2-build "${genome_version}.genome.fa" "${genome_version}_index"

# Change directory to the input directory
cd "${input_dir}"

# Loop to align with HISAT2 and sort the resulting BAM file
for f in `ls *.fq.gz | sed 's/_[12].fq.gz//g' | sort -u`; do
    echo "Aligning ${f}..."
    hisat2 -x "${genome_version}_index" -p 6 -1 "${f}_1.fq.gz" -2 "${f}_2.fq.gz" | samtools sort -o "${hisat2_output_dir}/${f}.hisat.sorted.bam"
done

# Move aligned data to the desired location
echo "Alignment completed. Aligned data is in ${hisat2_output_dir}."

