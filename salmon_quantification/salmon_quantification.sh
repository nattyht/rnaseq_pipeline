#!/bin/bash

# Check for the required Salmon executable
if ! command -v salmon &> /dev/null; then
    echo "Salmon is not installed or not in the system's PATH. Please install Salmon before running this script."
    exit 1
fi

# Check for user input
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <FASTQ_DIR> <ENSEMBL_RELEASE> <SINGLE_OR_PAIRED>"
    echo "  <FASTQ_DIR>: Directory containing FASTQ files (either single-end or paired-end)"
    echo "  <ENSEMBL_RELEASE>: Ensembl release version (e.g., 33)"
    echo "  <SINGLE_OR_PAIRED>: 'single' for single-end data, 'paired' for paired-end data"
    exit 1
fi

FASTQ_DIR="$1"
ENSEMBL_RELEASE="$2"
SINGLE_OR_PAIRED="$3"

# Downloading gencode transcriptome fasta
wget -nc ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${ENSEMBL_RELEASE}/gencode.v${ENSEMBL_RELEASE}.transcripts.fa.gz

# Indexing transcriptome reference
if [ ! -d "gencode.v${ENSEMBL_RELEASE}.transcripts.index" ]; then
    salmon index --gencode -t gencode.v${ENSEMBL_RELEASE}.transcripts.fa.gz -i gencode.v${ENSEMBL_RELEASE}.transcripts.index
fi

# Create the 'quants' directory if it doesn't exist
mkdir -p quants

# Perform Salmon quantification based on data type (single-end or paired-end)
if [ "$SINGLE_OR_PAIRED" == "single" ]; then
    for f in $FASTQ_DIR/*.fastq.gz; do
        SAMPLE_NAME=$(basename "$f" .fastq.gz)
        echo "Processing sample $SAMPLE_NAME"
        salmon quant -i gencode.v${ENSEMBL_RELEASE}.transcripts.index -l A -r "$f" -p 8 --validateMappings -o "quants/${SAMPLE_NAME}_quant"
    done
elif [ "$SINGLE_OR_PAIRED" == "paired" ]; then
    for f1 in $FASTQ_DIR/*_1.fastq.gz; do
        SAMPLE_NAME=$(basename "$f1" _1.fastq.gz)
        echo "Processing sample $SAMPLE_NAME"
        f2="${f1/_1.fastq.gz/_2.fastq.gz}"
        salmon quant -i gencode.v${ENSEMBL_RELEASE}.transcripts.index -l A -1 "$f1" -2 "$f2" -p 8 --validateMappings -o "quants/${SAMPLE_NAME}_quant"
    done
else
    echo "Invalid data type. Please specify 'single' or 'paired' for <SINGLE_OR_PAIRED>."
    exit 1
fi

echo "Salmon quantification completed."
