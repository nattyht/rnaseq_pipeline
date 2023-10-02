# Alignment Scripts Readme

This repository contains two scripts for aligning RNA-Seq data using HISAT2 and STAR. The scripts allow you to define options for input and output directories, reference genome versions, URLs for downloading reference files, and more.

## HISAT2 Alignment Script

### Usage

To use the HISAT2 alignment script, run the following command:

```bash
./hisat2_alignment.sh [OPTIONS]
```

### Options
- -i, --input: Path to the input directory containing raw FASTQ files.
- -o, --output: Path to the output directory for HISAT2 aligned data.
- -g, --genome_version: Version of the reference genome (e.g., GRCh38.p13).
- -u, --genome_url: URL to download the reference genome FASTA file.

### Example
```
./hisat2_alignment.sh -i /path/to/input -o /path/to/output -g GRCh38.p13 -u http://example.com/genome.fa.gz
```

## STAR Alignment Script

### Usage

To use the STAR alignment script, run the following command:

```
./star_alignment.sh [OPTIONS]
```

### Options
- -i, --input: Path to the input directory containing raw FASTQ files.
- -o, --output: Path to the output directory for STAR aligned data.
- -g, --genome_version: Version of the reference genome (e.g., GRCh38.p13).
- -u, --genome_url: URL to download the reference genome FASTA file.
- -t, --gtf_version: Version of the GTF annotation file (e.g., gencode.v33.primary_assembly.annotation.gtf).
- -v, --gtf_url: URL to download the GTF annotation file.

### Example

```
./star_alignment.sh -i /path/to/input -o /path/to/output -g GRCh38.p13 -u http://example.com/genome.fa.gz -t gencode.v33.primary_assembly.annotation.gtf -v http://example.com/gtf_annotation.gtf
```

These scripts will download the required reference files based on the provided options and perform the alignment of RNA-Seq data accordingly. Make sure to adjust the options and paths as needed for your specific use case.
