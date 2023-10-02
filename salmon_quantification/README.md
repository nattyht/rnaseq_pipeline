# Salmon Quantification Script

This Bash script automates the process of Salmon quantification for RNA-Seq data analysis. It allows you to specify the location of FASTQ files, the Ensembl release version, and whether the data is single-end or paired-end. Salmon quantification is performed using the Gencode transcriptome reference.

## Prerequisites

Before using this script, ensure that you have the following prerequisites installed:

- [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html)
- [wget](https://www.gnu.org/software/wget/)

## Usage

1. Clone this repository or download the `salmon_quant.sh` script.

2. Make the script executable by running the following command:

   ```bash
   chmod +x salmon_quant.sh

### Run the script with the following parameters:

- <FASTQ_DIR>: Directory containing FASTQ files (either single-end or paired-end).
- <ENSEMBL_RELEASE>: Ensembl release version (e.g., 33).
- <SINGLE_OR_PAIRED>: 'single' for single-end data, 'paired' for paired-end data.

Example usage for single-end data:
```
./salmon_quant.sh /path/to/fastq/files 33 single
```

Example usage for paired-end data:
```
./salmon_quant.sh /path/to/fastq/files 33 paired
```

The script will create a 'quants' directory to store the Salmon quantification results if it doesn't exist. It will then download the Gencode transcriptome reference, index it, and perform Salmon quantification for each sample in the specified directory. The results will be stored in the 'quants' directory.

## Output
The script generates Salmon quantification results in the 'quants' directory, with one subdirectory for each sample. The results include transcript-level abundance estimates.
