# Bulk RNASeq Data Downloader

This script allows you to efficiently download RNASeq SRA (Sequence Read Archive) files from NCBI and ENA (European Nucleotide Archive). It simplifies the process of bulk downloading SRA files associated with a list of accession numbers.

## Prerequisites

Before using this script, make sure you have the following prerequisites set up:

1. **Accession List**: Prepare a text file containing the list of SRA accession numbers you want to download. The file should contain one accession number per line. You can obtain this list from the SRA website.

2. **Download Tool**: Choose between two download tools:
   - `wget` (default): A command-line tool for downloading files from the web.
   - `aspera`: Aspera Connect, a high-speed data transfer tool for faster downloads (requires setup, see below).

3. **Aspera (Optional)**: If you choose to use Aspera for downloading, you need to have the Aspera SSH key (private key) ready. You'll specify the path to this key when using the script.

## Usage

Run the script with the following options:

```bash
./download_script.sh [-a <accession_list>] [-o <output_dir>] [-t <download_tool>] [-k <aspera_key_path>]
