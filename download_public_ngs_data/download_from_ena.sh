#!/bin/bash

# Initialize variables with default values
accession_list="SraAccList.txt"
output_dir="output"
download_tool="wget"  # Default download tool
aspera_key_path=""    # Path to Aspera SSH key (if using Aspera)

# Function to display script usage
usage() {
    echo "Usage: $0 [-a <accession_list>] [-o <output_dir>] [-t <download_tool>] [-k <aspera_key_path>]"
    echo "  -a <accession_list>: Specify the accession list file (default: SraAccList.txt)"
    echo "  -o <output_dir>: Specify the output directory (default: output)"
    echo "  -t <download_tool>: Specify the download tool (wget or aspera, default: wget)"
    echo "  -k <aspera_key_path>: Specify the path to the Aspera SSH key (required for aspera)"
    exit 1
}

# Parse command-line options
while getopts "a:o:t:k:" opt; do
    case "$opt" in
        a) accession_list="$OPTARG";;
        o) output_dir="$OPTARG";;
        t) download_tool="$OPTARG";;
        k) aspera_key_path="$OPTARG";;
        \?) usage;;
    esac
done

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Download SRA files based on the selected download tool
if [ "$download_tool" = "wget" ]; then
    while IFS= read -r accession; do
        echo "Downloading $accession..."
        $download_tool "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$accession/*" -P "$output_dir" &
        wait
    done < "$accession_list"
elif [ "$download_tool" = "aspera" ]; then
    if [ -z "$aspera_key_path" ]; then
        echo "Error: Specify the path to the Aspera SSH key (-k) when using Aspera."
        exit 1
    fi
    while IFS= read -r accession; do
        echo "Downloading $accession using Aspera..."
        ascp -QT -l 1000m -P33001 -i "$aspera_key_path" "era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/$accession" "$output_dir" &
        wait
    done < "$accession_list"
else
    echo "Error: Invalid download tool specified. Use 'wget' or 'aspera'."
    exit 1
fi

echo "Download completed."
