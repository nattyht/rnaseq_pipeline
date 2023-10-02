# FastQC and MultiQC Analysis Script

This Bash script simplifies the process of analyzing high-throughput sequence data using FastQC and MultiQC tools. It performs quality control analysis on fastq.gz files and generates summary reports.

## Prerequisites

Before using this script, make sure you have the following prerequisites set up:

1. **FastQC**: Install FastQC on your system. FastQC is a high throughput sequence quality control analysis tool. You can download and install it from the [FastQC website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

2. **MultiQC**: Install MultiQC on your system. MultiQC is a tool to create aggregate reports from multiple bioinformatics analyses. You can install it using `pip`:

```
pip install multiqc
```

## Input Data
Ensure you have the input data in the form of *.fastq.gz files in the current directory.

## Usage
Clone this repository or download the script.

Run the script using the following command:


```
./fastqc_multiqc_script.sh [-o <output_dir>] [-m <multiqc_options>]
```

-o <output_dir> (optional): Specify the directory for FastQC results (default: FastQCResults).

-m <multiqc_options> (optional): Additional options to pass to MultiQC. For example, you can specify -m "--cl_config general:read_count_file", where --cl_config is an option supported by MultiQC. Refer to the MultiQC documentation for more customization options.

FastQC will be run on all *.fastq.gz files in the current directory, and the results will be saved in the specified output directory.

MultiQC will generate a summary report for the FastQC results and save it in the current directory.

### Example Usage

Run the script with default options:

```
./fastqc_multiqc_script.sh
```

Run the script with custom output directory and MultiQC options:

```
./fastqc_multiqc_script.sh -o MyFastQCResults -m "--cl_config general:read_count_file"
```

## Additional Information
You can customize the FastQC and MultiQC analyses further by modifying the script or using additional FastQC and MultiQC options as needed.

If you encounter any issues or have questions about the analysis, please refer to the FastQC documentation and MultiQC documentation for guidance.
