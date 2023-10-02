# RNA-seq Analysis Script

## Overview

This R script is designed for the analysis of RNA-seq data, specifically for comparing gene expression between different conditions or groups. It utilizes the DESeq2 package for differential expression analysis and includes various visualization and enrichment analysis steps. The script assumes you have already obtained the necessary RNA-seq data and sample information.

## Prerequisites

Before using this script, ensure you have the following:

- R and RStudio installed on your computer.
- Required R packages installed. You can install these packages using the `install.packages()` function in R if they are not already installed. The required packages are listed within the script.

## Usage

1. **Set Up Directory Structure**:

   - Create a directory where you want to store your analysis results.
   - Organize your RNA-seq data files and sample information file within a directory structure similar to the one used in the script. Make sure to specify the paths to these files in the script.

2. **Configure Input and Output**:

   - Open the R script in RStudio.
   - Edit the following variables at the beginning of the script to match your data:
     - `input_dir`: Path to the directory containing your RNA-seq data and sample information.
     - `sample_info_file`: Filename of the sample information file.
     - `output_dir`: Path to the directory where you want to save the analysis results.

3. **Run the Script**:

   - Run the entire R script in RStudio. Make sure to save your work in RStudio before running the script.
   - The script will perform various analysis steps, including data normalization, differential expression analysis, visualization, and pathway enrichment analysis.

4. **Review Results**:

   - Check the specified `output_dir` for the generated output files, including differential expression results, heatmaps, and pathway enrichment results.

## Output

The script will generate several output files and visualizations in the specified `output_dir`. These may include:

- Normalized gene counts (in CSV format).
- Differential expression results (in CSV format).
- Volcano plots and MA plots visualizing differential expression.
- Heatmaps of significantly different genes.
- Pathway enrichment analysis results (in CSV format).
- Visualizations of enriched pathways and Gene Ontology terms.

## Troubleshooting

If you encounter any issues or errors while running the script, please refer to the R documentation or the documentation of the specific R packages used for assistance.

## License

This script is provided under an open-source license. You are free to use, modify, and distribute it in accordance with the terms of the license (if specified).

## Author

Dr. Nathaniel Hafford Tear
