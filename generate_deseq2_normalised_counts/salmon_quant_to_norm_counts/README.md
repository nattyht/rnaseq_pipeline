# Salmon Quants To DESeq2 Normalised Counts

This Rmarkdown document provides a workflow for analyzing Salmon quants data using R, with a focus on generating DESeq normalized counts at both the transcript and gene levels. This document allows the user to specify input and output file paths and names.

## Prerequisites

Before using this Rmarkdown file, ensure that you have the following prerequisites:

1. **R**: Make sure you have R installed on your system. You can download it from [here](https://www.r-project.org/).

2. **R Packages**: Install the required R packages by running the following commands in your R environment:

   ```
   install.packages("tximeta")
   install.packages("DESeq2")
## Getting Started

1. Clone the Repository: Clone or download this repository to your local machine.

2. Modify User Input: Open the SalmonQuantsAnalysis.Rmd file in an Rmarkdown editor or RStudio. In the "User Input" section, specify the following information:

  - dir: Path to the folder containing your Salmon quants data.
  - sample_info_file: The name of your sample information CSV file.
  - transcript_output_file: Desired output file name for transcript-level counts.
  - gene_output_file: Desired output file name for gene-level counts.
3. Knit the Document: Knit the Rmarkdown document (SalmonQuantsAnalysis.Rmd) to execute the analysis. This will generate the output files with normalized counts.

## Output Files

- norm_counts.transcript.csv: This file contains normalized counts at the transcript level.
- norm_counts.csv: This file contains normalized counts at the gene level.

## Notes
- Ensure that your Salmon quants data is organized in a folder structure compatible with the provided code.
- Modify the design formula in the DESeq analysis section if your experimental design differs from the example provided.

This document is meant as a template, and you should adapt it to your specific analysis needs.
