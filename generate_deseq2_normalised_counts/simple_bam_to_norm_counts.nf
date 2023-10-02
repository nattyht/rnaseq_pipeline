#!/usr/bin/env nextflow

// Define input parameters
params {
    ensemblRelease = 'release-93' // Default Ensembl release
    bamFilesDir = '.' // Default directory containing BAM files
    outputDir = '.' // Default output directory
    coldata = './coldata.csv' // Path to coldata.csv
    designFormula = '~ condition' // Default design formula
}

// Define the Ensembl GTF URL
def ensemblGtfUrl = "ftp://ftp.ensembl.org/pub/${params.ensemblRelease}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${params.ensemblRelease}.gtf.gz"

// Define the workflow
workflow {
    
    // Download and unzip the Ensembl GTF file
    input:
    file ensemblGtf from ensemblGtfUrl

    script:
    """
    wget $ensemblGtfUrl
    gunzip Homo_sapiens.GRCh38.${params.ensemblRelease}.gtf.gz
    """
    
    // List all BAM files in the specified directory
    script:
    """
    bamFiles=$(find ${params.bamFilesDir} -name "*.bam" -type f -printf '%p ')
    """
    
    // Create the output directory if it doesn't exist
    script:
    """
    mkdir -p ${params.outputDir}
    """
    
    // Run featureCounts
    input:
    file gtfFile from 'Homo_sapiens.GRCh38.${params.ensemblRelease}.gtf'
    
    output:
    file("${params.outputDir}/counts.txt") into countsTxtFile

    script:
    """
    featureCounts -T 4 -s 2 -p -t exon -g gene_id -a $gtfFile -o ${params.outputDir}/counts.txt $bamFiles
    """
    
    // Import to R and perform DESeq2 analysis
    script:
    """
    library(DESeq2)
    library(stringr)
    library(ggplot2)
    
    counts <- read.csv("${params.outputDir}/counts.txt", sep="", header=T, skip=1, row.names = "Geneid")
    counts <- counts[, -c(1:5)]
    
    colData <- read.csv("${params.coldata}")
    colData <- read.csv("${params.coldata}", row.names = 1, stringsAsFactors = FALSE)
    
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ${params.designFormula})
    dds <- DESeq(dds)
    
    res <- results(dds)
    
    nco <- counts(dds, normalized = TRUE)
    write.csv(nco, file="${params.outputDir}/norm_counts.csv")
    """
}

// Output the normalized counts
output:
file("${params.outputDir}/norm_counts.csv") into normCountsCsv
