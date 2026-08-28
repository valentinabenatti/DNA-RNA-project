# DNA Methylation Array Analysis

Bioinformatics analysis pipeline for DNA methylation data (Illumina Infinium Methylation Array) implemented in R / R Markdown.

## Overview
This repository hosts the final project developed by Group 3 for the DNA/RNA Dynamics course (MSc in Bioinformatics, University of Bologna, a.y. 2025/2026). It provides an R pipeline for Illumina HumanMethylation450K array data analysis. The workflow covers quality control, normalization and preprocessing via preprocessQuantile, principal component analysis (PCA), and statistical testing using dmpFinder to detect differentially methylated positions between control (CTRL) and disease (DIS) samples.


## Files
- `progetto.Rmd`: R Markdown notebook containing the full analysis and commented code.
- `progetto.html`: Rendered HTML report with interactive plots and outputs.
- `progetto.R`: Standalone R script for the core pipeline steps.

## Requirements

This project requires R (>= 4.0) and the following packages:
```r
install.packages("BiocManager")
install.packages("factoextra")
install.packages("qqman")
install.packages("gplots")

BiocManager::install("minfi")
BiocManager::install("IlluminaHumanMethylation450kmanifest")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
