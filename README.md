# COVID-19 Lung Biopsy RNA-seq Analysis

January 2024

This repository contains the analysis of differential gene expression in COVID-19 infected versus uninfected lung biopsy samples, using data from Blanco et al. 2020 (Cell).

## Project Overview
The analysis pipeline examines transcriptional changes in lung tissue following SARS-CoV-2 infection through:
1. Alignment-free gene quantification using kallisto
2. Differential expression analysis using DESeq2
3. Visualization of expression patterns
4. Interpretation of key regulated genes

## Data Sources 
- RNA-seq data from GEO accession GSE147507
- Uninfected samples: SRR11517725, SRR11517729
- Infected samples: SRR11517733-40 (pooled)

## Key Findings
- Identified significant up-regulation of cytokine genes (CCL4L2, CCL3)
- Found down-regulation of T-cell regulatory genes (CD83, CD274)
- Results suggest targeted immune response to viral infection

## Dependencies
- kallisto v0.46.1
- R v4.0+
- DESeq2
- ggplot2
- gplots
- viridis

## Usage
Detailed instructions for running the analysis pipeline are available in `code/README.md`.
