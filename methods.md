# Methods Documentation

## RNA-seq Analysis Pipeline

### 1. Transcript Quantification
Kallisto v0.46.1 was used for alignment-free transcript quantification:
- Created index from human RefSeq transcripts (hg38)
- Quantified single-end reads using parameters:
  - Fragment length: 200bp
  - SD: 50bp
  - Bootstrap samples: 30

### 2. Differential Expression Analysis
DESeq2 analysis workflow:
- Generated count matrix from kallisto output
- Defined experimental design:
  - 2 uninfected samples
  - 1 infected sample (8 technical replicates pooled)
- Performed differential expression testing
- Applied multiple testing correction (Benjamini-Hochberg)
- Significance threshold: adjusted p-value < 0.05

### 3. Visualization
Generated visualizations using R:
- Volcano plot showing log2 fold changes vs significance
- Expression heatmap of significantly regulated genes
- Log2 fold change heatmap for DEGs
- All plots use viridis color palette
- Highlighted IFI6 gene as reference

### 4. Gene Function Analysis
- Identified top differentially expressed genes
- Analyzed gene functions using UCSC Genome Browser
- Interpreted biological relevance to COVID-19 response

## Data Processing Notes
- Raw reads processed without downsampling
- Technical replicates merged for infected sample
- TPM normalization used for visualization
- Integer rounding applied to kallisto estimates for DESeq2
