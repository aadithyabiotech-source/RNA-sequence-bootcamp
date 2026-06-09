# RNA-Seq & Microarray Bootcamp Workflows

This repository contains reproducible pipelines for:
- Microarray analysis (GSE15852, limma)
- RNA-Seq analysis (airway dataset, DESeq2)

## Contents
- `scripts/` → R scripts for limma and DESeq2
- `results/` → Plots (volcano, heatmap, PCA)
- `notebooks/` → RMarkdown notebooks for reproducibility

## Workflow Summary
1. Data retrieval (GEOquery / airway dataset)
2. Normalization (log2 transform / DESeq2 vst)
3. Differential expression (limma / DESeq2)
4. Visualization (volcano plot, heatmap, PCA)
5. Functional analysis (GO, KEGG planned)

## How to Run
```r
source("scripts/microarray_limma.R")
source("scripts/rnaseq_deseq2.R")
