
# Differential Gene Expression Analysis using DESeq2

This repository contains an end-to-end RNA-seq differential gene expression
analysis performed in R using the **DESeq2** package.
The workflow demonstrates how raw gene count data can be processed, analyzed,
and visualized to identify genes that are differentially expressed between
**treated** and **untreated** samples.

---

## 📌 Overview of the Workflow

The analysis includes the following steps:

- Loading raw RNA-seq count data and sample metadata
- Filtering low-count genes
- Normalization and dispersion estimation using DESeq2
- Differential expression testing
- Visualization of results using:
  - MA plot
  - Volcano plot with labeled top genes
  - Heatmap of the top differentially expressed genes

---

## 📂 Files in This Repository

- **deseq2_analysis.R**  
  Main R script containing the full DESeq2 workflow.

- **counts.csv**  
  Gene-level raw count matrix (rows = genes, columns = samples).

- **metadata2.csv**  
  Sample metadata including treatment information.

- **outcome.csv**  
  Full DESeq2 differential expression results.

- **best_genes.csv**  
  Top 10 differentially expressed genes ranked by adjusted p-value.

- **volcano_plot_labeled.png**  
  Volcano plot highlighting significant genes.

- **heatmap_top30.png**  
  Heatmap of the top 30 differentially expressed genes (VST-transformed).

---

## 🧪 Methods Summary

- Differential expression analysis was performed using **DESeq2**
- Low-count genes were filtered prior to analysis
- Multiple testing correction was applied using the Benjamini–Hochberg method
- Variance-stabilizing transformation (VST) was used for heatmap visualization

---

## 🔧 Requirements

R packages used in this analysis:

- DESeq2
- dplyr
- ggplot2
- pheatmap

---

## ▶️ How to Run the Analysis

1. Place all input files (`counts.csv`, `metadata2.csv`) in the same directory
2. Open R or RStudio
3. Run the script:

```r
source("deseq2_analysis.R")
