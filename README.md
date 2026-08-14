# Jaccard Distance-based t-SNE Clustering for Molecular Classification

This repository contains the R code used for dimensionality reduction and clustering analysis of molecular data based on Jaccard distance and t-SNE, as described in our study on MRG-mutated AML.

## Overview

This pipeline performs the following analyses:

- Calculates Jaccard distance matrices between samples based on binary mutation data
- Performs t-SNE dimensionality reduction using Jaccard distance as the input
- Identifies clusters using DBSCAN
- Visualizes clustering results across different perplexity values
- Provides elbow method for optimal DBSCAN parameter selection

## Requirements

The following R packages are required:

| Package | Version | Purpose |
|---------|---------|---------|
| `readxl` | | Read Excel files |
| `openxlsx` | | Read/write Excel files |
| `stringr` | | String manipulation |
| `dplyr` | | Data manipulation |
| `proxy` | | Jaccard distance calculation |
| `Rtsne` | | t-SNE dimensionality reduction |
| `dbscan` | | DBSCAN clustering |
| `ggplot2` | | Data visualization |

### Installation

```r
install.packages(c("readxl", "openxlsx", "stringr", "dplyr", "proxy", "Rtsne", "dbscan", "ggplot2"))