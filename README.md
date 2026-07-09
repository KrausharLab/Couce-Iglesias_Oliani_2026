# Proteomics and multi-omics analysis scripts

This repository contains the R scripts used to generate the figures and processed objects in the paper "Recording proteome synthesis and turnover in the embryo at hour-timescale, spatial, and single-cell resolution with MEMBRYO", Couce-Iglesias, M., Oliani, E. et al, 2026.


## Repository structure

```text
scripts/            Scripts for reuse
README.md           This file
```

## Script order

Run the scripts in this order, because later scripts read objects exported by earlier scripts.

1. `scripts/01_Fig1b_FreeAA.R`
2. `scripts/02a_PreprocessingSILAC&controlSamples.R`
3. `scripts/02b_Fig1_QCproteomics.R`
4. `scripts/03a_IntegrationMultiomics.R`
5. `scripts/03b_Fig2_RiboseqProteomics2.R`
6. `scripts/03c_Fig3_MultiomicsAnalysis.R`
7. `scripts/04_Fig4_SpatialProteomics.R`

Edit the **User settings** block at the top of each script before running. The current paths are the original project paths and are kept to avoid silently changing the workflow.

## Required R packages

```r
install.packages(c(
  "dplyr", "tidyr", "stringr", "ggplot2", "ggpubr", "openxlsx",
  "readxl", "tibble", "purrr", "rstatix", "viridis", "RColorBrewer",
  "corrplot", "ggrepel", "ggpointdensity", "gridExtra", "scales",
  "fitdistrplus", "nortest", "MASS", "TOSTER", "protr", "sf"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "ComplexHeatmap", "circlize", "clusterProfiler", "org.Mm.eg.db", "proDA"
))
```

Exact package versions from the original environment were not available. For a paper repository, add a `renv.lock` file from the machine where the scripts are known to run.



## Reproducibility notes: These are cleaned analysis scripts, not a fully portable pipeline. For full reproducibility, add an input-file manifest, a `renv.lock` file, relative paths or a single `config.R`, and `sessionInfo()` from the successful run.
