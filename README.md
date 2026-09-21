# Developmental SILAC Proteomics Analysis

This repository contains the R scripts used to preprocess, quality-control, integrate, analyse, and visualise free-amino-acid mass spectrometry, bulk SILAC proteomics, RNA-seq, Ribo-seq, spatial SILAC proteomics, and single-cell SILAC data from the developing mouse brain.

The scripts are numbered according to their intended execution order and the main figures they support. Unless stated otherwise, they read data from `Data/`, write processed objects to `Output/Objects/`, figures to `Output/Figures/`, and timestamped logs to `Output/Logs/`.

The folder Data can be downloaded in the following link: https://nc.molgen.mpg.de/cloud/index.php/s/sfSRxTNMkLQF6dB

## Repository structure

```text
.
├── Data/
│   ├── Free_AminoAcid_MS/
│   ├── SILAC_Timecourse_Bulk_MS/
│   ├── SILAC_Titration_MS/
│   ├── SingleCell_DiBella2021/
│   ├── SILAC_SpatialProteomics_MS/
│   └── SILAC_SingleCellProteomics_MS/
├── Output/
│   ├── Figures/
│   ├── Logs/
│   └── Objects/
├── Scripts/
│   ├── 01_Fig1b_FreeAA.R
│   ├── 02a_PreprocessingSILAC&controlSamples.R
│   ├── 02b_Fig1_QCproteomics.R
│   ├── 02c_Subset_Proteins_Kinetics.R
│   ├── 03a_IntegrationMultiomics.R
│   ├── 03b_Fig2_RiboseqProteomics.R
│   ├── 03c_Fig3_MultiomicsAnalysis.R
│   ├── 03d_Fig3_CellTypeBias.R
│   ├── 04a_Fig4_SpatialProteomics_neocortex.R
│   ├── 04b_Fig4a_SpatialProteomics_hemisphere.R
│   └── 05_Fig5_scSILAC.R
```

Preserve the original filenames and table structures shown below, or update the corresponding **User settings** section in each script.


## Analysis workflow

### 1. Free-amino-acid analysis

#### `01_Fig1b_FreeAA.R`

Processes free-amino-acid MS peak tables and validates the measured transitions.

Main operations:

- standardises sample and metabolite annotations;
- joins sample metadata and BCA measurements;
- checks transition-area ratios against an amino-acid standard;
- validates retention-time agreement among transitions;
- selects one quantitative transition per amino acid;
- normalises signals by valine and total protein concentration;
- removes manually identified poor technical replicates;
- calculates heavy- and light-amino-acid percentages;
- interpolates heavy-labelling percentages at intermediate embryonic ages.

Required inputs in `Data/Free_AminoAcid_MS/`:

```text
Free_AA_results.xlsx
BCA.xlsx
Validation_Standard.xlsx
```

Principal output:

```text
Output/Objects/normValine_results_firstTransition_filterTechReplicates.rds
```

The script also generates Figure 1 and Extended Data Figure 1 PDF files.

---

### 2. Bulk SILAC preprocessing and quality control

#### `02a_PreprocessingSILAC&controlSamples.R`

Builds the processed bulk SILAC time-course, control, and titration datasets.

Main operations:

- attaches experimental metadata;
- calculates H/T and L/T ratios;
- removes manually selected poor-quality samples;
- normalises precursor-level heavy incorporation using age-specific free-amino-acid availability;
- filters proteins by replicate coverage;
- calculates iBAQ and riBAQ values;
- processes unlabelled control samples;
- combines three MaxQuant SILAC titration experiments.

Required inputs in `Data/SILAC_Timecourse_Bulk_MS/`:

```text
SILACtimecourse_ratios.csv
SILACtimecourse_light.csv
SILACtimecourse_pulse.csv
SILACtimecourse_total.csv
SILACtimecourse_precursors.csv
SILACtimecourse_samples_metadata.xlsx
Controls_light.csv
Controls_pulse.csv
evidence_ArgProConversion.txt
evidence_labellingEfficiency.txt
```

Required inputs in `Data/SILAC_Titration_MS/`:

```text
Exp173_proteinGroups.txt
Exp182_proteinGroups.txt
Exp195_proteinGroups.txt
```

Principal outputs:

```text
Output/Objects/ProteomicsTimecourse_SILAC.rds
Output/Objects/ProteomicsTimecourse_Control.rds
Output/Objects/ProteomicsTitration_SILAC.rds
```

#### `02b_Fig1_QCproteomics.R`

Generates quality-control plots for the bulk SILAC time course.

Main analyses include:

- SILAC titration intensity distributions;
- H/T and L/T incorporation across embryonic ages and incubation times;
- numbers of quantified proteins;
- PCA of log2 riBAQ values;
- Procrustes comparison before and after free-amino-acid normalisation;
- sample-level outlier detection using protein counts and median H/L ratios.

Required upstream objects:

```text
Output/Objects/ProteomicsTimecourse_SILAC.rds
Output/Objects/ProteomicsTitration_SILAC.rds
```

#### `02c_Subset_Proteins_Kinetics.R`

Prepares age-specific protein-intensity tables for downstream kinetic modelling.

The script:

- tests total protein abundance between 4 h and 8 h using `limma`;
- retains proteins that do not exceed the configured adjusted-P-value and log2-fold-change thresholds;
- requires heavy and light measurements at both time points;

Required upstream object:

```text
Output/Objects/ProteomicsTimecourse_SILAC.rds
```

Outputs:

```text
Output/Objects/Timecourse_Controls_t12_E12.5.xlsx
Output/Objects/Timecourse_Controls_t12_E13.5.xlsx
Output/Objects/Timecourse_Controls_t12_E14.5.xlsx
Output/Objects/Timecourse_Controls_t12_E15.5.xlsx
Output/Objects/Timecourse_Controls_t12_E16.5.xlsx
```

These workbooks are intended as inputs to the external half-life calculation step.

---

### 3. Multiomics integration and developmental analyses

#### `03a_IntegrationMultiomics.R`

Combines bulk SILAC proteomics, RNA-seq, Ribo-seq, and protein half-life estimates into a unified long-format dataset.

Main operations:

- reshapes all SILAC measurements into a common format;
- imports RNA-seq and Ribo-seq TPM data;
- imports age-specific protein half-life estimates;
- maps UniProt and Ensembl identifiers to mouse gene symbols with `org.Mm.eg.db`;
- generates replicate-level log2 values;
- computes developmental trajectory Z-scores;
- saves median Z-score and replicate-SD matrices for downstream clustering and plotting.

Additional required inputs in `Data/SILAC_Timecourse_Bulk_MS/`:

```text
RNAseq_Riboseq.xlsx
SILACtimecourse_t12_E12.5.xlsx
SILACtimecourse_t12_E13.5.xlsx
SILACtimecourse_t12_E14.5.xlsx
SILACtimecourse_t12_E15.5.xlsx
SILACtimecourse_t12_E16.5.xlsx
```

Outputs:

```text
Output/Objects/all_data.rds
Output/Objects/Zscore_RPF_RNAseq_riBAQH_T_HL_4h.rds
Output/Objects/Zscore_RPF_RNAseq_riBAQH_T_HL_4h_SD.rds
```

#### `03b_Fig2_RiboseqProteomics.R`

Compares Ribo-seq translation output with newly synthesised heavy proteome abundance.

Main analyses:

- stage-specific RPF-versus-riBAQ-H correlations;
- correlation of E15.5/E12.5 fold changes;
- identification of proteins outside ±2 residual standard deviations;
- GO enrichment of correlation outliers;
- k-means clustering of RPF and heavy-proteome developmental trajectories;
- GO enrichment of temporal clusters.

The number of temporal clusters is controlled by:

```r
k <- 10
```

A saved k-means object is reused when present to preserve the published clustering.

#### `03c_Fig3_MultiomicsAnalysis.R`

Performs integrated developmental clustering across RNA-seq, Ribo-seq, total proteome, heavy proteome, and H/L measurements.

Main analyses:

- developmental distribution plots;
- elbow-plot calculation;
- k-means clustering and heatmap generation;
- GO enrichment per multiomics cluster;
- cluster-level mean trajectories with standard deviations;
- trajectories and replicate variability for selected proteins.

The number of multiomics clusters is controlled by:

```r
k <- 11
```

#### `03d_Fig3_CellTypeBias.R`

Integrates published embryonic neocortex single-cell RNA-seq data and relates developmental cell-type composition to SILAC proteostasis measurements.

Main operations:

- reads five 10x HDF5 expression matrices;
- performs cell-level QC;
- integrates ages with Harmony;
- clusters cells and assigns marker-score-based cell types;
- calculates developmental cell-type proportions;
- plots SILAC trajectories for cell-type marker proteins;
- correlates cell-type proportions with SILAC metrics.

Required files in `Data/SingleCell_DiBella2021/`:

```text
GSM4635073_E12_5_filtered_gene_bc_matrices_h5.h5
GSM4635074_E13_5_filtered_gene_bc_matrices_h5.h5
GSM4635075_E14_5_filtered_gene_bc_matrices_h5.h5
GSM4635076_E15_5_S1_filtered_gene_bc_matrices_h5.h5
GSM4635077_E16_filtered_gene_bc_matrices_h5.h5
```

Principal output:

```text
Output/Objects/Fig3_CellTypeBias_Seurat.rds
```

---

### 4. Spatial proteomics

#### `04a_Fig4_SpatialProteomics_neocortex.R`

Processes microdissected neocortex and basal-ganglia spatial SILAC proteomics data.

Main operations:

- reshapes StackedLFQ ratio, total, light, and heavy output tables;
- attaches anatomical-region metadata;
- removes manually identified outlier samples;
- joins samples to GeoJSON spatial polygons;
- generates spatial protein-detection and intensity maps;
- performs proDA differential-enrichment analysis between ventricular zone and cortical plate;
- performs GO enrichment for differentially enriched proteins;
- displays selected protein spatial profiles and half-life positions.

Required files in `Data/SILAC_SpatialProteomics_MS/`:

```text
SILAC_Spatial_ratios.csv
SILAC_Spatial_total.csv
SILAC_Spatial_light.csv
SILAC_Spatial_pulse.csv
SILAC_Spatial_samples_metadata.xlsx
Geometry_coordinates.geojson
```

Principal model output:

```text
Output/Objects/fit_spatial_heavy_VZvsCP.rds
```

#### `04b_Fig4a_SpatialProteomics_hemisphere.R`

Processes the E14 full-slice DIA-NN spatial SILAC experiment.

Main operations:

- loads the DIA-NN report;
- annotates light, medium, and heavy SILAC peptide channels;
- filters contaminants and proteins supported by fewer than two precursors;
- constructs protein-level intensity matrices;
- calculates heavy/light ratios;
- maps detected-protein fractions onto the serial-section grid.

Required input:

```text
Data/SILAC_SpatialProteomics_MS/report_SpatialSILAC_Hemisphere.tsv
```

Principal output:

```text
Output/Objects/protMatrix_SpatialProteomics_hemisphere.rds
```


---

### 5. Single-cell SILAC

#### `05_Fig5_scSILAC.R`

Analyses single-cell SILAC measurements while retaining missing values and using NA-aware proDA distances.

Four assay representations are constructed:

| Assay | Source column | Interpretation |
|---|---|---|
| Total | `CS.Residual.Stacked.LandM.Intensity` | Total proteome |
| Light | `CS.Residual.L.Intensity` | L/H |
| Medium | `CS.Residual.M.Intensity` | M/H |
| Turnover | M/H − L/H | M/L turnover proxy |

Main operations:

- builds protein-by-cell matrices without imputation;
- filters proteins and cells by detection coverage;
- calculates NA-aware cell distances with `proDA::dist_approx`;
- performs multidimensional scaling, graph clustering, and UMAP;
- displays median total, M/L, and M/H values on the total-proteome UMAP for selected markers;
- compares total abundance with M/L and M/H values and residual boundaries for selected markers.

Required input in `Data/SILAC_SingleCellProteomics_MS/`:

```text
clean+TechNoiseSafety_unsafeRemoved_CSNorm.parquet
```

Principal outputs:

```text
Output/Objects/assay_matrices_all_assays.rds
Output/Objects/assay_object_total.rds
```

The main configurable parameters are grouped in the **User settings** section, including filtering thresholds, dimensionality, clustering resolution, and selected marker proteins.

## Recommended execution order

Run the scripts from the repository root:

```r
source("01_Fig1b_FreeAA.R")
source("02a_PreprocessingSILAC&controlSamples.R")
source("02b_Fig1_QCproteomics.R")
source("02c_Subset_Proteins_Kinetics.R")

# Run the external half-life calculation and place its output workbooks in
# Data/SILAC_Timecourse_Bulk_MS/ before continuing.

source("03a_IntegrationMultiomics.R")
source("03b_Fig2_RiboseqProteomics.R")
source("03c_Fig3_MultiomicsAnalysis.R")
source("03d_Fig3_CellTypeBias.R")
source("04a_Fig4_SpatialProteomics_neocortex.R")
source("04b_Fig4a_SpatialProteomics_hemisphere.R")
source("05_Fig5_scSILAC.R")
```

The spatial and single-cell branches are largely independent of each other, but `04a_Fig4_SpatialProteomics_neocortex.R` also reads `Output/Objects/all_data.rds` for the half-life density plot and therefore requires `03a_IntegrationMultiomics.R` to have completed.

## Software requirements

The scripts were written in R and use packages from CRAN, Bioconductor, and GitHub/package-specific repositories.

Core data manipulation and plotting packages:

```text
dplyr, tidyr, tibble, stringr, purrr, readxl, openxlsx, xlsx,
ggplot2, ggpubr, ggrepel, ggpointdensity, patchwork, viridis,
RColorBrewer, scales
```

Proteomics, statistics, and enrichment packages:

```text
proDA, limma, clusterProfiler, org.Mm.eg.db, AnnotationDbi,
diann, protr, enrichR, rstatix, TOSTER, vegan, MASS,
fitdistrplus, nortest, matrixStats
```

Heatmap and spatial packages:

```text
ComplexHeatmap, circlize, corrplot, ggcorrplot, sf, gridExtra
```

Single-cell and file-format packages:

```text
Seurat, harmony, Matrix, hdf5r, arrow, tidyverse
```

A typical installation pattern is:

```r
install.packages(c(
  "dplyr", "tidyr", "tibble", "stringr", "purrr", "readxl",
  "openxlsx", "xlsx", "ggplot2", "ggpubr", "ggrepel",
  "ggpointdensity", "patchwork", "viridis", "RColorBrewer",
  "scales", "rstatix", "TOSTER", "vegan", "MASS",
  "fitdistrplus", "nortest", "matrixStats", "corrplot",
  "ggcorrplot", "sf", "gridExtra", "hdf5r", "arrow"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "ComplexHeatmap", "limma", "clusterProfiler", "org.Mm.eg.db",
  "AnnotationDbi", "proDA"
))

install.packages("Seurat")
install.packages("harmony")
```


## Reproducibility

Several scripts set a random seed before clustering or dimensionality reduction. Published cluster assignments are additionally preserved through saved k-means RDS objects. Do not delete or overwrite these files when reproducing the published figures unless the input matrices or clustering parameters have intentionally changed.

Each script creates a timestamped log file in:

```text
Output/Logs/
```

Warnings, messages, uncaught errors, and tracebacks are written to the log. The scripts suppress warnings and messages in the interactive console, so inspect the log when execution fails or produces an unexpected result.



## Notes and known constraints

- All scripts assume they are launched from the repository root.
- Input column names and sample-name formats are tightly coupled to the original exports.
- Several QC exclusions and anatomical assignments are manually encoded in the scripts.
- UniProt-to-gene mapping depends on the installed version of `org.Mm.eg.db`.
- GO enrichment results may change slightly with annotation-database versions.
- The half-life estimation itself is external to this repository; the scripts prepare its input and import its output.
- Some scripts contain figure labels or filenames inherited from earlier versions of the analysis. The executed calculations, rather than the plot title text, should be used to interpret the data.
- `xlsx` requires Java. If Java-related errors occur, use a compatible Java installation or replace `xlsx::write.xlsx()` with `openxlsx::write.xlsx()`.



## Citation and data availability

This repository contains the analysis code associated with the manuscript:

**Couce-Iglesias M., Oliani E. et al. _Recording proteome dynamics in the late-stage embryo at hour-timescale, spatial, and single-cell resolution with MEMBRYO._ Manuscript in preparation.**

The final citation and repository accession numbers will be added once available.

External datasets used in the analysis must be cited separately, including the embryonic neocortex single-cell RNA-seq dataset from Di Bella et al. (2021).
