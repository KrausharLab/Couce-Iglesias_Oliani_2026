#########################################################################################################################################################################
# Multiomics integration: RNAseq, Riboseq and SILAC proteomics
#
# Purpose:
#   Integrate transcriptome, translatome and SILAC proteomics outputs for cross-omics comparison
#
# GitHub/reuse notes:
#   - Keep input files in the same format as the original analysis.
#   - Update the paths in "User settings" before running on another machine.
#########################################################################################################################################################################


# ---- Logging and warning control ----
log_dir <- "./Output/Logs/"
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  log_dir,
  paste0("03a_IntegrationMultiomics_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)

# Create the log file and write session information
cat(
  "Integration of multiomics data log\n",
  "Started: ", format(Sys.time()), "\n",
  "Working directory: ", getwd(), "\n\n",
  file = log_file,
  sep = ""
)

# Record warnings without printing them in the console
globalCallingHandlers(
  warning = function(w) {
    cat(
      "[WARNING] ", format(Sys.time()), "\n",
      conditionMessage(w), "\n\n",
      file = log_file,
      append = TRUE,
      sep = ""
    )

    invokeRestart("muffleWarning")
  },

  message = function(m) {
    cat(
      "[MESSAGE] ", format(Sys.time()), "\n",
      conditionMessage(m), "\n\n",
      file = log_file,
      append = TRUE,
      sep = ""
    )

    invokeRestart("muffleMessage")
  }
)

# Record uncaught errors before R stops
options(
  error = function() {
    error_message <- geterrmessage()

    cat(
      "[ERROR] ", format(Sys.time()), "\n",
      error_message, "\n\n",
      "Traceback:\n",
      file = log_file,
      append = TRUE,
      sep = ""
    )

    capture.output(
      traceback(2),
      file = log_file,
      append = TRUE
    )

    cat(
      "\nScript terminated: ", format(Sys.time()), "\n",
      file = log_file,
      append = TRUE,
      sep = ""
    )
  }
)

cat("Log file: ", normalizePath(log_file, mustWork = FALSE), "\n")


# ---- Loading packages ----
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(clusterProfiler)
  library(dplyr)
  library(ggplot2)
  library(ggpointdensity)
  library(ggpubr)
  library(ggrepel)
  library(openxlsx)
  library(org.Mm.eg.db)
  library(readxl)
  library(scales)
  library(stringr)
  library(tidyr)
  library(tidyverse)
  library(viridis)
})


# ---- User settings ----
path <- "./Output/"
dir.create(path, recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(path, "Objects/"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(path, "Figures/"), recursive = TRUE, showWarnings = FALSE)

data_path <- "./Data/SILAC_Timecourse_Bulk_MS/"

set.seed(123)


# ---- Prepare SILAC proteomics outputs ----
timecourse_silac <- readRDS(file.path(path, "Objects/ProteomicsTimecourse_SILAC.rds"))

timecourse_silac <- timecourse_silac %>%
  dplyr::select(-c("ID", "Experiment_Number","Embryo_Number","Sample_ID_Selbach","sample")) %>%
  pivot_longer(cols=c(ratios_HL,ratios_HT,ratios_LT,LFQ_L,LFQ_H,LFQ_T,ratios_HT_normFreeAA,H_normFreeAA,H,L,riBAQ_H,riBAQ_L,iBAQ_H,iBAQ_L), names_to = "Experiment", values_to = "Value") 


# ---- Prepare RNA-seq and Ribo-seq data ----
rna_ribo <- openxlsx::read.xlsx(file.path(data_path, "RNAseq_Riboseq.xlsx"), sheet="tpm_data")

# Convert RNA/Ribo replicate columns into long format.
rna_ribo <- rna_ribo %>%
  pivot_longer(
    cols = matches("_(ribo|total)_\\d+$"),
    names_to = c("Embryonic_Age", "Experiment", "Replicates"),
    names_pattern = "^(.*)_(ribo|total)_(\\d+)$",
    values_to = "Value"
  )

rna_ribo$Experiment <- gsub("ribo", "RPF", rna_ribo$Experiment)
rna_ribo$Experiment <- gsub("total", "RNAseq", rna_ribo$Experiment)

rna_ribo <- expand_grid(
  Incubation_Time = c("4h", "8h")
) %>%
  mutate(data = list(rna_ribo)) %>%
  unnest(data)


# ---- Prepare t12 data ----
t12_E125 <- read_excel(file.path(data_path, "SILACtimecourse_t12_E12.5.xlsx"))
t12_E135 <- read_excel(file.path(data_path, "SILACtimecourse_t12_E13.5.xlsx"))
t12_E145 <- read_excel(file.path(data_path, "SILACtimecourse_t12_E14.5.xlsx"))
t12_E155 <- read_excel(file.path(data_path, "SILACtimecourse_t12_E15.5.xlsx"))
t12_E165 <- read_excel(file.path(data_path, "SILACtimecourse_t12_E16.5.xlsx"))


clean_df <- function(df) {
    df <- as.data.frame(df)
    
    # Keep only columns whose first row is NOT NA
    df <- df[, !is.na(df[1, ]), drop = FALSE]

    # Use the first row as column names
    colnames(df) <- as.character(unlist(df[1, ], use.names = FALSE))
    df <- df[-1, , drop = FALSE]

    # Keep proteins and t1/2
    df <- df[, c("protein_group", "t1/2")]

    # Convert t1/2 to numeric
    df$`t1/2` <- as.numeric(as.character(df$`t1/2`))

    df
}

t12_E125 <- clean_df(t12_E125)
t12_E135 <- clean_df(t12_E135)
t12_E145 <- clean_df(t12_E145)
t12_E155 <- clean_df(t12_E155)
t12_E165 <- clean_df(t12_E165)

# Add a column to identify the stage
t12_E125$Embryonic_Age <- "E12.5"
t12_E135$Embryonic_Age <- "E13.5"
t12_E145$Embryonic_Age <- "E14.5"
t12_E155$Embryonic_Age <- "E15.5"
t12_E165$Embryonic_Age <- "E16.5"

# Combine all dataframes into one
t12 <- bind_rows(t12_E125, t12_E135, t12_E145, t12_E155, t12_E165)

# Prepare to bind
t12 <- t12 %>%
  dplyr::select(protein_group, Embryonic_Age, Value = `t1/2`) %>%
  unique() %>%
  mutate(Experiment = "t12", 
          Replicates = NA) 

t12 <- expand_grid(
  Incubation_Time = c("4h", "8h")
) %>%
  mutate(data = list(t12)) %>%
  unnest(data)


# ---- Join proteomics, t12, RNA-seq and Ribo-seq tables ----
# Map UniProt IDs to gene symbols for cross-omics joins.
mapped_ids <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = unique(timecourse_silac$protein_group),
  columns = c("ENSEMBL", "SYMBOL"),
  keytype = "UNIPROT"
)

timecourse_silac <- timecourse_silac %>%
  dplyr::left_join(mapped_ids, by = c("protein_group" = "UNIPROT")) %>%
  dplyr::rename(Ensembl_id = "ENSEMBL")

rna_ribo <- rna_ribo %>%
  dplyr::left_join(mapped_ids, by = c("gene_id" = "ENSEMBL")) %>%
  dplyr::rename(genes = "SYMBOL", protein_group = "UNIPROT", Ensembl_id = "gene_id")

t12 <- t12 %>%
  dplyr::left_join(mapped_ids, by = c("protein_group" = "UNIPROT")) %>%
  dplyr::rename(genes = "SYMBOL", Ensembl_id = "ENSEMBL")

timecourse_silac <- timecourse_silac[,c("protein_group", "genes", "Ensembl_id", "Embryonic_Age", "Incubation_Time", "Replicates", "Experiment", "Value")]
rna_ribo <- rna_ribo[,c("protein_group", "genes", "Ensembl_id", "Embryonic_Age", "Incubation_Time", "Replicates", "Experiment", "Value")]
t12 <- t12[,c("protein_group", "genes", "Ensembl_id", "Embryonic_Age", "Incubation_Time", "Replicates", "Experiment", "Value")]

all_data <- rbind(timecourse_silac, rna_ribo, t12)

saveRDS(all_data, file.path(path, "Objects/all_data.rds"))

rm(timecourse_silac, rna_ribo, t12, mapped_ids, t12_E125, t12_E135, t12_E145, t12_E155, t12_E165)


# ---- Z-score proteomics, t12, RNA-seq and Ribo-seq tables ----
data_4h <- all_data %>%
  filter(Incubation_Time == "4h" & Embryonic_Age != "P0" & Experiment %in% c("RPF", "RNAseq", "riBAQ_H", "LFQ_T", "ratios_HL")) %>%
  group_by(protein_group, genes, Embryonic_Age, Incubation_Time, Experiment) %>%
  filter(sum(!is.na(Value)) >= 2) %>%
  ungroup()

# Log2-transform replicate-level values
data_log <- data_4h %>%
  mutate(
    Experiment = gsub("riBAQ_|LFQ_|ratios_", "", Experiment),
    Age_numeric = as.numeric(sub("E", "", Embryonic_Age)),
    log2_value = log2(Value)
  )

# Median trajectory across replicates
trajectory_stats <- data_log %>%
  group_by(protein_group, genes, Experiment, Embryonic_Age, Age_numeric) %>%
  summarise(
    median_log2 = median(log2_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(protein_group, genes, Experiment) %>%
  mutate(
    trajectory_mean = mean(median_log2, na.rm = TRUE),
    trajectory_sd   = sd(median_log2, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  dplyr::select(
    protein_group, genes, Experiment, Embryonic_Age,
    trajectory_mean, trajectory_sd
  )

# Apply trajectory-derived Z-score parameters to each replicate
dat_z_replicates <- data_log %>%
  left_join(
    trajectory_stats,
    by = c("protein_group", "genes", "Experiment", "Embryonic_Age")
  ) %>%
  mutate(
    Zscore_log2_values = case_when(
      is.na(trajectory_sd) | trajectory_sd == 0 ~ 0,
      TRUE ~ (log2_value - trajectory_mean) / trajectory_sd
    )
  )

# Median trajectory and replicate SD at each age
dat_z_summary <- dat_z_replicates %>%
  group_by(protein_group, genes, Experiment, Embryonic_Age, Age_numeric) %>%
  summarise(
    Zscore_median = median(Zscore_log2_values, na.rm = TRUE),
    Zscore_mean   = mean(Zscore_log2_values, na.rm = TRUE),
    Zscore_SD     = sd(Zscore_log2_values, na.rm = TRUE),
    n             = sum(!is.na(Zscore_log2_values)),
    .groups = "drop"
  )

mat_wide_median <- dat_z_summary %>%
  unite(col = "Feat", Experiment, Embryonic_Age, sep = "_") %>%
  dplyr::select(protein_group, genes, Feat, Zscore_median) %>%
  pivot_wider(names_from = Feat, values_from = Zscore_median)

mat <- mat_wide_median %>%
    filter(!is.na(protein_group)) %>%
    mutate(protein_group = paste(protein_group, genes, sep = "_")) %>%
    column_to_rownames("protein_group") %>%
    dplyr::select(-genes) %>%
    as.matrix()

mat <- mat[complete.cases(mat), , drop = FALSE]

saveRDS(mat, file.path(path, "Objects/Zscore_RPF_RNAseq_riBAQH_T_HL_4h.rds"))


mat_wide_sd <- dat_z_summary %>%
  unite(col = "Feat", Experiment, Embryonic_Age, sep = "_") %>%
  dplyr::select(protein_group, genes, Feat, Zscore_SD) %>%
  pivot_wider(names_from = Feat, values_from = Zscore_SD)

mat_sd <- mat_wide_sd %>%
    filter(!is.na(protein_group)) %>%
    mutate(protein_group = paste(protein_group, genes, sep = "_")) %>%
    column_to_rownames("protein_group") %>%
    dplyr::select(-genes) %>%
    as.matrix()

mat_sd <- mat_sd[complete.cases(mat_sd), , drop = FALSE]

saveRDS(mat_sd, file.path(path, "Objects/Zscore_RPF_RNAseq_riBAQH_T_HL_4h_SD.rds"))

rm(list=ls())


cat(
  "\nScript completed successfully: ",
  format(Sys.time()),
  "\n",
  file = log_file,
  append = TRUE,
  sep = ""
)




#########################################################################################################################################################################
# END OF SCRIPT
#########################################################################################################################################################################
