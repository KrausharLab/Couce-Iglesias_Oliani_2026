#########################################################################################################################################################################
# SILAC and Controls time-course & SILAC titration samples preprocessing
#
# Purpose:
#   Build the long-format SILAC proteomics table, remove QC outliers, normalise heavy signal by free amino-acid availability, filter proteins, and calculate iBAQ/riBAQ.
#   Build the long-format control proteomics table
#   Build the long-format SILAC titration table
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
  paste0("02a_PreprocessingSILAC&controlSamples_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)

# Create the log file and write session information
cat(
  "Preprocessing of SILAC and control samples log\n",
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
  library(MASS)
  library(RColorBrewer)
  library(circlize)
  library(dplyr)
  library(fitdistrplus)
  library(ggplot2)
  library(ggpubr)
  library(nortest)
  library(protr)
  library(purrr)
  library(readxl)
  library(rstatix)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(viridis)
})

# ---- User settings ----
path <- "./Output/"
dir.create(path, recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(path, "Objects/"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(path, "Figures/"), recursive = TRUE, showWarnings = FALSE)

data_path <- "./Data/SILAC_Timecourse_Bulk_MS/"
data_path_titration <- "./Data/SILAC_Titration_MS/"


### SILAC timecourse #################################################
# Load data tables from the original analysis.
data_requantifyOn <- read.csv(file.path(data_path, "SILACtimecourse_ratios.csv"))
data_L_requantifyOn <- read.csv(file.path(data_path, "SILACtimecourse_light.csv"))
data_H_requantifyOn <- read.csv(file.path(data_path, "SILACtimecourse_pulse.csv"))
data_T_requantifyOn <- read.csv(file.path(data_path, "SILACtimecourse_total.csv"))


# ---- Preprocess input tables ----
samples <- read_excel(file.path(data_path, "SILACtimecourse_samples_metadata.xlsx"))
names(samples) <- gsub(" ", "_", names(samples))
names(samples)[1] <- "Sample_ID_Selbach"
samples$ID <- apply(samples[,-2], 1, paste, collapse="_")
samples$ID <- gsub(" ", "", samples$ID)
samples$Sample_ID_Selbach <- gsub(" ", "", samples$Sample_ID_Selbach)
samples <- drop_na(samples)

genes <- unique(rbind(data_requantifyOn[, c("protein_group", "genes")], data_L_requantifyOn[, c("protein_group", "genes")], 
                      data_H_requantifyOn[, c("protein_group", "genes")],  data_T_requantifyOn[, c("protein_group", "genes")]))

process_files <- function(data, value_name = "name") {
  data %>%
    dplyr::select(-c(X, genes)) %>%
    pivot_longer(!protein_group, names_to = "Run", values_to = value_name) %>%
    mutate(Run = gsub("Scooter_20241209_RJK_MMdia_dilutedKrausharLab_", "", Run),
           protein_group = gsub(";.*", "", protein_group))
}

data_requantifyOn_long <- process_files(data_requantifyOn, value_name = "ratios_HL")
data_L_requantifyOn <- process_files(data_L_requantifyOn, value_name = "LFQ_L")
data_H_requantifyOn <- process_files(data_H_requantifyOn, value_name =  "LFQ_H")
data_T_requantifyOn <- process_files(data_T_requantifyOn, value_name = "LFQ_T")

data_requantifyOn_long <- merge(samples, data_requantifyOn_long, by.x= "Sample_ID_Selbach", by.y= "Run", all.y = TRUE)
data_L_requantifyOn <- merge(samples, data_L_requantifyOn, by.x= "Sample_ID_Selbach", by.y= "Run", all.y = TRUE)
data_H_requantifyOn <- merge(samples, data_H_requantifyOn, by.x= "Sample_ID_Selbach", by.y= "Run", all.y = TRUE)
data_T_requantifyOn <- merge(samples, data_T_requantifyOn, by.x= "Sample_ID_Selbach", by.y= "Run", all.y = TRUE)

# Calculate heavy/total and light/total ratios from the H/L ratio.
data_requantifyOn_long$ratios_HT <- data_requantifyOn_long$ratios_HL/(data_requantifyOn_long$ratios_HL+1)
data_requantifyOn_long$ratios_LT <- 1/(data_requantifyOn_long$ratios_HL+1)

# Merge ratio, light, heavy and total intensity tables into one long table.
data_requantifyOn_long <- merge(data_requantifyOn_long, data_L_requantifyOn[, c("protein_group", "ID", "LFQ_L")], by=c("protein_group", "ID"), all = TRUE)
data_requantifyOn_long <- merge(data_requantifyOn_long, data_H_requantifyOn[, c("protein_group", "ID", "LFQ_H")], by=c("protein_group", "ID"), all = TRUE)
data_requantifyOn_long <- merge(data_requantifyOn_long, data_T_requantifyOn[, c("protein_group", "ID", "LFQ_T")], by=c("protein_group", "ID"), all = TRUE)

data_requantifyOn_long <- merge(data_requantifyOn_long, genes, by="protein_group", all.x = TRUE)

data_requantifyOn_long <- data_requantifyOn_long %>% # Eliminate rows with all NAs
  filter(if_any(-c(protein_group, ID), ~ !is.na(.))) %>%
  separate(ID, into = c("Sample_ID_Selbach", "Embryonic_Age", "Experiment_Number", "Embryo_Number", "Replicates", "Incubation_Time"), sep = "_", remove = FALSE)

rm(data_L_requantifyOn, data_H_requantifyOn, data_T_requantifyOn, genes)

# Remove manually selected poor-quality samples identified in the QC script.
outliers <- c("E12.5_8h_3", "E12.5_4h_5", "E13.5_4h_3", "E13.5_8h_4", "E13.5_8h_1", "E16.5_8h_1")

data_requantifyOn_long_noOutliers <- data_requantifyOn_long %>%
    mutate(sample = paste(Embryonic_Age, Incubation_Time, Replicates, sep = "_")) %>%
    filter(!sample %in% outliers)

rm(outliers)


# ---- Normalise by free amino-acid availability ----
# Load precursor-level quantities used for free-AA normalisation.
precursors <- read.csv(file.path(data_path, "SILACtimecourse_precursors.csv"))
precursors$Run <- gsub("Scooter_20241209_RJK_MMdia_dilutedKrausharLab_", "", precursors$Run)
precursors <- merge(samples, precursors, by.x= "Sample_ID_Selbach", by.y= "Run", all.y = TRUE)

# Apply the original age-specific free-AA percentages and summarise precursor values per protein/sample.
precursors <- precursors %>%
  dplyr::select(-X) %>%
  mutate(
    precursor_quantity_pulse_normFreeAA = precursor_quantity_pulse / case_when(
      Embryonic_Age == "E12.5" ~ 0.7594,
      Embryonic_Age == "E13.5" ~ 0.7184,
      Embryonic_Age == "E14.5" ~ 0.6796,
      Embryonic_Age == "E15.5" ~ 0.6429,
      Embryonic_Age == "E16.5" ~ 0.6082,
      TRUE ~ NA_real_
    ),
    precursors_ratios_HT_normFreeAA =  precursor_quantity_pulse_normFreeAA / (precursor_quantity_pulse + precursor_quantity_L)) %>%
    group_by(protein_group, Incubation_Time, Embryonic_Age,ID,Experiment_Number,Embryo_Number,Replicates,Sample_ID_Selbach) %>%
    summarise(ratios_HT_normFreeAA = median(precursors_ratios_HT_normFreeAA, na.rm = TRUE), 
              H_normFreeAA = median(precursor_quantity_pulse_normFreeAA, na.rm = TRUE),
              H = median(precursor_quantity_pulse, na.rm = TRUE),
              L = median(precursor_quantity_L, na.rm = TRUE), .groups = "drop")
precursors <- precursors[precursors$protein_group != "",]
precursors$protein_group <- gsub(":.*", "", precursors$protein_group)

# Add precursor-derived normalised quantities to the protein-group table.
data_requantifyOn_long_noOutliers_normfreeAA <- merge(data_requantifyOn_long_noOutliers, precursors[, c("protein_group", "ID", "ratios_HT_normFreeAA", "H_normFreeAA", "H", "L")], by=c("protein_group", "ID"), all = TRUE)

rm(precursors)


# ---- Filter proteins ----
data_requantifyOn_long_noOutliers_normfreeAA_filter <- data_requantifyOn_long_noOutliers_normfreeAA %>%
  group_by(protein_group, Incubation_Time, Embryonic_Age) %>%
  filter(sum(!is.na(ratios_HT)) >= 2) %>%
  ungroup()


# ---- Calculate iBAQ and riBAQ ----
precursors <- read.csv(file.path(data_path, "SILACtimecourse_precursors.csv"))

# Retrieve UniProt sequences for iBAQ peptide counting.
protein_groups <- unique(precursors$protein_group)
protein_groups <- gsub(":.*", "", protein_groups)  # Remove any suffixes after a colon

seq <- data.frame(protein_group=protein_groups, seq = rep(NA, length(protein_groups)))

for (i in 1:nrow(seq)) {
    print(i)
    tryCatch({
        seq$seq[i] <- getUniProt(seq$protein_group[i])
    }, error = function(e) {
        message("Error at iteration ", i, ": ", e$message)
        seq$seq[i] <- NA  # or any other placeholder indicating failure
    })
}

# Count theoretically observable tryptic peptides for iBAQ normalisation.
 tryptic_digest <- function(sequence) {
  # Split the sequence at K or R not followed by P
  peptides <- unlist(strsplit(sequence, "(?<=[KR])(?!P)", perl = TRUE))
  return(peptides)
}

count_tryptic_peptides <- function(sequence, min_length = 6, max_length = 30) {
  if (is.na(sequence) || sequence == "") {
    return(NA)  # or return(0) if you prefer
  }
  
  # Perform tryptic digestion
  peptides <- tryptic_digest(sequence)
  
  # Filter peptides by length
  peptides <- peptides[nchar(peptides) >= min_length & nchar(peptides) <= max_length]
  
  return(length(peptides))
}

# Apply the function using sapply
seq$observable_peptides <- sapply(seq$seq, count_tryptic_peptides)

seq <- drop_na(seq)

saveRDS(seq, file.path(path, "Objects/ProteomicsTimecourse_SILAC_ProteinSequences.rds"))


# Calculate iBAQ values for each protein group and sample
iBAQ <- precursors %>%
  filter(protein_group != "") %>%
  separate(protein_group, into = c("protein_group", "Gene"), sep = ":") %>%
  group_by(protein_group, Run) %>%
  summarise(iBAQ_H = sum(precursor_quantity_pulse, na.rm = TRUE) / unique(seq$observable_peptides[match(protein_group, seq$protein_group)]), 
          iBAQ_L = sum(precursor_quantity_L, na.rm = TRUE) / unique(seq$observable_peptides[match(protein_group, seq$protein_group)])) %>%
  ungroup()

iBAQ <- iBAQ %>%
  mutate(iBAQ_H = na_if(iBAQ_H, 0),
         iBAQ_L = na_if(iBAQ_L, 0))

# Add sample metadata back after iBAQ calculation.
iBAQ$Run <- gsub("Scooter_20241209_RJK_MMdia_dilutedKrausharLab_", "", iBAQ$Run)

iBAQ <- merge(samples, iBAQ, by.x= "Sample_ID_Selbach", by.y= "Run", all.y = TRUE)

iBAQ <- iBAQ[iBAQ$protein_group != "",]

# Convert iBAQ values to relative iBAQ within each sample.
riBAQ <- iBAQ %>%
  group_by(ID) %>%
  mutate(
    total_iBAQ_H = sum(iBAQ_H, na.rm = TRUE),
    total_iBAQ_L = sum(iBAQ_L, na.rm = TRUE),
    riBAQ_H = iBAQ_H / total_iBAQ_H,
    riBAQ_L = iBAQ_L / total_iBAQ_L
  ) %>%
  ungroup()

# Save the final filtered table used downstream.
data_requantifyOn_long_noOutliers_normfreeAA_filter <- merge(data_requantifyOn_long_noOutliers_normfreeAA_filter, riBAQ[, c("protein_group", "ID", "riBAQ_H", "riBAQ_L", "iBAQ_H", "iBAQ_L")], by=c("protein_group", "ID"), all.x = TRUE)

saveRDS(data_requantifyOn_long_noOutliers_normfreeAA_filter, file.path(path, "Objects/ProteomicsTimecourse_SILAC.rds"))



### Control timecourse #################################################
# ---- Preprocess control timecourse ----
data_L_requantifyOn <- read.csv(file.path(data_path, "Controls_light.csv"))
data_H_requantifyOn <- read.csv(file.path(data_path, "Controls_pulse.csv"))

genes <- unique(rbind(data_L_requantifyOn[, c("protein_group", "genes")], data_H_requantifyOn[, c("protein_group", "genes")]))

data_L_requantifyOn <- data_L_requantifyOn %>% 
    dplyr::select(-c(X, genes)) %>%
    pivot_longer(!protein_group, names_to = "Run", values_to = "LFQ_L")

data_H_requantifyOn <- data_H_requantifyOn %>% 
    dplyr::select(-c(X, genes)) %>%
    pivot_longer(!protein_group, names_to = "Run", values_to = "LFQ_H")

data_L_requantifyOn$protein_group <- gsub(";.*", "", data_L_requantifyOn$protein_group)
data_H_requantifyOn$protein_group <- gsub(";.*", "", data_H_requantifyOn$protein_group)

data_requantifyOn <- merge(data_H_requantifyOn, data_L_requantifyOn, by.x= c("protein_group", "Run"), all = TRUE)

data_requantifyOn_long <- data_requantifyOn %>% # Eliminate rows with all NAs
  filter(if_any(-c(protein_group, Run), ~ !is.na(.))) 

# Keep proteins detected in at least two heavy or light values per age, matching the original filter.
data_requantifyOn_long <- data_requantifyOn_long %>%
    separate(Run, into = c("Embryonic_Age", "Replicates"), sep = "_") %>%
    group_by(Embryonic_Age) %>%
    filter(sum(!is.na(LFQ_L)) >= 2 & sum(!is.na(LFQ_H)) >= 2) %>%
    ungroup()

saveRDS(data_requantifyOn_long, file.path(path, "Objects/ProteomicsTimecourse_Control.rds"))



### SILAC titration #################################################
# data_path_titration instead pathway
Exp195 = file.path(data_path_titration,"Exp195_proteinGroups.txt")
Exp182 = file.path(data_path_titration,"Exp182_proteinGroups.txt")
Exp173 = file.path(data_path_titration,"Exp173_proteinGroups.txt")

# Extract the LFQ-normalised protein matrix from one MaxQuant experiment
make_LFQ_matrix <- function(file, experiment) {
  read.delim(file, check.names = TRUE) %>%
    filter(Potential.contaminant != "+", Reverse != "+") %>%
    select(Protein.IDs, any_of("Gene.names"), starts_with("LFQ.intensity."), starts_with("Ratio.H.L.normalized")) %>%
    mutate(
      across(any_of(c("Protein.Names", "Gene.names")),~ gsub(";.*", "", .x))) %>%
    rename_with(
      ~ str_replace_all(
        .x,
        c(
          "^LFQ.intensity.H." = paste0(experiment, "_LFQ_Heavy_"),
          "^LFQ.intensity.L." = paste0(experiment, "_LFQ_Light_"),
          "^Ratio.H.L.normalized." = paste0(experiment, "_Ratios_HL_"),
          "_NC|_173" = "_4h",
          "x(?=\\d)" = "x_"
        )
      )
    )
}

Exp182_2 <- make_LFQ_matrix(Exp182, "182")
Exp173_2 <- make_LFQ_matrix(Exp173, "173")

Exp195_2 <- read.delim(Exp195, check.names = TRUE) %>%
  filter(Potential.contaminant != "+", Reverse != "+") %>%
  select(Protein.IDs, starts_with("LFQ.intensity."), starts_with("Ratio.H.L.normalized")) %>%
  select(!contains("SP")) %>%
  mutate(Protein.IDs = str_extract(Protein.IDs, "(?<=\\|)[^|;]+(?=\\|)")) %>%
  rename_with(
    ~ str_replace_all(
      .x,
      c(
        "^LFQ.intensity.H." = "195_LFQ_Heavy_",
        "^LFQ.intensity.L." = "195_LFQ_Light_",
        "^Ratio.H.L.normalized." = "195_Ratios_HL_",
        "NC" = "4h",
        "x(?=\\d)" = "x_"
      )
    )
)

titration <- merge(Exp173_2, Exp182_2, by= c("Protein.IDs", "Gene.names"), all=TRUE)
titration <- merge(titration, Exp195_2, by= "Protein.IDs", all=TRUE)

titration[titration == 0] <- NA

titration_long <- titration %>%
  mutate(Protein.IDs = gsub(";.*", "", Protein.IDs)) %>%
  filter(!is.na(Protein.IDs)) %>%
  pivot_longer(cols = !c(Protein.IDs,Gene.names), names_to = "Samples") %>%
  separate(Samples, into = c("Experiment", "Data_Type", "Analogue", "Incubation_Time", "Concentration", "Replicate"), sep="_", remove = FALSE) %>%
  mutate(Concentration = gsub("X", "x", Concentration))

# Filter by 2 ratios
protein_to_keep <- titration_long %>%
  group_by(Protein.IDs, Analogue, Concentration) %>%
  filter(sum(Data_Type == "Ratios" & !is.na(value)) >= 2) %>%
  ungroup() %>%
  drop_na(value)

titration_long_filtered <- titration_long %>%
  filter(
    Protein.IDs %in% unique(protein_to_keep$Protein.IDs),
    Data_Type == "LFQ"
  ) %>%
  mutate(
    Concentration = factor(
      Concentration,
      levels = c("noRoller", "Rep", "1x", "3x", "5x", "10x", "20x")
    )
  )


saveRDS(titration_long_filtered, file.path(path, "Objects/ProteomicsTitration_SILAC.rds"))



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
