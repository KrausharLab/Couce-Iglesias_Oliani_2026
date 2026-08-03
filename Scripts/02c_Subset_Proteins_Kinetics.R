#########################################################################################################################################################################
# Proteomics time-course subset preparation for kinetics
#
# Purpose:
#   Select proteins whose total abundance is stable between 4 h and 8 h.
#   Add the 0 h control samples.
#   Retain proteins measured at 0 h, 4 h and 8 h for each embryonic age.
#   Export one wide-format Excel file per embryonic age.
#
# GitHub/reuse notes:
#   - Run 02a_PreprocessingSILAC&controlSamples.R first.
#   - Update the paths in "User settings" before running on another machine.
#########################################################################################################################################################################


# ---- Logging and warning control ----
log_dir <- "./Output/Logs/"
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  log_dir,
  paste0("02c_Subset_Proteins_Kinetics_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)

# Create the log file and write session information
cat(
  "Data subset for protein kinetics calculation log\n",
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
  library(dplyr)
  library(limma)
  library(stringr)
  library(tidyr)
  library(xlsx)
})

# ---- User settings ----
path <- "./Output/"
objects_path <- file.path(path, "Objects")
figures_path <- file.path(path, "Figures")

padj_cutoff <- 0.05
log2fc_cutoff <- 0.6


# ---- Load data ----
prot <- readRDS(file.path(objects_path, "ProteomicsTimecourse_SILAC.rds"))

# Total abundance calculated from precursor-derived heavy and light intensities.
# Keep NA when both measurements are missing.
prot <- prot %>%
  mutate(
    sample = paste(Embryonic_Age, Incubation_Time, Replicates, sep = "_"),
    T = if_else(is.na(H) & is.na(L), NA_real_, coalesce(H, 0) + coalesce(L, 0))
  )


# ---- Select proteins in steady state between 4 h and 8 h ----
prot_total <- prot %>%
  filter(Incubation_Time %in% c("4h", "8h")) %>%
  dplyr::select(protein_group, Embryonic_Age, sample, T) %>%
  pivot_wider(names_from = sample, values_from = T)

ages <- sort(unique(prot$Embryonic_Age))
steady_proteins <- setNames(vector("list", length(ages)), ages)

for (age in ages) {
  age_data <- prot_total %>%
    filter(Embryonic_Age == age) %>%
    dplyr::select(-Embryonic_Age) %>%
    filter(if_any(-protein_group, ~ !is.na(.)))

  expression_matrix <- log2(as.matrix(age_data[, -1]))
  rownames(expression_matrix) <- age_data$protein_group

  timepoint <- str_extract(colnames(expression_matrix), "4h|8h") %>%
    factor(levels = c("4h", "8h"))

  design <- model.matrix(~ timepoint)
  fit <- eBayes(lmFit(expression_matrix, design))

  results <- topTable(
    fit,
    coef = "timepoint8h",
    number = Inf,
    sort.by = "none"
  )

  changed <- results$adj.P.Val < padj_cutoff &
    abs(results$logFC) >= log2fc_cutoff

  steady_proteins[[age]] <- rownames(results)[!changed | is.na(changed)]
}

# Keep each protein only in ages where it passed the steady-state criterion.
steady_table <- bind_rows(
  lapply(names(steady_proteins), function(age) {
    tibble(
      Embryonic_Age = age,
      protein_group = steady_proteins[[age]]
    )
  })
)

prot_filtered <- prot %>%
  inner_join(steady_table, by = c("Embryonic_Age", "protein_group"))

rm(prot_total, age_data, expression_matrix, timepoint, design, fit, results, changed, steady_proteins, steady_table)


# ---- Require measurements at 4 h and 8 h ----
# A protein must have at least one light and one heavy value at every time point.
prot_filtered <- prot_filtered %>%
  group_by(protein_group, Embryonic_Age) %>%
  filter(
    all(c("4h", "8h") %in% Incubation_Time[!is.na(L)]),
    all(c("4h", "8h") %in% Incubation_Time[!is.na(H)])
  ) %>%
  ungroup()


# ---- Export one table per embryonic age ----
prot_filtered_wide <- prot_filtered %>%
  pivot_longer(
    cols = c(L, H),
    names_to = "Analogue",
    values_to = "Intensity"
  ) %>%
  mutate(sample = paste(Embryonic_Age, Incubation_Time, Analogue, Replicates, sep = "_")) %>%
  dplyr::select(protein_group, sample, Intensity) %>%
  pivot_wider(names_from = sample, values_from = Intensity)

for (age in ages) {
  age_table <- prot_filtered_wide %>%
    dplyr::select(protein_group, contains(age)) %>%
    filter(if_any(-protein_group, ~ !is.na(.))) %>%
    as.data.frame()

  write.xlsx(
    age_table,
    file.path(objects_path, paste0("Timecourse_Controls_t12_", age, ".xlsx")),
    row.names = FALSE
  )
}

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