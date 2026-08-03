#########################################################################################################################################################################
# Full-slice spatial proteomics preprocessing and QC
#
# Purpose:
#   Load DIA-NN output from the E14 full-slice spatial proteomics experiment, annotate SILAC channels,
#   construct peptide/protein matrices, calculate heavy-to-light ratios, and generate protein-detection QC plots.
#
# GitHub/reuse notes:
#   - Keep the DIA-NN report and associated input files in the format used for the original analysis.
#   - Update the paths in "User settings" before running on another machine.
#########################################################################################################################################################################


# ---- Logging and warning control ----
log_dir <- "./Output/Logs/"
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  log_dir,
  paste0("04b_Fig4_SpatialProteomics_fullslice_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)

cat(
  "Full-slice spatial proteomics analysis log\n",
  "Started: ", format(Sys.time()), "\n",
  "Working directory: ", getwd(), "\n\n",
  file = log_file,
  sep = ""
)

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

options(
  error = function() {
    cat(
      "[ERROR] ", format(Sys.time()), "\n",
      geterrmessage(), "\n\n",
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
  library(diann)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
})


# ---- User settings ----
path <- "./Output/"
data_path <- "./Data/SILAC_SpatialProteomics_MS/"

dir.create(path, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(path, "Objects"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(path, "Figures"), recursive = TRUE, showWarnings = FALSE)

report_file <- file.path(data_path, "report_SpatialSILAC_Hemisphere.tsv")


# ---- Load and annotate DIA-NN output ----
data <- diann_load(report_file)

# Simplify sample names.
data$File.Name <- gsub(
  paste0(
    "\\\\\\\\10\\.64\\.161\\.152\\\\share\\\\mass_data\\\\Exploris480\\\\2024_10\\\\",
    "20241029_EO246_10_100um_serial\\\\20241025_EO246_10_100um_serial_|",
    "\\.raw$|",
    "\\\\\\\\10\\.64\\.161\\.152\\\\share\\\\mass_data\\\\Exploris480\\\\2024_10\\\\",
    "20241023_EO_bestBrain\\\\20241022_EO_"
  ),
  "",
  data$File.Name
)

data$File.Name <- gsub("SILAC", "bulk", data$File.Name)

# Remove bulk reference samples.
data <- data %>%
  filter(!File.Name %in% c("bulk_A1", "bulk_A2"))

# Map well names to their position in the serial-section grid.
grid <- data.frame(
  Well = c(
    "C11", "C12", "D1", "D2", "D3", "D5", "D6", "D7", "D8", "D9",
    "D10", "D11", "D12", "E1", "E2", "E3", "E4", "E7", "E8", "E9",
    "E10", "E11", "E12", "F3", "F4", "F5", "F6", "F7", "F8", "F9",
    "F11", "F12", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8",
    "G9", "G10", "G11", "G12", "H1", "H2", "H3", "H4", "H5", "H6",
    "H7", "H8", "H9", "H10", "H11", "H12", "A1", "A2", "A3", "A4",
    "A5", "A6", "A8", "A9", "A10", "A11", "A12", "B1", "B2", "B6",
    "B7", "B8", "B9", "B10", "C9"
  ),
  Number = 1:75
)

data$Number_grid <- grid$Number[match(data$File.Name, grid$Well)]

# Assign SILAC analogue and sample identifiers.
peptides <- data %>%
  mutate(
    Analogue = case_when(
      str_detect(Modified.Sequence, "\\(SILAC-(R|K)-M\\)") ~ "medium",
      str_detect(Modified.Sequence, "\\(SILAC-(R|K)-H\\)") ~ "heavy",
      TRUE ~ "light"
    ),
    ID = paste(File.Name, Number_grid, Analogue, sep = "_")
  )

# Treat different precursor charge states as separate peptide features.
peptides$Stripped.Sequence <- paste(
  peptides$Stripped.Sequence,
  peptides$Precursor.Charge,
  sep = "_"
)

# Remove contaminants.
peptides <- peptides %>%
  filter(!grepl("Cont|cont", Protein.Group))


# ---- Filter proteins supported by at least two precursors ----
peptides_min2 <- peptides %>%
  add_count(ID, Protein.Names, name = "Repetitions") %>%
  filter(Repetitions > 1) %>%
  dplyr::select(-Repetitions)


# ---- Create unnormalised protein matrices and heavy/light ratios ----
create_protein_matrix <- function(peptide_data) {
  protein_matrix <- peptide_data %>%
    group_by(ID, Protein.Names) %>%
    summarise(
      Precursor.Translated = median(Precursor.Translated, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = ID,
      values_from = Precursor.Translated
    ) %>%
    as.data.frame()

  protein_matrix[protein_matrix == 0] <- NA

  sample_ids <- unique(
    vapply(
      str_split(names(protein_matrix)[-1], "_"),
      function(x) paste(x[-length(x)], collapse = "_"),
      character(1)
    )
  )

  for (sample_id in sample_ids) {
    heavy_col <- paste(sample_id, "heavy", sep = "_")
    light_col <- paste(sample_id, "light", sep = "_")
    ratio_col <- paste(sample_id, "ratio", sep = "_")

    if (all(c(heavy_col, light_col) %in% names(protein_matrix))) {
      protein_matrix[[ratio_col]] <- protein_matrix[[heavy_col]] / protein_matrix[[light_col]]
    }
  }

  protein_matrix
}

prot_matrix_min2 <- create_protein_matrix(peptides_min2)

saveRDS(
  prot_matrix_min2,
  file.path(path, "Objects", "protMatrix_SpatialProteomics_hemisphere.rds")
)

rm(data, grid, pep_matrix_min2, peptides_min2)


# ---- Fig 4a ----
# % of detected proteins in the heavy channel (QC heatmap)
intensity_data <- prot_matrix_min2[, grepl("light|heavy", names(prot_matrix_min2)), drop = FALSE]
intensity_data[intensity_data == 0] <- NA

non_na_counts <- data.frame(
  Sample = names(colSums(!is.na(intensity_data))),
  Count = as.numeric(colSums(!is.na(intensity_data)))
) %>%
  separate(
    Sample,
    into = c("Well", "Number_grid", "Analogue"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(Number_grid = as.numeric(Number_grid))


# Spatial mapping
map_template <- data.frame(
  col1 = c(NA, NA, 13, 19, 25, 32, 40, 48, 55, NA, NA),
  col2 = c(NA, 6, 14, 20, 26, 33, 41, 49, 57, 63, NA),
  col3 = c(NA, 7, 15, 21, 27, 34, 42, 50, 57, 64, NA),
  col4 = c(1, 8, 16, 22, 28, 35, 43, 51, 58, 65, 70),
  col5 = c(2, 9, 17, 23, 29, 36, 44, 52, 59, 66, 71),
  col6 = c(3, 10, NA, NA, 30, 37, 45, 53, 60, 67, 72),
  col7 = c(4, 11, NA, NA, NA, 38, 46, 54, 61, 68, 73),
  col8 = c(5, 12, 18, 24, 31, 39, 47, 55, 62, 69, 74)
)

map_counts_to_template <- function(map_df, count_df) {
  map_df[] <- lapply(
    map_df,
    function(column) {
      vapply(
        column,
        function(value) {
          if (is.na(value)) {
            return(NA_real_)
          }

          matched_count <- count_df$Count[match(value, count_df$Number_grid)]
          ifelse(length(matched_count) == 0 || is.na(matched_count), NA_real_, matched_count)
        },
        numeric(1)
      )
    }
  )

  map_df
}

non_na_counts_light <- non_na_counts %>% filter(Analogue == "light")
non_na_counts_heavy <- non_na_counts %>% filter(Analogue == "heavy")

map_matrix_light <- map_counts_to_template(map_template, non_na_counts_light)
map_matrix_heavy <- map_counts_to_template(map_template, non_na_counts_heavy)


# Fraction of detected proteins in the heavy channel
map_matrix_ratio <- map_matrix_heavy / (map_matrix_light + map_matrix_heavy)

map_matrix_ratio_long <- map_matrix_ratio %>%
  pivot_longer(everything(), names_to = "Column", values_to = "Value") %>%
  group_by(Column) %>%
  mutate(Row = row_number()) %>%
  ungroup() %>%
  mutate(
    Row = factor(Row, levels = rev(1:11)),
    Column = gsub("col", "", Column)
  ) %>%
  # H5/position 49 was absent from the original DIA-NN data.
  filter(!(Column == "2" & Row == "8"))

p <- ggplot(map_matrix_ratio_long, aes(x = Column, y = Row, fill = Value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "white",
    high = "red4",
    na.value = "grey50",
    limits = c(0, max(map_matrix_ratio_long$Value, na.rm = TRUE))
  ) +
  labs(
    title = "Fraction of detected proteins in the heavy channel",
    x = "Column",
    y = "Row",
    fill = "Heavy / total"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = file.path(path, "Figures", "Fig4a_FractionDetectedProteins_heavyChannel.pdf"),
  plot = p,
  width = 10,
  height = 7
)

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