#########################################################################################################################################################################
# SILAC time-course quality control
#
# Purpose:
#   Generate QC plots for incorporation, intensity distributions, protein counts, outlier detection, PCA, replicate correlations, and free amino-acid normalisation.
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
  paste0("02b_Fig1_QCproteomics_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)

# Create the log file and write session information
cat(
  "Spatial proteomics analysis log\n",
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
  library(factoextra)
  library(fitdistrplus)
  library(ggplot2)
  library(ggpubr)
  library(nortest)
  library(patchwork)
  library(purrr)
  library(readxl)
  library(rstatix)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(tidyverse)
  library(vegan)
  library(viridis)
  library(openxlsx)
})

# ---- User settings ----
path <- "./Output/"
dir.create(path, recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(path, "Objects/"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(path, "Figures/"), recursive = TRUE, showWarnings = FALSE)

data <- readRDS(file.path(path, "Objects/ProteomicsTimecourse_SILAC.rds"))


# ---- Quality-control plots ----

# ---- Figure 1c ----
# Log2 LFQ titration
titration_long_filtered <- readRDS(file.path(path, "Objects/ProteomicsTitration_SILAC.rds"))

p <- ggplot(titration_long_filtered, aes(x = Concentration, y = log2(value))) +
      geom_boxplot(position="dodge") +
      labs(title = "Titration", x = "Samples", y = "Log2LFQ") +
      theme_minimal() +
      facet_grid(~ Analogue)

pdf(paste(path,"Figures/Fig1c_Boxplot_TitrationSILAC_Log2LFQ.pdf", sep = ""), width=10, height=6)
print(p)
dev.off()


# ---- Figure 1d ----
# heavy/total ratios over developmental time
df_plot <- data %>%
  dplyr::select(protein_group, Embryonic_Age, Incubation_Time, Replicates, ratios_HT, ratios_LT) %>%
  pivot_longer(cols = c(ratios_HT, ratios_LT), names_to = "Analogue", values_to = "Ratio") %>%
  mutate(
    Embryonic_Age = factor(Embryonic_Age, levels = c("E12.5","E13.5","E14.5","E15.5","E16.5")),
    Incubation_Time = factor(Incubation_Time, levels = c("4h","8h"))
  ) %>%
  filter(is.finite(Ratio))  %>%
  group_by(Embryonic_Age) %>%
  filter(n_distinct(Incubation_Time[!is.na(Ratio)]) == 2) %>%
  ungroup()

# Test 4h versus 8h incorporation within each embryonic age
stats <- df_plot %>%
  group_by(Embryonic_Age, Analogue) %>%
  wilcox_test(reformulate("Incubation_Time", response = "Ratio")) %>% 
  adjust_pvalue(method = "BH") %>%            # optional
  add_significance("p.adj") %>%               # uses p.adj.signif
  # y-position for labels
  left_join(df_plot %>% group_by(Embryonic_Age, Analogue) %>%
              summarise(y.position = max(Ratio, na.rm = TRUE), .groups="drop"),
            by = c("Embryonic_Age", "Analogue")) %>%
  mutate(Embryonic_Age = factor(Embryonic_Age, levels = c("E12.5","E13.5","E14.5","E15.5","E16.5"))) %>%
  mutate(y.position = y.position + 0.1) %>%
  # columns required by stat_pvalue_manual for 2-group comparisons:
  mutate(group1 = "4h", group2 = "8h")

p <- ggplot(df_plot, aes(x = Embryonic_Age, y = Ratio, fill = Incubation_Time)) + #
  geom_boxplot(position = position_dodge(0.8)) +
  theme_minimal() +
  labs(title = "Boxplot proteins", x = "Samples", y = "Ratio") +
   stat_pvalue_manual(
     stats,
     label = "p.adj.signif",   # or "p.signif" if you skip adjust_pvalue()
     xmin = "Embryonic_Age", xmax = "Embryonic_Age",
     y.position = "y.position",
      tip.length = 0
   ) + 
  facet_grid(~ Analogue, scales = "free_y") 

pdf(paste(path,"Figures/Fig1d_Boxplot_TimecourseSILAC_4h8h_HT&LT.pdf", sep = ""), width=10, height=6)
print(p)
dev.off()

rm(df_plot, stats, p)


# ---- Figure 1e ----
# quantified proteins per sample
data_wide <- data %>%
  dplyr::select(protein_group, Embryonic_Age, Incubation_Time, Replicates, LFQ_L, LFQ_H) %>%
  pivot_longer(cols = c(LFQ_L, LFQ_H), names_to = "Analogue", values_to = "LFQ", names_prefix = "LFQ_") %>%
  drop_na(LFQ) %>%
  pivot_wider(names_from= c("Embryonic_Age", "Incubation_Time", "Replicates", "Analogue"), values_from = LFQ, id_cols = "protein_group") 

# Number of proteins
non_na_counts <- colSums(!is.na(data_wide[, -1])) %>% as.data.frame()
names(non_na_counts) <- "Count"
non_na_counts$Sample <- rownames(non_na_counts)
non_na_counts <- non_na_counts %>%
    separate(Sample, into=c("Embryonic_Age", "Incubation_Time", "Replicates", "Analogue"), sep = "_") #"Experiment_Number", 

sum_counts <- non_na_counts %>%
  group_by(Embryonic_Age, Incubation_Time, Analogue) %>% #Experiment_Number
  summarise(
    mean_count = mean(Count, na.rm = TRUE),
    sd_count   = sd(Count, na.rm = TRUE),
    n          = sum(!is.na(Count)),
    se_count   = sd_count / sqrt(n),
    .groups = "drop"
  ) 

p <- ggplot(sum_counts, aes(x = reorder(Embryonic_Age, -mean_count), y = mean_count, fill=Analogue)) + #+, col=Experiment_Number
      geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
      geom_errorbar(
        aes(ymin = mean_count - sd_count, ymax = mean_count + sd_count),
        width = 0.25, stat = "identity", position = position_dodge(width = 0.9)
      ) +
      geom_text(aes(label = round(mean_count)), vjust = -0.5, stat = "identity",position = position_dodge(width = 0.9)) +
      labs(
        title = "Number of non-NA proteins per Sample (mean ± SD across replicates)",
        x = "Samples",
        y = "Number of proteins"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
            text = element_text(size = 14)) +
      facet_grid(~Incubation_Time)

pdf(paste(path,paste("Figures/Fig1e_NumberProt_TimecourseSILAC.pdf", sep = ""), sep = ""), width=10)
print(p)
dev.off()

rm(data_wide, sum_counts, p)


# ---- Extended figure 1e ----
# Arginine to proline conversion
arg_pro <- read_tsv("/project/KrausharLab/Cortex_team/Marta/Paper/Scripts/Github_structure/Data/SILAC_Timecourse_Bulk_MS/evidence_ArgProConversion.txt") %>%
  select(Sequence, Pro6, Type, `Raw file`, `Intensity H`, Reverse,
         `Potential contaminant`)

ev <- arg_pro %>%
  filter(is.na(Reverse), is.na(`Potential contaminant`)) %>%
  mutate(
    nP = str_count(Sequence, "P"),
    run = as.integer(str_extract(`Raw file`, "\\d+$")),
    intensity = replace_na(as.numeric(`Intensity H`), 0),
    Pro6 = replace_na(as.integer(Pro6), 0))
  # ) %>%
  # filter(run == 6) #run 6

n_identified_pro <- n_distinct(ev$Sequence)

# Remove isotope-envelope coincidences
ev <- ev %>%
  filter(!(Pro6 > 0 & str_starts(Type, "MULTI"))) %>%
  mutate(form = if_else(Pro6 == 0, "I_unconverted", "I_converted"))

peptides <- ev %>%
  filter(intensity > 0) %>%
  group_by(Sequence, nP, form, run) %>%
  summarise(intensity = sum(intensity), .groups = "drop") %>%
  pivot_wider(names_from = form, values_from = intensity, values_fill = 0)

# Add missing columns if necessary
if (!"I_unconverted" %in% names(peptides)) peptides$I_unconverted <- 0
if (!"I_converted" %in% names(peptides)) peptides$I_converted <- 0

observable_pro <- peptides %>%
  filter(nP > 0, I_unconverted > 0)

observed_pro <- observable_pro %>%
  filter(I_converted > 0) %>%
  mutate(
    r = I_converted / I_unconverted,
    C = 100 * r / (nP + r)
  )

n_observable_pro <- nrow(observable_pro)
n_observed_pro <- nrow(observed_pro)

cat("identified", n_identified_pro,
    ", observable", n_observable_pro,
    ", observed", n_observed_pro, "\n")
cat(sprintf("median C = %.2f%%  IQR %.2f-%.2f%%\n",
            median(observed_pro$C),
            quantile(observed_pro$C, 0.25),
            quantile(observed_pro$C, 0.75)))

# Add metadata
metadata <- read.xlsx("/project/KrausharLab/Cortex_team/Marta/Paper/Scripts/Github_structure/Data/SILAC_Timecourse_Bulk_MS/SILACtimecourse_samples_metadata.xlsx", sheet = 1)
observed_pro <- merge(observed_pro, metadata[, c("Sample.ID", "Embryonic.Age", "Replicates", "Incubation.Time")], by.x = "run", by.y = "Sample.ID")

# Plot Arg -> Pro
med <- median(C)

p <- ggplot(observed_pro, aes(x = C, col = Embryonic.Age, fill = Embryonic.Age)) +
  geom_density(alpha = 0.7) +
  labs(x = "Arg→Pro conversion (% per proline)", y = "Count") +
  theme_minimal(base_size = 8) +
  facet_grid(~Incubation.Time) +
  scale_color_manual(values = c("E12.5"="#00aeef", "E13.5"="#007169", "E14.5"="#2bb573", "E15.5"="#f6931d", "E16.5"="#eb008b")) +
  scale_fill_manual(values = c("E12.5" = "#00aeef","E13.5" = "#007169","E14.5" = "#2bb573","E15.5" = "#f6931d","E16.5" = "#eb008b")) + 
  # geom_vline(xintercept = med, linewidth = 0.4) +
  # annotate("text", x = med / 1.13, y = max(d_df$y) * 1.27,
  #          label = sprintf("median %.1f%%", med), hjust = 1, size = 2.6) +
    annotate("text", x = 15, y = 2,
           label = sprintf("%s peptides identified\n%s observable\n%s observed",
                           n_identified_lab, n_observable_lab, n_observed_lab),
           hjust = 1, vjust = 1, size = 4, color = "#6b6b6b") +
  xlim(0, 20) 

pdf(paste(path,paste("Figures/ExtFig1e_ArgToPro_conversion.pdf", sep = ""), sep = ""), width=10)
print(p)
dev.off()


# ---- Extended figure 1f ----
# Labeling efficiency estimation
lab_efficiency <- read_tsv("/project/KrausharLab/Cortex_team/Marta/Paper/Scripts/Github_structure/Data/SILAC_Timecourse_Bulk_MS/evidence_labellingEfficiency.txt") %>%
  select(Sequence, Arg10norm, Lys8norm, `Raw file`, Intensity, Reverse,
         `Potential contaminant`)

vm <- lab_efficiency %>%
  filter(is.na(Reverse), is.na(`Potential contaminant`)) %>%
  mutate(
    S = str_count(Sequence, "K") + str_count(Sequence, "R"),
    n_heavy = replace_na(as.integer(Arg10norm), 0) +
              replace_na(as.integer(Lys8norm), 0),
    intensity = replace_na(as.numeric(Intensity), 0),
    run = as.integer(str_extract(`Raw file`, "\\d+$"))
  ) 
  # filter(run == RUN)

n_identified_lab <- n_distinct(vm$Sequence)
vm <- vm %>% filter(intensity > 0)

per_peptide <- list()
n_observable_lab <- 0
n_observed_lab <- 0

for (S in sort(unique(vm$S[vm$S >= 2]))) {

  peptides <- vm %>%
    filter(.data$S == S) %>%
    group_by(Sequence, n_heavy, run) %>%
    summarise(intensity = sum(intensity), .groups = "drop") %>%
    pivot_wider(names_from = n_heavy, values_from = intensity, values_fill = 0)

  all_heavy <- as.character(S)
  one_light <- as.character(S - 1)

  if (!all_heavy %in% names(peptides)) next

  observable_lab <- peptides %>%
    filter(.data[[all_heavy]] > 0)

  n_observable_lab <- n_observable_lab + nrow(observable_lab)

  if (!one_light %in% names(peptides)) next

  observed_lab <- observable_lab %>%
    filter(.data[[one_light]] > 0)

  n_observed_lab <- n_observed_lab + nrow(observed_lab)

  r <- observed_lab[[one_light]] / observed_lab[[all_heavy]]
   per_peptide[[as.character(S)]] <- observed_lab %>%
    transmute(
      Sequence,
      run,
      labelling_efficiency = 100 * S / (S + r)
    )
}

a <- bind_rows(per_peptide, .id = "S")

cat("identified", n_identified_lab,
    ", observable", n_observable_lab,
    ", observed", n_observed_lab, "\n")
cat(sprintf("median a = %.2f%%  IQR %.2f-%.2f%%\n",
            median(a$labelling_efficiency), quantile(a$labelling_efficiency, 0.25), quantile(a$labelling_efficiency, 0.75)))

# Add metadata
metadata <- read.xlsx("/project/KrausharLab/Cortex_team/Marta/Paper/Scripts/Github_structure/Data/SILAC_Timecourse_Bulk_MS/SILACtimecourse_samples_metadata.xlsx", sheet = 1)
a <- merge(a, metadata[, c("Sample.ID", "Embryonic.Age", "Replicates", "Incubation.Time")], by.x = "run", by.y = "Sample.ID")

# Plot labelling efficiency
med <- median(a$labelling_efficiency)

p <- ggplot(a, aes(x = labelling_efficiency, col = Embryonic.Age, fill = Embryonic.Age)) +
  geom_density(alpha = 0.7) +
  labs(x = "Labelling efficiency (%)", y = "Count") +
  theme_minimal(base_size = 8) +
  scale_color_manual(values = c("E12.5"="#00aeef", "E13.5"="#007169", "E14.5"="#2bb573", "E15.5"="#f6931d", "E16.5"="#eb008b")) +
  scale_fill_manual(values = c("E12.5" = "#00aeef","E13.5" = "#007169","E14.5" = "#2bb573","E15.5" = "#f6931d","E16.5" = "#eb008b")) + 
  facet_grid(~Incubation.Time) +
  # geom_vline(xintercept = med, linewidth = 0.4) +
  # annotate("text", x = med / 1.13, y = max(d_df$y) * 1.27,
  #          label = sprintf("median %.1f%%", med), hjust = 1, size = 2.6) +
  annotate("text", x = 60, y = 0.12,
           label = sprintf("%s peptides identified\n%s observable\n%s observed",
                           n_identified_lab, n_observable_lab, n_observed_lab),
           hjust = 1, vjust = 1, size = 4, color = "#6b6b6b") +
  xlim(0, 100) 

pdf(paste(path,paste("Figures/ExtFig1f_Labelling_efficiency.pdf", sep = ""), sep = ""), width=10)
print(p)
dev.off()


# ---- Figure 2b ----
# PCA of log2 riBAQ values
data_plot <- data %>%
  ungroup() %>%
  dplyr::select(protein_group, Embryonic_Age, Incubation_Time, Replicates, riBAQ_H, riBAQ_L) %>%
  pivot_longer(cols = c(riBAQ_H, riBAQ_L), names_to = "Analogue", values_to = "riBAQ", names_prefix = "riBAQ_") %>%
  filter(!is.na(riBAQ)) %>%
  mutate(log2_riBAQ = log2(riBAQ)) %>%
  pivot_wider(names_from = c(Embryonic_Age, Incubation_Time, Replicates, Analogue), values_from = log2_riBAQ, id_cols = c("protein_group")) #Experiment_Number, // , "genes" 

data_plot <- data_plot %>%
  dplyr::select(-protein_group) %>% #, -genes
  t()

data_plot <- data_plot[, colSums(is.na(data_plot)) == 0, drop = FALSE]

# PCA 4H
pca_result_4h <- prcomp(data_plot[grepl("4h", rownames(data_plot)),], center = TRUE, scale. = TRUE)
pca_data_4h <- as.data.frame(pca_result_4h$x)[, 1:2] # Keep PC1 and PC2
pca_data_4h$Sample <- rownames(pca_data_4h)
pca_data_4h <- pca_data_4h %>%
  separate(Sample, into=c("Embryonic_Age", "Incubation_Time", "Replicates", "Analogue"), sep = "_") %>% #Experiment_Number, 
  mutate(
    Embryonic_Age = factor(Embryonic_Age, levels = c("E12.5","E13.5","E14.5","E15.5","E16.5")),
    Analogue = factor(Analogue, levels = c("H", "L"))
  )

# Calculate variance explained by each principal component
pca_var_4h <- pca_result_4h$sdev^2
pca_var_explained_4h <- round(100 * pca_var_4h / sum(pca_var_4h), 1)  # in percentage

# Create the PCA plot
p_4h <- ggplot(pca_data_4h, aes(x = PC1, y = PC2, colour = Embryonic_Age, shape=Analogue)) +
  geom_point(size = 5) +
  labs(
    title = paste("PCA of log2 riBAQ"),
    x = paste0("PC1 (", pca_var_explained_4h[1], "%)"),
    y = paste0("PC2 (", pca_var_explained_4h[2], "%)")
  ) +
  theme_minimal() +
  theme(text = element_text(size = 14))


# PCA 8H
pca_result_8h <- prcomp(data_plot[grepl("8h", rownames(data_plot)),], center = TRUE, scale. = TRUE)
pca_data_8h <- as.data.frame(pca_result_8h$x)[, 1:2] # Keep PC1 and PC2
pca_data_8h$Sample <- rownames(pca_data_8h)
pca_data_8h <- pca_data_8h %>%
  separate(Sample, into=c("Embryonic_Age", "Incubation_Time", "Replicates", "Analogue"), sep = "_") %>% #Experiment_Number, 
  mutate(
    Embryonic_Age = factor(Embryonic_Age, levels = c("E12.5","E13.5","E14.5","E15.5","E16.5")),
    Analogue = factor(Analogue, levels = c("H", "L"))
  )

# Calculate variance explained by each principal component
pca_var_8h <- pca_result_8h$sdev^2
pca_var_explained_8h <- round(100 * pca_var_8h / sum(pca_var_8h), 1)  # in percentage

# Create the PCA plot
p_8h <- ggplot(pca_data_8h, aes(x = PC1, y = PC2, colour = Embryonic_Age, shape=Analogue)) +
  geom_point(size = 5) +
  labs(
    title = paste("PCA of log2 riBAQ"),
    x = paste0("PC1 (", pca_var_explained_8h[1], "%)"),
    y = paste0("PC2 (", pca_var_explained_8h[2], "%)")
  ) +
  theme_minimal() +
  theme(text = element_text(size = 14))

pdf(paste(path,"Figures/Fig2b_PCA_TimecourseSILAC_Log2riBAQ.pdf", sep = ""))
print(p_4h | p_8h)
dev.off()

rm(data_plot, pca_result, pca_data, pca_var, pca_var_explained, p)


# ---- Extended Figure 2a ----
# effect of free-AA normalisation
df_clean <- data %>%
  filter(!is.na(ratios_HT), !is.na(ratios_HT_normFreeAA)) %>%
  ungroup()

# Create wide matrices for PCA (Proteins as columns, Samples as rows)
# We handle missing values by selecting proteins present in most samples or using simple imputation
# Here we keep proteins present in all samples for a fair Procrustes comparison
common_proteins <- df_clean %>%
  group_by(protein_group) %>%
  tally() %>%
  filter(n == max(n)) %>%
  pull(protein_group) 

pca_data_1 <- df_clean %>%
  filter(protein_group %in% common_proteins) %>%
  dplyr::select(sample, protein_group, ratios_HT) %>%
  pivot_wider(names_from = protein_group, values_from = ratios_HT) %>%
  column_to_rownames("sample")

pca_data_2 <- df_clean %>%
  filter(protein_group %in% common_proteins) %>%
  dplyr::select(sample, protein_group, ratios_HT_normFreeAA) %>%
  pivot_wider(names_from = protein_group, values_from = ratios_HT_normFreeAA) %>%
  column_to_rownames("sample")

# Run PCAs
pca1 <- prcomp(pca_data_1, scale. = TRUE)
pca2 <- prcomp(pca_data_2, scale. = TRUE)

# Procrustes Analysis: Rotating PCA2 to fit PCA1
pro_analysis <- procrustes(X = pca1, Y = pca2, symmetric = TRUE)
prot_test <- protest(X = pca1, Y = pca2, scores = "sites", permutations = 999) 

prot_test$t0 # Correlation between the two configurations (0 to 1, higher means more similar)
prot_test$signif # p-value from permutation test (lower means more significant similarity)
prot_test$ss # Sum of squared distances between the two configurations (lower means more similar)

# Plot Procrustes results
# X is the target (Normalization 1), Yrot is the rotated/scaled version (Normalization 2)
points_norm1 <- as.data.frame(pro_analysis$X[, 1:2])
points_norm2 <- as.data.frame(pro_analysis$Yrot[, 1:2])

colnames(points_norm1) <- c("PC1", "PC2")
colnames(points_norm2) <- c("PC1", "PC2")

# Add metadata back to the points
metadata <- df_clean %>%
  dplyr::select(sample, Embryonic_Age, Incubation_Time, Replicates) %>%
  distinct() %>%
  column_to_rownames("sample")

# Ensure the metadata matches the order of the PCA points
metadata <- metadata[rownames(points_norm1), ]

plot_data <- bind_rows(
  points_norm1 %>% mutate(Method = "HT", SampleID = rownames(points_norm1)),
  points_norm2 %>% mutate(Method = "HT_normFreeAA", SampleID = rownames(points_norm2))
) %>%
  mutate(
    Embryonic_Age = metadata$Embryonic_Age[match(SampleID, rownames(metadata))],
    Incubation_Time = metadata$Incubation_Time[match(SampleID, rownames(metadata))],
    Method = factor(Method, levels = c("HT", "HT_normFreeAA"))
  ) %>%
  arrange(SampleID, Method)

# Create the Plot
p <- ggplot(plot_data, aes(x = PC1, y = PC2, color = Embryonic_Age, shape = Incubation_Time)) +
  # Draw Arrows connecting the same sample (from HT to HT_normFreeAA)
  geom_path(aes(group = SampleID), 
            arrow = arrow(length = unit(0.20, "cm"), type = "closed"), 
            alpha = 0.4, 
            linetype = "solid") +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "Procrustes Analysis: HT vs HT_normFreeAA",
    subtitle = "Dashed lines show how much individual samples 'moved' due to normalization",
    x = "Procrustes Dim 1",
    y = "Procrustes Dim 2"
  ) +
  scale_color_brewer(palette = "Spectral")

pdf(paste(path,"Figures/ExtFig2a_Procrustes_PCA_HT_vs_HTnormFreeAA.pdf", sep = ""), width=10)
print(p)
dev.off()

rm(df_clean, common_proteins, pca_data_1, pca_data_2, pca1, pca2, pro_analysis, prot_test, points_norm1, points_norm2, metadata, plot_data, p)


# ---- Detect sample-level QC outliers ----
iqr_per_sample_nonNAcounts <- non_na_counts %>%
  group_by(Embryonic_Age, Incubation_Time) %>%
  summarise(Q1 = quantile(Count, 0.25, na.rm = TRUE),
    Q3 = quantile(Count, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    Lower_Limit = Q1 - 1.5 * IQR,
    Upper_Limit = Q3 + 1.5 * IQR, .groups = "drop")

outliers_nonNAcounts <- non_na_counts %>%
  tibble::rownames_to_column("ID") %>%
  left_join(iqr_per_sample_nonNAcounts %>% select(Embryonic_Age, Incubation_Time, Lower_Limit, Upper_Limit),
            by = c("Embryonic_Age", "Incubation_Time")) %>%
  filter(Count < Lower_Limit | Count > Upper_Limit) %>%
  pull(ID)

## Outlier sample IDs based on median fraction_H, by sample
median_fractionH_per_sample <- data %>%
  mutate(Sample = paste(Embryonic_Age, Incubation_Time, Replicates, sep = "_")) %>%
  group_by(Embryonic_Age, Incubation_Time, Replicates, Sample) %>%
  summarise(Median = median(ratios_HL, na.rm = TRUE), .groups = "drop")

iqr_fractionH_per_sample <- median_fractionH_per_sample %>%
  group_by(Embryonic_Age, Incubation_Time) %>%
  summarise(
    Q1 = quantile(Median, 0.25, na.rm = TRUE),
    Q3 = quantile(Median, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    Lower_Limit = Q1 - 1.5 * IQR,
    Upper_Limit = Q3 + 1.5 * IQR,
    .groups = "drop"
  )

# Flag sample as "outlier sample" if it contains at least one outlier fraction_H value
outliers_fractionH <- median_fractionH_per_sample %>%
  left_join(iqr_fractionH_per_sample, by = c("Embryonic_Age", "Incubation_Time")) %>%
  filter(!is.na(Median),
         Median < Lower_Limit | Median > Upper_Limit) %>%
  distinct(Sample) %>%
  pull(Sample)

## Optional: samples that are outliers by either criterion
outliers_any <- union(outliers_nonNAcounts, outliers_fractionH)

cat("outliers nonNA counts: ",outliers_nonNAcounts, "\n")
cat("outliers log2 ratios: ",outliers_fractionH, "\n")

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
