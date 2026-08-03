#########################################################################################################################################################################
# Free amino-acid MS processing and validation plots
#
# Purpose:
#   Clean amino-acid peak tables, validate transitions, normalise by valine/BCA, and export validation plots and processed RDS objects.
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
  paste0("01_Fig1b_FreeAA_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)

# Create the log file and write session information
cat(
  "Free amino acid analysis log\n",
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
  library(ggplot2)
  library(TOSTER)
  library(ggpubr)
  library(openxlsx)
  library(stringr)
  library(tidyr)
})


# ---- User settings ----
path <- "./Output/"
dir.create(path, recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(path, "Objects/"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(path, "Figures/"), recursive = TRUE, showWarnings = FALSE)

data_path <- "./Data/Free_AminoAcid_MS/"

# Load amino-acid peak table exported from the MS processing software.
results_0 = read.xlsx(paste0(data_path, "Free_AA_results.xlsx"))

# Standardise column names, convert decimal commas and parse precursor/fragment masses.
results = results_0
results$Component.Name <- gsub(" ", "", results$Component.Name)
names(results) <- c("Sample_Name", "Component_Name", "Area", "Height", "Retention_Time", "Width_at_50", "Signal_to_Noise")

results <- results %>%
    mutate(across(everything(), ~na_if(.x, "N/A"))) %>%
    mutate(Component_Name = gsub("(\\([^()]*)\\([^()]*\\)([^()]*\\))", "\\1\\2", Component_Name)) %>%
    extract(
        Component_Name,
        into = c("mass_AA", "mass_fragment"),
        regex = "\\(([^/]+)/(.*)\\)",
        remove = FALSE
    ) %>%
    mutate(
        mass_AA = as.numeric(sub(",", ".", mass_AA, fixed = TRUE)),
        mass_fragment = as.numeric(sub(",", ".", mass_fragment, fixed = TRUE)),
        Retention_Time = as.numeric(sub(",", ".", Retention_Time, fixed = TRUE)),
        Area = as.numeric(sub(",", ".", Area, fixed = TRUE)),
        Height = as.numeric(sub(",", ".", Height, fixed = TRUE)),
        Width_at_50 = as.numeric(sub(",", ".", Width_at_50, fixed = TRUE)),
        Signal_to_Noise = as.numeric(sub(",", ".", Signal_to_Noise, fixed = TRUE))
    ) %>%
    mutate(Component_Name = sub("\\s*\\(.*\\)\\s*$", "", Component_Name))

results$Sample_Name <- gsub("NEO", "Neo", results$Sample_Name)
results$Sample_Name <- gsub("SILAC_Ligation", "Ligation", results$Sample_Name)
results$Sample_Name <- gsub("Basal_ganglia", "BasalGanglia", results$Sample_Name)

# Change age that was incorrectly annotated
results$Sample_Name <- gsub("Midbrain_EO275_SILAC_1_rep3_E12-5", "Midbrain_EO275_SILAC_1_rep3_E16-5", results$Sample_Name)

# Add metadata columns
metadata <- read.xlsx(paste0(data_path, "BCA.xlsx"))

results <- results %>%
    separate(Sample_Name, into = c("Tissue", "Experiment_Number", "Treatment", "Number", "Technical_Replicate", "Age"), sep = "_", remove = FALSE)

results$Sample_Name <- paste(results$Tissue, results$Experiment_Number, results$Treatment, results$Number, results$Age, sep = "_")
results$Sample_Name <- gsub("_NA", "", results$Sample_Name)

results <- merge(results, metadata, by.x = "Sample_Name", by.y = "ID_MS", all.x = TRUE)

results <- results[,c("Sample_Name", "Tissue", "Experiment_Number", "Treatment", "Stage", "Duration", "Replicate", "Technical_Replicate", "Order_MS",
                     "Component_Name", "mass_AA", "mass_fragment", "Area", "Height", "Retention_Time", "Width_at_50", "Signal_to_Noise")]

names(results) <- c("Sample_Name", "Tissue", "Experiment_Number", "Treatment", "Embryonic_Age", "Incubation_Time", "Biological_Replicate", "Technical_Replicate", "Order_MS",
                     "Component_Name", "mass_AA", "mass_fragment", "Area", "Height", "Retention_Time", "Width_at_50", "Signal_to_Noise")

results$Sample_Name <- paste(results$Tissue, results$Experiment_Number, results$Embryonic_Age, results$Treatment, results$Biological_Replicate, results$Technical_Replicate, sep = "_")

results$Sample_Name <- gsub("Midbrain_EO277_NA_SILAC_NA_rep2", "Midbrain_EO277_E12.5_SILAC_1_rep2", results$Sample_Name)
results$Sample_Name <- gsub("Midbrain_EO277_NA_SILAC_NA_rep1", "Midbrain_EO277_E12.5_SILAC_1_rep2", results$Sample_Name)
results$Sample_Name <- gsub("_NA", "", results$Sample_Name)


# ---- Validate transition area ratios ----
# Compare measured transition ratios with the amino-acid standard.
results <- results[is.na(results$Area) == FALSE, ] # eliminate peaks with NAs

results <- results %>%
    group_by(Sample_Name, Component_Name) %>%
    mutate(Area_ratios = Area / Area[1])

validation_standard <- read.xlsx(paste0(data_path, "Validation_Standard.xlsx"))
names(validation_standard) <- validation_standard[1,]
validation_standard <- validation_standard[-1,]
validation_standard <- validation_standard[, c(1,5,12)]
validation_standard$Metabolite <- gsub(" ", "", validation_standard$Metabolite)
validation_standard$Metabolite <- gsub("-d8-", "-d8", validation_standard$Metabolite)
names(validation_standard) <- c("Component_Name", "mass_fragment", "Area_ratios_standard")
validation_standard <- validation_standard %>%
  mutate(
    mass_fragment = as.numeric(sub(",", ".", mass_fragment, fixed = TRUE)),
    Area_ratios_standard = as.numeric(sub(",", ".", Area_ratios_standard, fixed = TRUE))
  )

results <- merge(results, validation_standard, by = c("Component_Name", "mass_fragment"), all.x = TRUE)

results <- results %>%
  mutate(Diff_Area_ratios = Area_ratios_standard - Area_ratios,
        Pass_validation_Area_ratios = dplyr::between(Area_ratios_standard - Area_ratios, -0.1, 0.1))


# ---- Validate retention times ----
# Check that each transition elutes within the expected retention-time window.
results <- results %>%
  group_by(Component_Name, Sample_Name) %>%
  mutate(
    # TRUE if this row is within 3 (absolute) of at least one *other* row
    Pass_validation_retention_time = if (n() == 1) {
      FALSE
    } else {
      sapply(seq_along(Retention_Time), function(i) {
        any(abs(Retention_Time[i] - Retention_Time[-i]) < 3, na.rm = TRUE)
      })
    }
  ) %>%
  ungroup()


# ---- Select quantitative transition per amino acid ----
# Keep the first transition, excluding leucine (3rd) and isoleucine (2nd)
pick_transition <- function(.x, .y){
  comp <- tolower(.y$Component_Name)

  .x <- .x %>%
    mutate(
      mass_fragment = as.numeric(as.vector(mass_fragment)),
      mass_AA       = as.numeric(as.vector(mass_AA))
    )

  # drop pseudotransitions (mass_fragment == mass_AA) for ranking
  x_real <- .x %>% filter(mass_fragment != mass_AA)

  # if everything is pseudotransition, fall back to original group
  if (nrow(x_real) == 0) return(.x %>% dplyr::slice_max(mass_fragment, n = 1, with_ties = FALSE))

  k <- dplyr::case_when(
    comp == "leucine"     ~ 3L,
    comp == "isoleucine"  ~ 2L,
    TRUE                  ~ 1L
  )

  # if not enough real transitions exist, take the largest available
  k <- min(k, nrow(x_real))

  x_real %>%
    arrange(desc(mass_fragment)) %>%
    dplyr::slice(k)
}

results_firstTransition <- results %>%
  group_by(Sample_Name, Component_Name) %>%
  group_modify(~ pick_transition(.x, .y)) %>%
  ungroup()


# ---- Normalise by valine ----
# Use valine as the internal reference
normValine_results_firstTransition <- results_firstTransition %>%
  group_by(Sample_Name) %>%
  mutate(Area_norm_Valine = Area / Area[Component_Name == "Valine,2,3,4,4,4,5,5,5-d8"]) %>%
  ungroup()

rm(results_firstTransition)


# ---- Normalise by BCA ----
# Scale amino-acid signal by total protein content 
normValine_results_firstTransition <- merge(normValine_results_firstTransition, metadata[, c("Order_MS", "BCA")], by = "Order_MS", all.x = TRUE)
normValine_results_firstTransition <- normValine_results_firstTransition %>%
  mutate(Area_norm_Valine_BCA = Area_norm_Valine / BCA)

# Eliminate samples with BCA == 0 (it is 1 in the column BCA) or NA
normValine_results_firstTransition <- normValine_results_firstTransition[is.na(normValine_results_firstTransition$BCA) == FALSE & normValine_results_firstTransition$BCA != 1, ]

rm(metadata, results)


# ---- Remove poor technical replicates ----
# Exclusions follow the original manual QC decisions.
normValine_results_firstTransition_filtered <- normValine_results_firstTransition[!(normValine_results_firstTransition$Sample_Name %in% c("Heart_EO279_E12.5_Replenished_1_rep1", "Heart_EO279_E12.5_Replenished_1_rep2")),]

saveRDS(normValine_results_firstTransition_filtered, file.path(path, "Objects/normValine_results_firstTransition_filterTechReplicates.rds"))


# ---- Export plots ----
data_plot <- normValine_results_firstTransition_filtered %>%
    filter(Component_Name %in% c("L-Lysine", "L-Arginine", "L-Lysine13C615N2", "L-Arginine13C615N4")) %>%
    mutate(Analogue = ifelse(grepl("13C615N2|13C615N4", Component_Name), "heavy", "light"), 
            AminoAcid = str_remove(Component_Name, "13C615N2|13C615N4"))   


### Figure 1d
data_plot2 <-  data_plot %>%
  group_by(AminoAcid, Analogue, Tissue, Treatment, Embryonic_Age, Biological_Replicate) %>% #collapse technical replicates
  summarise(value = median(Area_norm_Valine_BCA, na.rm = TRUE), .groups = "drop") %>%
  ungroup() %>%
  group_by(AminoAcid, Tissue, Treatment, Embryonic_Age, Biological_Replicate) %>% #convert to percentage
  mutate(total = sum(value, na.rm = TRUE),
         pct   = if_else(total > 0, 100 * value / total, NA_real_)) %>%
  ungroup()

# Mean % and SE across biological replicates for each component in each condition
data_plot_mean <- data_plot2 %>%
  group_by(AminoAcid, Analogue, Tissue, Treatment, Embryonic_Age) %>%
  summarise(
    mean_pct = median(pct, na.rm = TRUE),
    sd_pct   = sd(pct, na.rm = TRUE),
    n        = sum(!is.na(pct)),
    se_pct   = sd_pct / sqrt(n),
    .groups  = "drop"
  )

# Plot
data_plot_mean$mean_pct <- data_plot_mean$mean_pct/100
data_plot_mean$sd_pct <- data_plot_mean$sd_pct/100
data_plot_mean <- data_plot_mean[data_plot_mean$Treatment != "Ligation",]

p <- ggplot(data_plot_mean, aes(fill=Analogue, y=mean_pct, x=Embryonic_Age)) + 
 geom_col(position = "fill") +
  geom_pointrange(
    data = data_plot_mean,
    aes(x = Embryonic_Age, y = mean_pct, ymin = pmax(0, mean_pct - sd_pct), ymax = pmin(100, mean_pct + sd_pct)), color = "black",
    position = position_dodge(width = 0.8),
    inherit.aes = FALSE
  ) +
  facet_wrap(AminoAcid+Treatment ~ Tissue) +
  labs(x = "Embryonic age | Treatment", y = "Composition (%)", fill = "Component") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0(path, "Figures/Fig1b_Percentage_HL_perSample_NormValineBCA_filtered_noLigation_Median.pdf"), width = 15, height = 20)
print(p)
dev.off()

t <- data_plot_mean %>%
  group_by(Analogue, Treatment, Embryonic_Age) %>%
  summarise(mean_pct = mean(mean_pct, na.rm = TRUE), .groups = "drop")

### Extended figure 1c
data_plot_medianE12 <- data_plot %>%
  filter(Embryonic_Age=="E12.5") %>%
  group_by(AminoAcid, Tissue, Treatment, Analogue) %>%
  summarise(median = median(Area_norm_Valine_BCA, na.rm = TRUE))

data_plot2 <- merge(data_plot, data_plot_medianE12, by=c("AminoAcid", "Tissue", "Treatment", "Analogue"))

data_plot2$Area_norm_Valine_BCA_normE12 <- data_plot2$Area_norm_Valine_BCA/data_plot2$median

p <- ggplot(data_plot2[data_plot2$Treatment != "Ligation", ],
    aes(x = Analogue, y = log2(Area_norm_Valine_BCA_normE12), fill = Embryonic_Age)) +
    geom_boxplot() +
    labs(x = "Samples", y = "Log2FC of E16.5 vs E12.5 Area normalized by Valine and BCA", title="Log2FC of E16.5 vs E12.5 Area per aminoacid per sample normalized by Valine and BCA") +
    theme(axis.text.x = element_text(size = 5)) +
    theme_minimal() +
    facet_grid(AminoAcid~Treatment+Tissue)

pdf(paste0(path, "Figures/ExtFig1c_Log2FC_E16vsE12_L_perSample_NormValineBCA_noLigation_WithE12.pdf"), width = 15)
print(p)
dev.off()


### Extended figure 1d
# data_plot_L <- data_plot[data_plot$Analogue=="light",]
data_plot2 <- data_plot %>%
  group_by(Tissue, Treatment, AminoAcid, Analogue) %>%
  mutate(
    median_Area_norm_Valine_BCA_E12 = median(
      Area_norm_Valine_BCA[Embryonic_Age == "E12.5"],
      na.rm = TRUE
    ),
    log2FC_E16vsE12 = if_else(
      Embryonic_Age == "E16.5" & !is.na(median_Area_norm_Valine_BCA_E12) & median_Area_norm_Valine_BCA_E12 > 0,
      log2(Area_norm_Valine_BCA / median_Area_norm_Valine_BCA_E12),
      NA_real_
    )
  ) %>%
  ungroup()

data_wide <- data_plot2 %>%
  filter(Embryonic_Age == "E16.5" & Treatment == "SILAC") %>%
  dplyr::select(AminoAcid, Tissue, Biological_Replicate, Technical_Replicate, Analogue, log2FC_E16vsE12) %>%
  distinct() %>%
  pivot_wider(
    names_from = Analogue, 
    values_from = log2FC_E16vsE12
  ) %>%
  # Remove rows where either light or heavy is NA to avoid plotting errors
  filter(!is.na(light) & !is.na(heavy))

p <- ggplot(data_wide, aes(x = light, y = heavy)) +
  geom_point(aes(col = AminoAcid), alpha = 0.6) + # Points colored by Treatment
  geom_smooth(method = "lm", color = "black", se = TRUE) + # Linear regression line
  stat_cor(method = "pearson", label.x = min(data_wide$heavy), label.y = max(data_wide$light)) + # Adds R and p-value
  facet_wrap(~Tissue) + # Creates separate plots "per condition" (Treatment)
  theme_minimal() +
  labs(
    title = "Correlation of log2FC (E16 vs E12) between Analogues",
    x = "log2FC (Light Analogue)",
    y = "log2FC (Heavy Analogue)",
    color = "Amino Acid"
  )

pdf(paste0(path, "Figures/ExtFig1d_Scatterplot_Log2FC_E16vsE12_HvsL_perSample_SILAC_NormValineBCA_noLigation.pdf"))
print(p)
dev.off()


# ---- Interpolate labelling percentages across ages ----
age_to_num <- function(x) as.numeric(str_remove(x, "^E"))

target_ages <- c("E13.5", "E14.5", "E15.5")
target_t    <- age_to_num(target_ages)

data_plot_interp <- data_plot_mean %>%
  mutate(t = age_to_num(Embryonic_Age)) %>%
  # keep only ages that define the interpolation endpoints
  filter(Embryonic_Age %in% c("E12.5", "E16.5")) %>%
  # define "series" keys: adjust as needed
  group_by(Embryonic_Age, Analogue, Treatment, Tissue, AminoAcid) %>%
  mutate(pct_mean_analogues = mean(mean_pct, na.rm = TRUE), .groups = "drop") %>%
  dplyr::select(-c(sd_pct, n, se_pct, mean_pct, AminoAcid)) %>%
  unique() %>%
  group_by(Analogue, Treatment, Tissue) %>%
  summarise(
    y12_5 = pct_mean_analogues[t == 12.5][1],
    y16_5 = pct_mean_analogues[t == 16.5][1],
    .groups = "drop"
  ) %>%
  # guardrails (log needs positive values)
  filter(is.finite(y12_5), is.finite(y16_5), y12_5 > 0, y16_5 > 0) %>%
  # compute per-series exponential rate k
  mutate(k = (log(y16_5) - log(y12_5)) / (16.5 - 12.5)) %>%
  # generate the intermediate ages
  tidyr::expand_grid(Embryonic_Age = target_ages) %>%
  mutate(
    t = age_to_num(Embryonic_Age),
    pct_interpolated = y12_5 * exp(k * (t - 12.5))
  )

data_plot_interp_sum <- data_plot_interp %>%
  filter(Analogue=="heavy" & Tissue == "Neocortex" & Treatment == "SILAC") %>%
  group_by(Embryonic_Age) %>%
  summarise(mean_pct = mean(pct_interpolated, na.rm = TRUE), 
            mean_y12_5 = mean(y12_5, na.rm = TRUE),
            mean_y16_5 = mean(y16_5, na.rm = TRUE), .groups = "drop")


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
