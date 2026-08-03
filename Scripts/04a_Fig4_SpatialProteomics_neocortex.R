#########################################################################################################################################################################
# Spatial proteomics preprocessing and DEP
#
# Purpose:
#   Preprocess spatial proteomics tables, create QC/spatial plots, run proDA differential enrichment, and perform GO enrichment.
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
  paste0("04_Fig4_SpatialProteomics_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
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
  library(RColorBrewer)
  library(circlize)
  library(clusterProfiler)
  library(corrplot)
  library(dplyr)
  library(ggplot2)
  library(ggpointdensity)
  library(ggrepel)
  library(gridExtra)
  library(org.Mm.eg.db)
  library(proDA)
  library(readxl)
  library(sf)
  library(stringr)
  library(tidyverse)
  library(viridis)
})


# ---- User settings ----
path <- "./Output/"
data_path <- "./Data/SILAC_SpatialProteomics_MS/"
dir.create(path, recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(path, "Objects/"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(path, "Figures/"), recursive = TRUE, showWarnings = FALSE)


# ---- Load StackedLFQ output tables ----
prot_ratios <- read.table(file.path(data_path, "SILAC_Spatial_ratios.csv"), sep = ",",header = TRUE)
prot_total <- read.table(file.path(data_path, "SILAC_Spatial_total.csv"), sep = ",",header = TRUE)
prot_light <- read.table(file.path(data_path, "SILAC_Spatial_light.csv"), sep = ",",header = TRUE)
prot_pulse <- read.table(file.path(data_path, "SILAC_Spatial_pulse.csv"), sep = ",",header = TRUE)

preprocess_sample <- function(df, analogue_name) {
  df <- df[, -1]
  df <- df %>% pivot_longer(cols = -c(genes, protein_group), names_to = "Sample", values_to = "Intensity")
  
  df[["Analogue"]] <- rep(analogue_name, nrow(df))
  
  df[["Sample"]] <- gsub("basal_ganglia", "basalganglia", df[["Sample"]])
  df[["Sample"]] <- unlist(lapply(str_split(df[["Sample"]], "_"), function(x) paste(x[5:7], collapse="_")))
  
  df <- df %>%
  separate(Sample, into=c("Area", "Tissue", "Replicate"), sep = "_", remove = FALSE)

  df$ID <- paste(df$Sample, df$Analogue, sep = "_")
  df$Area <- factor(df$Area, levels = c("1k", "2.5k", "5k", "7.5k", "10k", "15k", "20k", "30k"))

  return(df)
}

ratios <- preprocess_sample(prot_ratios, "ratios")
total <- preprocess_sample(prot_total, "total")
H <- preprocess_sample(prot_pulse, "heavy")
L <- preprocess_sample(prot_light, "light")

spatial_prot <- rbind(ratios, total, L, H)

# Add anatomical region metadata for each microdissected area
areas <- read_excel(file.path(data_path, "SILAC_Spatial_samples_metadata.xlsx"))
names(areas)[3] <- "Area_brain" 
areas$Sample <- gsub("BG", "basalganglia", areas$Sample)
areas$Sample <- gsub("ctx", "cortex", areas$Sample)
spatial_prot <- merge(spatial_prot, areas[, -1], by="Sample", all.x = TRUE)

spatial_prot[spatial_prot$Sample %in% c("15k_cortex_4","15k_cortex_3"), "Area_brain"] <- "CP" # Manual annotation based on image
spatial_prot[spatial_prot$Sample %in% c("15k_cortex_2","15k_cortex_1", "20k_cortex_4"), "Area_brain"] <- "VZ" # Manual annotation based on image

# Filter out outlier samples
spatial_prot <- spatial_prot %>%
  filter(!Sample %in% c("10k_basalganglia_3",
  "15k_basalganglia_3",
  "1k_cortex_4",
  "2.5k_cortex_5",
  "7.5k_cortex_3",
  "10k_cortex_5",
  "20k_cortex_1")) #"7.5k_cortex_3", "1k_cortex_4"


# ---- Map intensities to spatial coordinates ----
spatial_coordinates = st_read(file.path(data_path, "Geometry_coordinates.geojson"))

# Add plotting coordinates/colours for spatial visualisation.
spatial_coordinates <- spatial_coordinates %>%
  mutate(
    name = str_extract(classification, '(?<="name": ")[^"]+'),
    color = str_extract(classification, '(?<="color": \\[ )[\\d, ]+(?= ])')
  )

rgb_color_map <- data.frame(
  rgb = c("0, 0, 255", "255, 0, 0", "255, 0, 255", "0, 255, 0", 
          "255, 255, 0", "0, 255, 255", "255, 255, 255"),
  color_name = c("dark blue", "red", "pink", "green", 
                 "yellow", "light blue", "white") # Example names
)

spatial_coordinates <- spatial_coordinates %>%
  left_join(rgb_color_map, by = c("color" = "rgb")) %>%
  drop_na(objectType)

spatial_coordinates$name <- gsub("basal_ganglia", "basalganglia", spatial_coordinates$name)

rm(rgb_color_map)



# ---- Figure 4b ----
# Labelled-protein percentage by region 
data_wider <- spatial_prot %>%
  pivot_wider(id_cols=protein_group, names_from = ID, values_from = Intensity) 

# Count non-NA values in each column
non_na_counts <- colSums(!is.na(data_wider[, -1]))
non_na_counts <- as.data.frame(non_na_counts)
names(non_na_counts) <- "Count"

non_na_counts <- non_na_counts %>%
  rownames_to_column(var = "Sample") %>% 
  separate(Sample, into=c("Area", "Tissue", "Replicate", "Analogue"), sep = "_", remove = FALSE)

# Calculate % labelled protein per area
non_na_counts_percentage <- non_na_counts %>%
  group_by(Area, Tissue, Replicate) %>%
  mutate(
    total_proteins = sum(Count),
    percentage = (Count / total_proteins) * 100
  ) %>%
  ungroup()

# Heatmap of the % labelled proteins per area
non_na_counts_percentage$Sample <- gsub("_heavy|_light", "", non_na_counts_percentage$Sample)

merged_data_counts <- non_na_counts_percentage %>%
  left_join(spatial_coordinates, by = c("Sample" = "name")) %>%
  drop_na(Count)

data_plot <- merged_data_counts[merged_data_counts$Analogue=="heavy",]

p <- ggplot(data_plot) +
  geom_sf(aes(fill = percentage, geometry = geometry), color="black", size=0.5) +
  scale_fill_gradient2(low="white", high="red4", na.value="grey50", limits=c(0, max(data_plot$percentage))) +
  theme_minimal() +
  theme(legend.position="bottom",
        legend.text = element_text(angle=45, vjust=0.5, hjust=0.5)) +
  labs(title="Percentage of number of labelled proteins per area", fill="Percentage of Number of Labelled Protein")

pdf(paste(path, "Figures/Fig4b_Heatmap_PercentageNumberLabelledProteins.pdf", sep = ""))
print(p)
dev.off()



# ---- Figure 4c ----
# Number of non-NA proteins per sample with error bars
summary_counts <- non_na_counts %>%
  filter(Analogue %in% c("heavy", "light")) %>%
  group_by(Area, Tissue, Analogue) %>%
  summarise(
    mean_count = mean(Count, na.rm = TRUE),
    sd_count = sd(Count, na.rm = TRUE),
    .groups = "drop"
  )

summary_counts$Area <- factor(summary_counts$Area, levels = c("1k", "2.5k", "5k", "7.5k", "10k", "15k", "20k", "30k"))

p <- ggplot(summary_counts, aes(x = Area, y = mean_count, col = Analogue, fill = Analogue)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
      geom_errorbar(
        aes(ymin = mean_count - sd_count, ymax = mean_count + sd_count),
        position = position_dodge(width = 0.9),
        width = 0.3) +
      geom_point(data = non_na_counts, aes(x = Area, y = Count), color = "black",
                position = position_jitterdodge(jitter.width = 0, dodge.width = 0.9),
                size = 1.8, alpha = 0.7, inherit.aes = TRUE) + 
      geom_text(aes(label = round(mean_count)), position = position_dodge(width = 0.9), vjust = -0.5) +
      scale_color_manual(values = c("heavy" = "#b70309", "light" = "grey50")) +
      scale_fill_manual(values = c("heavy" = "#b70309", "light" = "grey50")) +
      labs(
        title = "Number of non-NA proteins per Sample",
        x = "Area",
        y = "Mean number of proteins ± SD") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
        text = element_text(size = 14)) +
      facet_wrap(~Tissue, ncol = 2)

# Create the bar plot using ggplot2 (per Area)
pdf(paste(path,paste("Figures/Fig4c_NumberProt_ErrorBars.pdf", sep = ""), sep = ""), width=10)
print(p)
dev.off()



# ---- Figure 4d ----
# Differential enrichment analysis in heavy between VZ and CP regions using proDA
data_deg <- spatial_prot %>%
            filter(Area_brain %in% c("VZ", "CP")) #& Analogue == "heavy"

# Add replicates of the VZ and CP
df <- unique(data_deg[, c("Area", "Replicate", "Area_brain")])
df <- df[order(df$Area_brain),]
df <- drop_na(df)
df <- df %>%
  group_by(Area_brain) %>%
  mutate(Replicate_area = row_number()) %>%
  ungroup()

data_deg <- merge(data_deg, df, by=c("Area", "Replicate", "Area_brain"), all.x = TRUE)

data_deg <- data_deg %>%
            mutate(ID_AreaBrain = paste(Area_brain, Replicate_area, Analogue, sep = "_"))

data_deg <- drop_na(data_deg) 

data_wide <- data_deg %>%
    pivot_wider(id_cols = protein_group, 
                names_from = ID_AreaBrain, 
                values_from = Intensity)

data_wide[data_wide==0] <- NA

#Gene names
genes_names <- data_wide[,1]

#Log2 transform the data for proDA analysis
log2 <- log2(data_wide[, -1])

#Create dataframe with the samples and additional information (condition and replicates)
sample_info_df <- data.frame(name = colnames(log2),stringsAsFactors = FALSE)
a <- str_split(sample_info_df$name, "_")
for (i in 1:length(a)){
  sample_info_df$condition[i] <- paste(a[[i]][1], a[[i]][3], sep = "_")
  sample_info_df$replicate[i] <- a[[i]][2]} 

#Fit the model with proDA
log2 <- as.matrix(log2)

fit_path <- file.path(path, "Objects/fit_spatial_heavy_VZvsCP.rds")

# if (file.exists(fit_path)) {
#   fit <- readRDS(fit_path)
# } else {
  fit <- proDA(
    log2,
    design = sample_info_df$condition,
    col_data = sample_info_df
    # reference_level = "Control"
  )

  saveRDS(fit, fit_path)
# }

rm(a, sample_info_df, log2, df)

# Run pairwise proDA contrasts between VZ and CP regions.
#result_names(fit) #to check the names of the conditions
test_res <- test_diff(fit, VZ_heavy - CP_heavy) 

volcano_data <- cbind(Genes=genes_names, log2FoldChange=test_res$diff, padj=test_res$adj_pval)
volcano_data <- as.data.frame(volcano_data)
volcano_data$log2FoldChange <- as.double(volcano_data$log2FoldChange)
volcano_data$padj <- as.double(volcano_data$padj)
Genes <- bitr(volcano_data$protein_group, fromType = "UNIPROT", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
volcano_data <- merge(volcano_data, Genes, by.x="protein_group", by.y="UNIPROT", all.x = TRUE)

# add a column of NAs
volcano_data$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
volcano_data$diffexpressed[volcano_data$log2FoldChange > 0.6 & volcano_data$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
volcano_data$diffexpressed[volcano_data$log2FoldChange<(-0.6) & volcano_data$padj < 0.05] <- "DOWN"

# Test GO enrichment among proteins classified as UP/DOWN.
# Background genes = all genes in this age
bg_genes <- volcano_data$protein_group
bg_entrez <- bitr(bg_genes, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

for (direction in c("UP", "DOWN")) {

  message("\nProcessing GO for ", direction, " regulated proteins")

  de_genes <- volcano_data %>%
    filter(diffexpressed == direction) %>%
    pull(protein_group) %>%
    unique() %>%
    na.omit()

  message("Proteins selected: ", length(de_genes))

  de_entrez <- bitr(de_genes, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

  message("Proteins mapped to ENTREZID: ", nrow(de_entrez))

  if (nrow(de_entrez) == 0) {
    message("No ENTREZ IDs found for ", direction)
    next
  }

  for (ont in c("BP", "MF", "CC")) {

    message("  Running ontology: ", ont)

    ego <- enrichGO(
      gene = unique(de_entrez$ENTREZID),
      universe = unique(bg_entrez$ENTREZID),
      OrgDb = org.Mm.eg.db,
      keyType = "ENTREZID",
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )

    ego_df <- as.data.frame(ego)

    message("  Enriched terms: ", nrow(ego_df))

    if (nrow(ego_df) == 0) {
      message("  No significant ", ont, " terms for ", direction)
      next
    }

    csv_file <- file.path(path, paste0("Objects/Fig4d_GO_", direction, "_", ont, "_VZvsCP_heavy.csv"))

    plot_file <- file.path(path, paste0( "Figures/Fig4d_GO_", direction, "_", ont, "_VZvsCP_heavy.pdf" ))

    write.csv(ego_df, file = csv_file, row.names = FALSE)

    p <- dotplot(ego,showCategory = min(10, nrow(ego_df))) +
      ggtitle(paste(ont,"—",direction,"regulated proteins, VZ vs CP heavy"))

    ggsave(filename = plot_file, plot = p, width = 8, height = 6)

    message("  CSV saved to: ", csv_file)
    message("  Plot saved to: ", plot_file)
  }
}

rm(plot_file, ego, de_genes, de_entrez, bg_genes, bg_entrez, direction, ont, p, volcano_data, test_res)



# ---- Figure 4e ----
# Heatmap of selected proteins in the VZ and CP regions
# Calculate H/T and H/L
spatial_prot2 <- spatial_prot %>%
  pivot_wider(
    id_cols = c(protein_group, genes, Sample, Area, Tissue, Area_brain, Replicate),
    names_from = Analogue,
    values_from = Intensity
  ) %>%
  mutate(
    heavyTotal = if_else(!is.na(total) & total != 0, heavy / total, NA_real_)
  ) %>%
  pivot_longer(
    cols = c(total, heavy, light, ratios, heavyTotal),
    names_to = "Analogue",
    values_to = "Intensity"
  ) %>%
  mutate(
    ID = paste(Sample, Analogue, sep = "_")
  )


# Median of all proteins per area, tissue, and analogue
all_prot <- spatial_prot2 %>%
              group_by(ID) %>%
              summarise(Intensity_median = median(Intensity, na.rm = TRUE), .groups = "drop") %>%
              separate(ID, into=c("Area", "Tissue", "Replicate", "Analogue"), sep = "_", remove = FALSE) %>%
              group_by(Analogue, Tissue) %>%
              mutate(
                Intensity_median_scaled_AreaBrainAnalogue = Intensity_median / max(Intensity_median, na.rm = TRUE), 
                Sample = paste(Area, Tissue, Replicate, sep = "_")) %>%
              ungroup() %>%
              left_join(spatial_coordinates, by = c("Sample" = "name"))

p <- ggplot(all_prot) +
  geom_sf(aes(fill = Intensity_median_scaled_AreaBrainAnalogue, geometry = geometry), color="black", size=0.5) +
  scale_fill_gradient2(low="white", high="red4", na.value="grey50", limits=c(0, max(all_prot$Intensity_median_scaled_AreaBrainAnalogue))) +
  theme_minimal() +
  theme(legend.position="bottom",
        legend.text = element_text(angle=45, vjust=0.5, hjust=0.5)) +
  labs(title="Median Protein Intensity Overlay", fill="Median LFQ intensity/ratio") +
  facet_grid(~Analogue)

pdf(paste(path, "Figures/Fig4e_Heatmap_MedianAllProteinIntensityRatios.pdf", sep = ""))
print(p)
dev.off()


# Selected protein per area, tissue, and analogue
protein = c("Neurod2", "Dnmt1", "Rbm3", "Ncam1")

protein_spatial <- spatial_prot2 %>%
              filter(genes %in% protein) %>%
              group_by(Analogue, Tissue, genes) %>%
              mutate(
                Intensity_median_scaled_AreaBrainAnalogue = Intensity / max(Intensity, na.rm = TRUE)) %>%
              ungroup() %>%
              left_join(spatial_coordinates, by = c("Sample" = "name"))

p <- ggplot(protein_spatial) +
  geom_sf(aes(fill = Intensity_median_scaled_AreaBrainAnalogue, geometry = geometry), color="black", size=0.5) +
  scale_fill_gradient2(low="white", high="red4", na.value="grey50", limits=c(0, max(protein_spatial$Intensity_median_scaled_AreaBrainAnalogue))) +
  theme_minimal() +
  theme(legend.position="bottom",
        legend.text = element_text(angle=45, vjust=0.5, hjust=0.5)) +
  labs(title="Median Protein Intensity Overlay", fill="Median LFQ intensity/ratio") +
  facet_grid(genes~Analogue)

pdf(paste(path, "Figures/Fig4e_Heatmap_", protein, "_IntensityRatios_2.pdf", sep = ""))
print(p)
dev.off()



# ---- Figure 4f ----
# Proteins halflife analysis
all_data <- readRDS(file.path(path, "Objects/all_data.rds"))

t12 <- all_data %>%
  filter(Experiment == "t12" & Embryonic_Age == "E14.5" & Incubation_Time == "8h") %>%
  unique()

# Extract the proteins and calculate the density height at each Value
proteins_to_mark <- c("Neurod2", "Dnmt1", "Rbm3", "Ncam1")

marked_proteins <- t12 %>%
  filter(genes %in% proteins_to_mark) %>%
  distinct(genes, Value) %>%
  mutate(
    density_y = approx(
      x = density(t12$Value, na.rm = TRUE)$x,
      y = density(t12$Value, na.rm = TRUE)$y,
      xout = Value
    )$y
  )

median_value <- median(t12$Value, na.rm = TRUE)

p <- ggplot(t12, aes(x = Value)) +
      geom_density(fill = "lightblue", alpha = 0.5) +

      # Median
      geom_vline(
        xintercept = median_value,
        color = "red",
        linetype = "dashed"
      ) +
      annotate(
        "text",
        x = median_value,
        y = 0,
        label = paste("Median:", round(median_value, 2)),
        vjust = -1,
        hjust = 1.1,
        color = "red"
      ) +

      # Selected proteins
      geom_vline(
        data = marked_proteins,
        aes(xintercept = Value),
        linetype = "dotted",
        color = "black"
      ) +
      geom_point(
        data = marked_proteins,
        aes(x = Value, y = density_y),
        size = 2.5
      ) +
      geom_text_repel(
        data = marked_proteins,
        aes(
          x = Value,
          y = density_y,
          label = paste0(genes, "\n", round(Value, 2), " h")
        ),
        direction = "y",
        nudge_y = 0.01,
        min.segment.length = 0,
        show.legend = FALSE
      ) +

      labs(
        title = "Density plot of protein half-lives",
        x = "Half-life (hours)",
        y = "Density"
      ) +
      theme_minimal() +
      xlim(0,500)

pdf(paste(path, "Figures/Fig4f_DensityPlot_halfLives_MarkedProteins.pdf", sep = ""))
print(p)
dev.off()


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
