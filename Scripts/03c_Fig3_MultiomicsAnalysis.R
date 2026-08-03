#########################################################################################################################################################################
# Multiomics integration: RNAseq, Riboseq and SILAC proteomics
#
# Purpose:
#   Clustering and GO enrichment for Multiomics.
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
  paste0("03c_Fig3_MultiomicsAnalysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)

# Create the log file and write session information
cat(
  "Multiomics analysis log\n",
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


k = 11 # Multiomics clusters after elbow plot

all_data <- readRDS(file.path(path, "Objects/all_data.rds"))
mat <- readRDS(file.path(path, "Objects/Zscore_RPF_RNAseq_riBAQH_T_HL_4h.rds"))
sd <- readRDS(file.path(path, "Objects/Zscore_RPF_RNAseq_riBAQH_T_HL_4h_SD.rds"))


# ---- Figure 3a and 3b ----
# Boxplots for total, H/T, H/L and t12
data_plot <- all_data[all_data$Experiment=="ratios_HT",]

p <- ggplot(data_plot, aes(x = Embryonic_Age, y = log2(Value))) + 
    geom_boxplot() +
    facet_grid(~Incubation_Time)

pdf(paste0(path, "Figures/Fig3a_Boxplot_Log2RatiosHT.pdf"))
print(p)
dev.off()


data_plot <- all_data[all_data$Experiment=="ratios_HL",]

p <- ggplot(data_plot, aes(x = Embryonic_Age, y = log2(Value))) + 
    geom_boxplot() +
    facet_grid(~Incubation_Time)

pdf(paste0(path, "Figures/Fig3a_Boxplot_Log2RatiosHL.pdf"))
print(p)
dev.off()


data_plot <- all_data[all_data$Experiment=="LFQ_T",]

p <- ggplot(data_plot, aes(x = Embryonic_Age, y = log2(Value))) + 
    geom_boxplot() +
    facet_grid(~Incubation_Time)

pdf(paste0(path, "Figures/Fig3a_Boxplot_Log2LFQtotal.pdf"))
print(p)
dev.off()


data_plot <- all_data[all_data$Experiment=="t12",]

p <- ggplot(data_plot, aes(x = Embryonic_Age, y = log2(Value))) + 
    geom_boxplot()

pdf(paste0(path, "Figures/Fig3b_Boxplot_Log2t12.pdf"))
print(p)
dev.off()



# ---- Elbow plot ----
k_max <- 30                  # elbow curve max k

k_max <- min(k_max, nrow(mat))  # safety
wss <- numeric(k_max)

for (kk in 1:k_max) {
  set.seed(123)               # reproducibility
  km_tmp <- kmeans(mat, centers = kk, nstart = 50)
  wss[kk] <- km_tmp$tot.withinss
}

elbow_df <- data.frame(k = 1:k_max, tot_withinss = wss)

p_elbow <- ggplot(elbow_df, aes(k, tot_withinss)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = 1:k_max) +
  labs(
    title = "Elbow plot (k-means on row-zscored log2 RPF, RNAseq, t12, LFQ T, H/L & 4h riBAQH)",
    x = "k",
    y = "Total within-cluster sum of squares (WSS)"
  ) +
  theme_minimal(base_size = 11)

pdf(paste0(path, "Figures/ElbowPlot_log2HL_log2T_log2RPF_log2RNAseq_log2riBAQH4h.pdf", sep=""))
print(p_elbow)
dev.off()


# ---- Extended Figure 3a ----
# Multiomics heatmap per cluster
col_df <- data.frame(
  assay = sub("_.*$", "", colnames(mat)),                
  Embryonic_Age   = sub("^.*_E", "E", colnames(mat)),               
  stringsAsFactors = FALSE
)

col_df$age <- factor(col_df$Embryonic_Age)
col_df$assay <- factor(col_df$assay)

mat <- mat[, order(col_df$assay, col_df$age), drop = FALSE]
col_df <- col_df[order(col_df$assay, col_df$age), , drop = FALSE]

set.seed(123)               # reproducibility

kmeans_path <- file.path(paste0(path,"/Objects/Fig3a_fixed_kmeans_k", k, ".rds"))

if (file.exists(kmeans_path)) {
  km <- readRDS(kmeans_path)
} else {
  km <- kmeans(mat, centers = k, nstart = 50, iter.max = 100)

  saveRDS(km, kmeans_path)
}

clusters <- km$cluster
centroids <- km$centers

# Check that the current matrix contains exactly the same proteins
if (!setequal(rownames(mat), names(clusters))) {
  stop(
    "The proteins in mat do not match the proteins used for the ",
    "saved clustering."
  )
}

# Put mat_wide in exactly the same protein order as the saved clustering
mat <- mat[names(clusters), , drop = FALSE]

if (!setequal(rownames(mat), names(clusters))) {
  stop(
    "The proteins in mat_wide do not match the proteins used for the ",
    "published clustering."
  )
}

# order within each cluster by distance to centroid (cleaner blocks)
dist_to_centroid <- vapply(seq_len(nrow(mat)), function(i) {
cl <- clusters[i]
sum((mat[i, ] - centroids[cl, ])^2)
}, numeric(1))

ord <- order(clusters, dist_to_centroid)
mat_z_ord <- mat[ord, , drop = FALSE]

row_ha <- rowAnnotation(
KMeans = factor(clusters[ord], levels = sort(unique(clusters))),
annotation_name_side = "top"
)

# optional: split heatmap by kmeans cluster (nice for readability)
row_split <- factor(clusters[ord], levels = sort(unique(clusters)))

ht <- Heatmap(
mat_z_ord,
name = "Row z-score",
left_annotation = row_ha,
row_split = row_split,               # show blocks per cluster
cluster_rows = FALSE,                # we already ordered by kmeans
cluster_columns = FALSE,             # keep Embryonic_Age order
show_row_names = FALSE,
column_split = col_df$assay,
column_title = "Embryonic Age",
row_title = paste0("Proteins (k-means k = ", k, ")"),
heatmap_legend_param = list(title = "z")
)

pdf(paste0(path, "Figures/ExtFig3c_Heatmap_k", k, "log2HL_log2T_log2RPF_log2RNAseq_log2riBAQH4h2.pdf", sep=""))
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()


# ---- Figure 3e ----
# GO per multiomics cluster
all_genes <- names(clusters)
all_genes <- gsub("_.*$", "", all_genes)  # remove suffixes after underscore

bg_map <- bitr(
  all_genes,
  fromType = "UNIPROT",
  toType   = "ENTREZID",
  OrgDb    = org.Mm.eg.db
)

background_entrez <- unique(bg_map$ENTREZID)

# split genes by cluster
cluster_genes <- split(names(clusters), clusters)
cluster_genes <- lapply(
  cluster_genes,
  \(x) gsub("_.*$", "", x)
)
cluster_genes_df <- data.frame(
  Protein = names(clusters),
  Cluster = clusters
) 

# map each cluster to ENTREZ IDs
cluster_entrez <- lapply(cluster_genes, function(g) {
  m <- bitr(
    g,
    fromType = "UNIPROT",
    toType   = "ENTREZID",
    OrgDb    = org.Mm.eg.db
  )
  unique(m$ENTREZID)
})

# Remove very small clusters
cluster_entrez <- cluster_entrez[sapply(cluster_entrez, length) >= 20]

# run GO enrichment across clusters
# for (i in c("BP", "MF", "CC")) {
i = "CC"  # or "MF", "BP"
  cc <- compareCluster(
    geneCluster   = cluster_entrez,
    fun           = "enrichGO",
    universe      = background_entrez,
    OrgDb         = org.Mm.eg.db,
    keyType       = "ENTREZID",
    ont           = i,   # or "MF", "CC", "BP"
    pAdjustMethod = "BH",
    pvalueCutoff  = 1, #0.05,
    qvalueCutoff  = 1, #0.2,
    readable      = TRUE
  )

  cc@compareClusterResult <- cc@compareClusterResult %>%
  dplyr::filter(
    p.adjust < 0.05,
    Count >= 3
  ) %>%
  dplyr::group_by(Cluster, Description) %>%
  dplyr::slice_min(p.adjust, n = 100, with_ties = FALSE) %>%
  dplyr::ungroup()

  cc_simplified <- clusterProfiler::simplify(
    cc, 
    cutoff = 0.7,    # Similarity threshold (lower = more grouping/fewer terms)
    by = "p.adjust", # Keep the term with the best p-value in the group
    select_fun = min
  )

  # results table
  cc_df <- as.data.frame(cc)
  write.csv(cc_df, paste0(path, "Objects/GO_k",k,"_log2HL_log2T_log2RPF_log2RNAseq_log2riBAQ4h_",i,"_noPvalueFilter.csv", sep=""), row.names = FALSE)

  # dotplot with all clusters
  
  p <- dotplot(cc, showCategory = 5) +
  scale_color_continuous(
      low = "red", high = "blue", 
      limits = c(0, 0.05),  
      name = "p.adjust"
    ) +
    # 2. Fix the Size Scale (GeneRatio)
    scale_size_continuous(
      limits = c(0, 0.2), 
      range = c(2, 8),        # Adjust range to control actual dot diameters
      name = "GeneRatio"
    ) +
  theme(axis.text.x = element_text(size = 12, hjust = 1), # Resize gene names
          axis.text.y = element_text(size = 8))  

  pdf(paste0(path, "Figures/Fig3e_GO_k",k,"_log2HL_log2T_log2RPF_log2RNAseq_log2riBAQ4h_",i,"2.pdf", sep=""))
  print(p)
  dev.off()
#}


# ---- Figure 3d and Extended Figure 3b ----
# Multiomics line plot per cluster
cluster_genes_df <- data.frame(
  Protein = names(clusters),
  Cluster = unname(clusters),
  stringsAsFactors = FALSE
) %>%
  filter(
    !is.na(Protein),
    Protein != "",
    Protein != "NA"
  )

mat2_df <- mat %>% 
  as.data.frame() %>%
  rownames_to_column("Protein") %>%
  left_join(cluster_genes_df, by = "Protein") %>%
  pivot_longer(cols = -c(Protein, Cluster), names_to = "Feat", values_to = "Z") %>%
  separate(Feat, into = c("Experiment", "Embryonic_Age"), sep = "_") %>%
  mutate(Cluster = factor(Cluster, levels = sort(unique(Cluster))),
          Embryonic_Age_numeric = as.numeric(str_remove(Embryonic_Age, "^E")))

df_summary <- mat2_df %>%
  group_by(Cluster, Experiment, Embryonic_Age_numeric) %>%
  summarise(
    mean_val = mean(Z, na.rm = TRUE),
    sd_val = sd(Z, na.rm = TRUE),
    .groups = "drop"
  )

p <- ggplot(df_summary, aes(x = Embryonic_Age_numeric, y = mean_val, color = Experiment, fill = Experiment, group = Experiment)) +
  # Draw the ribbon (Mean +/- SD)
  geom_ribbon(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val), 
              alpha = 0.2, color = NA) +
  # Draw the mean line
  geom_line() +
  # Add points for clarity
  geom_point(size = 2) +
  # Aesthetics
  theme_minimal() +
  labs(
    title = "Trends of HL, LFQ T, RPF, riBAQ H 4h, RNAseq (Mean ± SD)",
    x = "Embryonic Age",
    y = "Z-score Intensity",
    color = "Assay Type",
    fill = "Assay Type"
  ) +
  scale_color_manual(values = c("HL" = "#00AEEF", "T" = "#231F20", "RPF" = "#F7941D", "H" = "#ED1C24", "RNAseq" = "#92D050")) +
  scale_fill_manual(values = c("HL" = "#00AEEF", "T" = "#231F20", "RPF" = "#F7941D", "H" = "#ED1C24", "RNAseq" = "#92D050")) +
  theme(legend.position = "bottom") +
  facet_grid(Cluster ~ ., switch = "both")  # separate plot per cluster

pdf(paste0(path, "Figures/ExtFig3d_LinePlots_k",k,"_log2HL_log2LFQT_log2RPF_log2RNAseq_log2riBAQH4h.pdf", sep=""), width = 6, height = 15)
print(p)
dev.off()


# ---- Figure 3e ----
# Multiomics line plot per cluster for specific proteins of interest
proteins <- c("Q9D1C9", "Q640M1", "Q0V8M0", "Q8BW10", # "Rrp7a", "Utp14a", "Kri1", "Nob1",
  "Q64336", "Q60632", "P35922", "Q9D3A8",             # "Tbr1", "Nr2f1", "Fmr1", "Nos1ap",
  "Q8K209", "P63015", "P53783", "P13864",             # "Adgrg1", "Pax6", "Sox1", "Dnmt1",
  "O89086", "Q9JLN9", "Q8K4Q0", "Q8BSK8")             # "Rbm3", "Mtor", "Rptor", "Rps6kb1"
     

cluster_genes_df <- data.frame(
  Protein = names(clusters),
  Cluster = unname(clusters),
  stringsAsFactors = FALSE
) %>%
  filter(
    !is.na(Protein),
    Protein != "",
    Protein != "NA"
  ) 

mat2_df <- mat %>% 
  as.data.frame() %>%
  rownames_to_column("Protein") %>%
  separate(Protein, into = c("Protein", "Gene"), sep = "_") %>%
  left_join(cluster_genes_df, by = "Protein") %>%
  pivot_longer(cols = -c(Protein, Cluster, Gene), names_to = "Feat", values_to = "Z") %>%
  separate(Feat, into = c("Experiment", "Embryonic_Age"), sep = "_") %>%
  filter(Protein %in% proteins)


mat_sd_df <- mat_sd %>% 
  as.data.frame() %>%
  rownames_to_column("Protein") %>%
  separate(Protein, into = c("Protein", "Gene"), sep = "_") %>%
  left_join(cluster_genes_df, by = "Protein") %>%
  pivot_longer(cols = -c(Protein, Cluster, Gene), names_to = "Feat", values_to = "SD") %>%
  separate(Feat, into = c("Experiment", "Embryonic_Age"), sep = "_") %>%
  filter(Protein %in% proteins)

plot_data <- merge(mat2_df, mat_sd_df, by = c("Protein", "Cluster", "Gene", "Experiment", "Embryonic_Age")) 
plot_data <- plot_data%>%
                mutate(Cluster = factor(Cluster, levels = sort(unique(Cluster))),
                        Embryonic_Age_numeric = as.numeric(str_remove(Embryonic_Age, "^E")))


p <- ggplot(plot_data, aes(x = Embryonic_Age_numeric, y = Z, color = Experiment, fill = Experiment, group = Experiment)) +
  # Draw the mean line
  geom_line() +
  # Draw the ribbon (Mean +/- SD)
  geom_ribbon(aes(ymin = Z - SD, ymax = Z + SD), 
              alpha = 0.2, color = NA) +
  # Add points for clarity
  geom_point(size = 2) +
  # Aesthetics
  theme_minimal() +
  labs(
    title = "Trends of HL, LFQ T, RPF, riBAQ H 4h, RNAseq (Mean ± SD) for selected proteins",
    x = "Embryonic Age",
    y = "Z-score Intensity",
    color = "Assay Type"
  ) +
  scale_color_manual(values = c("HL" = "#00AEEF", "T" = "#231F20", "RPF" = "#F7941D", "H" = "#ED1C24", "RNAseq" = "#92D050")) +
  scale_fill_manual(values = c("HL" = "#00AEEF", "T" = "#231F20", "RPF" = "#F7941D", "H" = "#ED1C24", "RNAseq" = "#92D050")) +
  theme(legend.position = "bottom") +
  facet_grid(Protein ~ ., switch = "both")  # separate plot per cluster

pdf(paste0(path, "Figures/ExtFig3e_LinePlots_k",k,"_log2HL_log2LFQT_log2RPF_log2RNAseq_log2riBAQH4h_SelectedProteins.pdf", sep=""), width = 6, height = 20)
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
