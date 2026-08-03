#########################################################################################################################################################################
# Riboseq and SILAC proteomics comparison
#
# Purpose:
#   Linear correlation, clustering and GO enrichment for Riboseq and riBAQ heavy 4h.
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
  paste0("03b_Fig2_RiboseqProteomics_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)

# Create the log file and write session information
cat(
  "Riboseq-Proteomics analysis log\n",
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

k = 10 # Riboseq + proteomics clusters after elbow plot

mat <- readRDS(file.path(path, "Objects/all_data.rds"))

mat_ribo <- mat[mat$Experiment %in% c("RPF", "riBAQ_H") & mat$Incubation_Time == "4h", ]
mat_ribo <- unique(mat_ribo)


# ---- Figure 2c ----
# Compare Riboseq and heavy proteomics layers by stage-specific correlations
plot_data <- mat_ribo %>%
  filter(Embryonic_Age %in% c("E12.5", "E15.5")) %>% 
  mutate(Value = na_if(Value, 0)) %>% # Change 0 to NA for log transformation
  pivot_wider(id_cols = c(protein_group, genes, Ensembl_id, Embryonic_Age, Replicates), names_from = Experiment, values_from = Value) %>%
  drop_na(RPF, riBAQ_H) %>%
  mutate(log2_RPF = log2(RPF), log2_H = log2(riBAQ_H))

stats_df <- plot_data %>%
  group_by(Embryonic_Age) %>%
  summarise(
    model = list(lm(log2_RPF ~ log2_H)),
    intercept = coef(model[[1]])[1],
    slope = coef(model[[1]])[2],
    resid_sd = sd(residuals(model[[1]])),
    .groups = "drop"
  )

p <- ggplot(plot_data, aes(x = log2_H, y = log2_RPF)) +
  geom_pointdensity(adjust = 0.1) + 
  scale_color_viridis_c() +
  geom_smooth(method = "lm", color = "black", se = FALSE, size = 1) + # Main Regression Line
  geom_abline(data = stats_df,                                        # Add the +2 SD line
              aes(intercept = intercept + (2 * resid_sd), slope = slope), 
              color = "red", linetype = "dashed") +
  geom_abline(data = stats_df,                                        # Add the -2 SD line
              aes(intercept = intercept - (2 * resid_sd), slope = slope), 
              color = "red", linetype = "dashed") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  facet_wrap(~Embryonic_Age) +
  theme_minimal() +
  labs(
    title = "RPF vs H with Outlier Boundaries",
    subtitle = "Red dashed lines represent ±2 Standard Deviations of the residuals",
    x = "log2(H)",
    y = "log2(RPF)",
    color = "Density"
  )

pdf(paste(path, "Figures/Fig2c_Correlation_log2RPF_log2riBAQH_E12.5_E15.5_SD2.pdf", sep = ""))
print(p)
dev.off()


# ---- Fig 2d ----
# Correlation of fold changes E15.5 vs E12.5 between Riboseq and heavy proteomics
plot_data_FC <- mat_ribo %>%
  filter(Embryonic_Age %in% c("E12.5", "E15.5")) %>%
  mutate(Value = na_if(Value, 0)) %>% # Change 0 to NA for log transformation
  pivot_wider(names_from = Experiment, values_from = Value) %>%
  mutate(log2_RPF = log2(RPF), log2_H = log2(riBAQ_H)) %>%
  drop_na(log2_RPF, log2_H) %>%
  group_by(protein_group) %>%
  summarise(Log2FC_E15.5_vs_E12.5_RPF = log2(RPF[Embryonic_Age == "E15.5"] / RPF[Embryonic_Age == "E12.5"]),
        Log2FC_E15.5_vs_E12.5_H = log2(riBAQ_H[Embryonic_Age == "E15.5"] / riBAQ_H[Embryonic_Age == "E12.5"]), .groups = "drop")

p <- ggplot(plot_data_FC, aes(x = Log2FC_E15.5_vs_E12.5_H, y = Log2FC_E15.5_vs_E12.5_RPF)) +
  geom_pointdensity(adjust = 0.1, method = "neighbors") + 
  scale_color_viridis_c() +
  geom_smooth(method = "lm", color = "black", se = TRUE, linetype = "dashed") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_minimal() +
  labs(
    title = "Correlation between log2FC E15.5 vs E12.5 RPF and H",
    subtitle = "Spearman correlation for log2FC E15.5 and E12.5",
    x = "Log2 Fold Change (H)",
    y = "Log2 Fold Change (RPF)",
    color = "Density"
  )

pdf(paste(path, "Figures/Fig2d_Correlation_Log2FCE15vsE12_RPF_riBAQH.pdf", sep = ""))
print(p)
dev.off()


# ---- Figure 2e ----
# Run GO enrichment for genes/proteins identified as outliers in the correlation plots.
outliers_extracted <- plot_data %>%
  group_by(Embryonic_Age) %>%
  mutate(
    model_fit = predict(lm(log2_RPF ~ log2_H)),     # Fit the linear model
    residual = log2_RPF - model_fit,                # Calculate the residual (Actual - Predicted)
    sd_threshold = sd(residual)                     # Calculate the SD of residuals
  ) %>%
  filter(abs(residual) > (2 * sd_threshold)) %>%
  mutate(Status = if_else(residual > 0, "Up", "Down")) %>%
  arrange(Embryonic_Age, desc(abs(residual)))

# GO of the outliers 
all_genes <- unlist(unique(plot_data$protein_group))
bg_map <- bitr(
  all_genes,
  fromType = "UNIPROT",
  toType   = "ENTREZID",
  OrgDb    = org.Mm.eg.db
)
background_entrez <- unique(bg_map$ENTREZID)

# split genes by cluster
outliers_extracted$GO_cluster <- paste(outliers_extracted$Embryonic_Age, outliers_extracted$Status, sep = "_")
cluster_genes <- split(outliers_extracted[, "protein_group"], outliers_extracted$GO_cluster)

# map each cluster to ENTREZ IDs
cluster_entrez <- lapply(cluster_genes, function(g) {
  ids_to_map <- as.character(g$protein_group)
  
  m <- bitr(
    ids_to_map,
    fromType = "UNIPROT",
    toType   = "ENTREZID",
    OrgDb    = org.Mm.eg.db
  )
  return(unique(m$ENTREZID))
})

# Remove very small clusters
cluster_entrez <- cluster_entrez[sapply(cluster_entrez, length) >= 5]

# run GO enrichment across clusters
cc <- compareCluster(
  geneCluster   = cluster_entrez,
  fun           = "enrichGO",
  universe      = background_entrez,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF", #BP, CC, MF
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

# results table
cc_df <- as.data.frame(cc)
write.csv(cc_df, paste(path, "Objects/Fig2e_GO_outliers_log2RPF_log2riBAQH4h_MF.csv", sep = ""), row.names = FALSE) #BP, CC, MF

# dotplot with all clusters
p <- dotplot(cc, showCategory = 5, ) +
 theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1), # Resize gene names
        axis.text.y = element_text(size = 8))    

pdf(paste(path, "Figures/Fig2e_GO_outliers_log2RPF_log2riBAQH4h_MF.pdf", sep = "")) #BP, CC, MF
print(p)
dev.off()


# ---- Figure 2f ----
# Z-score proteomics and Ribo-seq tables 
per_protein_age <- mat_ribo %>%
  group_by(protein_group, Experiment,Embryonic_Age) %>%
  filter(sum(!is.na(Value)) >= 2) %>%
  summarise(
    Value = median(Value, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(Age_numeric = as.numeric(sub("E|P", "", Embryonic_Age))) %>%
  filter(Embryonic_Age != "P0")  # Remove E19.5 for clustering, as it is only available for Ribo-seq

# Z-score each protein/gene across ages so temporal shapes are comparable across assays.
zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))  # constant/1-point series -> 0
  (x - mean(x, na.rm = TRUE)) / s
}

dat_z <- per_protein_age %>%
  group_by(protein_group, Experiment) %>%
  arrange(Embryonic_Age, .by_group = TRUE) %>%
  mutate(Zscore_log2_values = zscore(log2(Value))) %>%
  ungroup()

mat_wide <- dat_z %>%
  unite(col = "Feat", Experiment, Embryonic_Age, sep = "_") %>%
  dplyr::select(protein_group, Feat, Zscore_log2_values) %>%
  pivot_wider(names_from = Feat, values_from = Zscore_log2_values)

mat_wide <- mat_wide %>%
    filter(!is.na(protein_group)) %>%
    column_to_rownames("protein_group") %>%
    as.matrix()

mat_wide <- mat_wide[complete.cases(mat_wide), , drop = FALSE]


# Cluster proteins by temporal profile: Heatmap
col_df <- data.frame(
  assay = sub("_.*$", "", colnames(mat_wide)),                
  Embryonic_Age   = sub("^.*_E", "E", colnames(mat_wide)),               
  stringsAsFactors = FALSE
)

col_df$age <- factor(col_df$Embryonic_Age)
col_df$assay <- factor(col_df$assay)

mat_wide <- mat_wide[, order(col_df$assay, col_df$age), drop = FALSE]
col_df <- col_df[order(col_df$assay, col_df$age), , drop = FALSE]

# Clustering: Use this for reproducibility
set.seed(123)

kmeans_path <- file.path(paste0(path,"/Objects/Fig2f_fixed_clusters_k", k, ".rds"))

if (file.exists(kmeans_path)) {
  km <- readRDS(kmeans_path)
} else {
  km <- kmeans(mat_wide, centers = k, nstart = 50, iter.max = 100)

  saveRDS(km, kmeans_path)
}

clusters <- km$cluster
centroids <- km$centers

if (!setequal(rownames(mat_wide), names(clusters))) {
  stop(
    "The proteins in mat_wide do not match the proteins used for the ",
    "published clustering."
  )
}


# order within each cluster by distance to centroid (cleaner blocks)
dist_to_centroid <- vapply(seq_len(nrow(mat_wide)), function(i) {
cl <- clusters[i]
sum((mat_wide[i, ] - centroids[cl, ])^2)
}, numeric(1))

ord <- order(clusters, dist_to_centroid)
mat_z_ord <- mat_wide[ord, , drop = FALSE]

row_ha <- rowAnnotation(
KMeans = factor(clusters[ord], levels = sort(unique(clusters))),
annotation_name_side = "top"
)

# Split heatmap by kmeans cluster (nice for readability)
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

pdf(paste(path, "Figures/Fig2f_Heatmap_k", k, "_log2RPF_log2riBAQH4h2.pdf", sep = ""))
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()


# ---- Figure 2g ----
# GO enrichment for temporal clusters
all_genes <- names(clusters)

bg_map <- bitr(
  all_genes,
  fromType = "UNIPROT",
  toType   = "ENTREZID",
  OrgDb    = org.Mm.eg.db
)

background_entrez <- unique(bg_map$ENTREZID)

# split genes by cluster
cluster_genes <- split(names(clusters), clusters)
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
cluster_entrez <- cluster_entrez[sapply(cluster_entrez, length) >= 5]

# run GO enrichment across clusters
# for (i in c("BP", "MF", "CC")) {
i = "BP"  # or "MF", "CC"
  cc <- compareCluster(
    geneCluster   = cluster_entrez,
    fun           = "enrichGO",
    universe      = background_entrez,
    OrgDb         = org.Mm.eg.db,
    keyType       = "ENTREZID",
    ont           = i,   # or "MF", "CC", "BP"
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )

  cc@compareClusterResult <- cc@compareClusterResult %>%
  dplyr::filter(
    p.adjust < 0.05,
    Count >= 3
  ) %>%
  dplyr::group_by(Cluster, ONTOLOGY) %>%
  dplyr::slice_min(p.adjust, n = 100, with_ties = FALSE) %>%
  dplyr::ungroup()

  cc_simplified <- clusterProfiler::simplify(
    cc,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )

  # results table
  cc_df <- as.data.frame(cc)
  
  # write.csv(cc_df, paste(path, "Objects/Fig2g_GO_k",k,"_log2RPF_log2riBAQH4h_",i,".csv", sep = ""), row.names = FALSE)

  # dotplot with all clusters
  
  p <- dotplot(cc_simplified, showCategory = 5, ) +
  theme(axis.text.x = element_text(size = 12, hjust = 1), # Resize gene names
          axis.text.y = element_text(size = 8))  

  pdf(paste(path, "Figures/Fig2g_GO_k",k,"_log2RPF_log2riBAQH4h_",i,".pdf", sep = ""))
  print(p)
  dev.off()
#}


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
