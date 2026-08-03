#########################################################################################################################################################################
# Single-cell SILAC analysis
#
# Purpose:
#   Preprocess four scSILAC-derived data types, cluster cells using NA-aware proDA distances,
#   identify cluster markers and generate Figure 5 plots.
#
# GitHub/reuse notes:
#   - Input values are assumed to be log10-transformed and cell-size-normalized.
#   - Missing values are retained and are not imputed.
#   - Update the paths and analysis settings before running on another machine.
#########################################################################################################################################################################


# ---- Logging and warning control ----
log_dir <- "./Output/Logs/"
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  log_dir,
  paste0("05_Fig5_scSILAC_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)

# Create the log file and write session information
cat(
  "Single cell SILAC analysis log\n",
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
  library(arrow)
  library(dplyr)
  library(enrichR)
  library(forcats)
  library(ggplot2)
  library(ggpubr)
  library(ggpointdensity)
  library(Matrix)
  library(matrixStats)
  library(patchwork)
  library(proDA)
  library(Seurat)
  library(tibble)
  library(tidyr)
})


# ---- User settings ----
path <- "./Output/"
dir.create(path, recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(path, "Objects/"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(path, "Figures/"), recursive = TRUE, showWarnings = FALSE)

data_path <- "./Data/SILAC_SingleCellProteomics_MS/"

# Filtering
min_cells_per_protein <- 5
min_proteins_per_cell <- 100

# Dimensional reduction and clustering
n_dimensions <- 20
resolution <- 0.92
target_n_clusters <- 3

# Differential-abundance thresholds
n_top_markers <- 20
marker_log2fc_threshold <- 0.6
marker_adj_p_threshold <- 0.05

# Proteins displayed on UMAPs
selected_markers <- c("Ncam1", "Plch1", "Eif3g", "Rpl18a")

set.seed(123)


# ---- Assay definitions ----
assay_config <- list(
  total = list(column = "log10_total_intensity", label = "Total proteome"),
  light = list(column = "log10_LH_ratio", label = "Light proteome (L/H)"),
  medium = list(column = "log10_MH_ratio", label = "Medium-labelled proteome (M/H)"),
  turnover = list(column = "log10_ML_ratio", label = "Turnover proxy (M/L)")
)


# ---- Minimal helper functions ----
# These operations are repeated for all four assays and are kept as functions
# to avoid duplicating long blocks of analysis code.

# Construct a protein-by-cell matrix while retaining missing values.
build_na_matrix <- function(data, value_column, all_cells) {
  long_data <- data %>%
    transmute(
      protein = Genes,
      cell = Well,
      value = .data[[value_column]]
    ) %>%
    filter(!is.na(protein), protein != "", !is.na(cell)) %>%
    group_by(protein, cell) %>%
    summarise(
      value = if (all(is.na(value))) NA_real_ else median(value, na.rm = TRUE),
      .groups = "drop"
    )

  proteins <- sort(unique(long_data$protein))
  values <- matrix(
    NA_real_,
    nrow = length(proteins),
    ncol = length(all_cells),
    dimnames = list(proteins, all_cells)
  )

  observed <- !is.na(long_data$value)
  values[cbind(
    match(long_data$protein[observed], proteins),
    match(long_data$cell[observed], all_cells)
  )] <- long_data$value[observed]

  keep_proteins <- rowSums(!is.na(values)) >= min_cells_per_protein
  values[keep_proteins, , drop = FALSE]
}


# Perform one-versus-rest differential abundance with proDA.
run_proda_markers <- function(values, clusters) {
  clusters <- factor(clusters)
  names(clusters) <- colnames(values)

  bind_rows(lapply(levels(clusters), function(cluster_id) {
    group <- factor(ifelse(clusters == cluster_id, "in", "out"), levels = c("out", "in"))
    sample_data <- data.frame(group = group, row.names = colnames(values))

    fit <- proDA(
      values,
      design = ~ group,
      col_data = sample_data,
      data_is_log_transformed = TRUE,
      verbose = FALSE
    )

    test_results <- test_diff(
      fit,
      contrast = "groupin",
      pval_adjust_method = "BH",
      sort_by = NULL
    ) %>%
      transmute(
        gene = name,
        cluster = as.character(cluster_id),
        p_val = pval,
        p_val_adj = adj_pval,
        mean_diff_log10 = diff,
        log2FC = diff * log2(10),
        proDA_t = t_statistic,
        proDA_se = se,
        proDA_df = df,
        proDA_avg_abundance = avg_abundance,
        proDA_n_approx = n_approx,
        proDA_n_obs = n_obs
      )

    cells_in <- names(clusters)[clusters == cluster_id]
    cells_out <- names(clusters)[clusters != cluster_id]
    values_in <- values[, cells_in, drop = FALSE]
    values_out <- values[, cells_out, drop = FALSE]

    descriptive_statistics <- tibble(
      gene = rownames(values),
      pct_in = rowMeans(!is.na(values_in)),
      pct_out = rowMeans(!is.na(values_out)),
      mean_log10_detected_in = rowMeans(values_in, na.rm = TRUE),
      mean_log10_detected_out = rowMeans(values_out, na.rm = TRUE)
    ) %>%
      mutate(
        mean_log10_detected_in = ifelse(is.nan(mean_log10_detected_in), NA_real_, mean_log10_detected_in),
        mean_log10_detected_out = ifelse(is.nan(mean_log10_detected_out), NA_real_, mean_log10_detected_out),
        detected_log2FC = (mean_log10_detected_in - mean_log10_detected_out) * log2(10)
      )

    left_join(test_results, descriptive_statistics, by = "gene")
  })) %>%
    arrange(cluster, p_val_adj, desc(abs(log2FC)))
}


# Plot selected protein values on an existing UMAP.
plot_selected_proteins <- function(object, values, cluster_column, umap_name, assay_label, output_file) {
  genes <- intersect(selected_markers, rownames(values))
  cells <- intersect(colnames(object), colnames(values))

  if (length(genes) == 0 || length(cells) == 0) return(invisible(NULL))

  object_subset <- subset(object, cells = cells)
  coordinates <- Embeddings(object_subset, umap_name)[cells, 1:2, drop = FALSE]

  cluster_plot <- DimPlot(
    object_subset,
    reduction = umap_name,
    group.by = cluster_column,
    label = TRUE,
    repel = TRUE,
    pt.size = 1
  ) +
    ggtitle(paste0(assay_label, ": clusters")) +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  marker_plots <- lapply(genes, function(gene) {
    plot_data <- data.frame(
      UMAP_1 = coordinates[, 1],
      UMAP_2 = coordinates[, 2],
      value = as.numeric(values[gene, cells])
    )

    ggplot(plot_data, aes(UMAP_1, UMAP_2)) +
      geom_point(color = "grey80", size = 0.9) +
      geom_point(
        data = filter(plot_data, !is.na(value)),
        aes(color = value),
        size = 0.9
      ) +
      scale_color_viridis_c(option = "inferno", na.value = "grey80") +
      coord_equal() +
      theme_classic() +
      labs(title = gene, color = "log10 value", x = "UMAP 1", y = "UMAP 2") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
  })

  combined_plot <- cluster_plot | wrap_plots(marker_plots, ncol = 2)

  ggsave(
    output_file,
    combined_plot,
    width = 16,
    height = max(7, 3.2 * ceiling(length(genes) / 2)),
    limitsize = FALSE
  )
}





# ---- Load and prepare the data ----
raw_data <- read_parquet(file.path(data_path, "clean+TechNoiseSafety_unsafeRemoved_CSNorm.parquet"))

required_columns <- c(
  "Well",
  "Genes",
  "CS.Residual.M.Intensity",
  "CS.Residual.L.Intensity",
  "CS.Residual.Stacked.LandM.Intensity"
)

missing_columns <- setdiff(required_columns, colnames(raw_data))
if (length(missing_columns) > 0) {
  stop("The following required columns are missing: ", paste(missing_columns, collapse = ", "))
}

data <- raw_data %>%
  transmute(
    Well,
    Genes,
    log10_MH_ratio = CS.Residual.M.Intensity,
    log10_LH_ratio = CS.Residual.L.Intensity,
    log10_total_intensity = CS.Residual.Stacked.LandM.Intensity,
    log10_ML_ratio = CS.Residual.M.Intensity - CS.Residual.L.Intensity
  )

all_cells <- sort(unique(data$Well))


# ---- Construct and filter assay matrices ----
assay_matrices <- lapply(assay_config, function(config) {
  values <- build_na_matrix(data, config$column, all_cells)
  list(values = values, presence = !is.na(values))
})

# Retain cells passing the protein-count threshold in every assay.
cells_pass_by_assay <- lapply(assay_matrices, function(x) {
  colnames(x$values)[colSums(x$presence) >= min_proteins_per_cell]
})

common_cells <- sort(Reduce(intersect, cells_pass_by_assay))
if (length(common_cells) < 3) stop("Fewer than three cells pass filtering in every assay.")

assay_matrices <- lapply(assay_matrices, function(x) {
  x$values <- x$values[, common_cells, drop = FALSE]
  x$presence <- x$presence[, common_cells, drop = FALSE]
  x
})

message("Cells retained in every assay: ", length(common_cells))
saveRDS(assay_matrices, file.path(path, "Objects/assay_matrices_all_assays.rds"))



# ---- Figure 5a ----
# Number of quantified proteins per retained cell and data type.
counts_per_cell <- bind_rows(lapply(names(assay_matrices), function(assay_name) {
  tibble(
    cell = colnames(assay_matrices[[assay_name]]$presence),
    assay = assay_name,
    n_proteins = colSums(assay_matrices[[assay_name]]$presence)
  )
})) %>%
  mutate(
    assay = factor(
      assay,
      levels = c("light", "medium", "turnover", "total"),
      labels = c("Light", "Medium", "Turnover", "Total")
    )
  )

summary_counts <- counts_per_cell %>%
  group_by(assay) %>%
  summarise(
    mean_n_proteins = mean(n_proteins),
    sd_n_proteins = sd(n_proteins),
    sem_n_proteins = sd(n_proteins) / sqrt(n()),
    n_cells = n(),
    .groups = "drop"
  )

p <- ggplot(summary_counts, aes(assay, mean_n_proteins)) +
  geom_col(width = 0.7, fill = "grey80", color = "black") +
  geom_errorbar(
    aes(
      ymin = pmax(0, mean_n_proteins - sd_n_proteins),
      ymax = mean_n_proteins + sd_n_proteins
    ),
    width = 0.15
  ) +
  geom_jitter(
    data = counts_per_cell,
    aes(assay, n_proteins),
    width = 0.12,
    size = 1.8,
    alpha = 0.7,
    inherit.aes = FALSE
  ) +
  theme_classic() +
  labs(
    title = "Protein detection after cell filtering",
    x = NULL,
    y = "Detected proteins per cell"
  )

ggsave(file.path(path, "Figures/Fig5a_ProteinCounts.pdf"), p, width = 6, height = 5)


# ---- Create the common Seurat container ----
# Seurat stores cells, reductions, graphs and cluster assignments. Each assay
# starts from an independent copy of the same cell container.
base_counts <- Matrix(
  assay_matrices$total$presence * 1,
  sparse = TRUE,
  dimnames = dimnames(assay_matrices$total$presence)
)

base_object <- CreateSeuratObject(
  counts = base_counts,
  assay = "detection",
  project = "scSILAC_E145_neocortex",
  min.cells = 0,
  min.features = 0
)
base_object$Well <- colnames(base_object)



# ---- Figures 5c ----
# Clustering on total values
assay_name = "total"
assay_label <- assay_config[[assay_name]]$label

message("Processing ", assay_name, " (", assay_label, ")")

values <- assay_matrices[[assay_name]]$values
assay_object <- base_object

# Calculate NA-aware sample distances with proDA.
distance_estimate <- dist_approx(values)
distance_mean <- as.matrix(distance_estimate$mean)
distance_mean <- (distance_mean + t(distance_mean)) / 2
diag(distance_mean) <- 0

# Classical multidimensional scaling.
n_mds_dimensions <- min(n_dimensions, ncol(values) - 1)
mds <- cmdscale(as.dist(distance_mean), k = n_mds_dimensions, eig = TRUE)
mds_embedding <- as.matrix(mds$points)
colnames(mds_embedding) <- paste0("MDS_", seq_len(ncol(mds_embedding)))

reduction_name <- paste0("mds_", assay_name)
umap_name <- paste0("umap_", assay_name)
nn_name <- paste0(assay_name, "_nn")
snn_name <- paste0(assay_name, "_snn")
cluster_column <- paste0("clusters_", assay_name)

assay_object[[reduction_name]] <- CreateDimReducObject(
  embeddings = mds_embedding,
  key = paste0("MDS", assay_name, "_"),
  assay = DefaultAssay(assay_object)
)

dimensions_used <- seq_len(ncol(mds_embedding))
assay_object <- FindNeighbors(
  assay_object,
  reduction = reduction_name,
  dims = dimensions_used,
  graph.name = c(nn_name, snn_name),
  verbose = FALSE
)

assay_object <- FindClusters(
  assay_object,
  graph.name = snn_name,
  resolution = resolution,
  cluster.name = cluster_column,
  random.seed = 1,
  verbose = FALSE
)

assay_object <- RunUMAP(
  assay_object,
  reduction = reduction_name,
  dims = dimensions_used,
  reduction.name = umap_name,
  reduction.key = paste0("UMAP", assay_name, "_"),
  seed.use = 1,
  verbose = FALSE
)


# Display median total, M/L and M/H values per cell on the total-proteome UMAP.
total_object <- assay_object

total_values <- assay_matrices$total$values
ML_values <- assay_matrices$turnover$values
MH_values <- assay_matrices$medium$values

cells <- Reduce(
  intersect,
  list(colnames(total_object), colnames(total_values), colnames(ML_values), colnames(MH_values))
)
if (length(cells) == 0) stop("No common cells were found between the Seurat object and the matrices.")

total_object_subset <- subset(total_object, cells = cells)
umap_coordinates <- Embeddings(total_object_subset, "umap_total")
cells <- rownames(umap_coordinates)
umap_coordinates <- umap_coordinates[cells, 1:2, drop = FALSE]

median_total <- colMedians(as.matrix(total_values[, cells, drop = FALSE]), na.rm = TRUE)
median_ML <- colMedians(as.matrix(ML_values[, cells, drop = FALSE]), na.rm = TRUE)
median_MH <- colMedians(as.matrix(MH_values[, cells, drop = FALSE]), na.rm = TRUE)

median_total[is.nan(median_total)] <- NA_real_
median_ML[is.nan(median_ML)] <- NA_real_
median_MH[is.nan(median_MH)] <- NA_real_

plot_data <- data.frame(
  cell = cells,
  UMAP_1 = umap_coordinates[, 1],
  UMAP_2 = umap_coordinates[, 2],
  Median_total = median_total[cells],
  Median_ML = median_ML[cells],
  Median_MH = median_MH[cells]
)

cluster_plot <- DimPlot(
  total_object_subset,
  reduction = "umap_total",
  group.by = "clusters_total",
  label = TRUE,
  repel = TRUE,
  pt.size = 1
) +
  ggtitle("Total-proteome clusters") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# Convert the three median columns to long format so one ggplot call creates all panels.
median_plot_data <- plot_data %>%
  pivot_longer(
    cols = c(Median_total, Median_ML, Median_MH),
    names_to = "measurement",
    values_to = "value"
  ) %>%
  mutate(
    measurement = factor(
      measurement,
      levels = c("Median_total", "Median_ML", "Median_MH"),
      labels = c("Median total", "Median M/L", "Median M/H")
    )
  )

median_plots <- lapply(levels(median_plot_data$measurement), function(measurement_name) {
  current_data <- filter(median_plot_data, measurement == measurement_name)

  ggplot(current_data, aes(UMAP_1, UMAP_2)) +
    geom_point(color = "grey80", size = 0.9) +
    geom_point(
      data = filter(current_data, !is.na(value)),
      aes(color = value),
      size = 0.9
    ) +
    scale_color_viridis_c(option = "inferno", na.value = "grey80") +
    coord_equal() +
    theme_classic() +
    labs(title = measurement_name, color = "Median value", x = "UMAP 1", y = "UMAP 2") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
})

combined_plot <- (cluster_plot | median_plots[[1]] | median_plots[[2]] | median_plots[[3]]) +
  plot_layout(ncol = 2, nrow = 2) +
  plot_annotation(
    title = "Median total, M/L and M/H on total-proteome UMAP",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5))
  )

ggsave(
  file.path(path, "Figures/Fig5c_Median_total_ML_MH_on_total_UMAP.pdf"),
  combined_plot,
  width = 14,
  height = 12,
  limitsize = FALSE
)

saveRDS(assay_object, file.path(path, paste0("Objects/assay_object_",assay_name,".rds")))



# ---- Figures 5d ----
# Plotting selected markers M/L, M/H and total intensities on total-proteome UMAP
# Differential abundance: each cluster versus all remaining cells.
# clusters <- assay_object[[cluster_column, drop = TRUE]]
# names(clusters) <- colnames(assay_object)
# markers <- run_proda_markers(values, clusters)

# # Selected protein values on the assay-specific UMAP.
total_object <- assay_object

plot_selected_proteins(
  total_object,
  assay_matrices$turnover$values,
  "clusters_total",
  "umap_total",
  "M/L on total-proteome clusters",
  file.path(path, "Figures/Fig5d_SelectedMarkers_ML_on_total_clusters.pdf")
)

plot_selected_proteins(
  total_object,
  assay_matrices$medium$values,
  "clusters_total",
  "umap_total",
  "M/H on total-proteome clusters",
  file.path(path, "Figures/Fig5d_SelectedMarkers_MH_on_total_clusters.pdf")
)

plot_selected_proteins(
  total_object,
  assay_matrices$total$values,
  "clusters_total",
  "umap_total",
  "total on total-proteome clusters",
  file.path(path, "Figures/Fig5d_SelectedMarkers_total_on_total_clusters.pdf")
)



# ---- Figure 5e ----
# Display M/L vs M/H values per cell
MH_values_long <- MH_values %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "protein") %>%
  tidyr::pivot_longer(
    cols = -protein,
    names_to = "Cell",
    values_to = "MH"
  )

ML_values_long <- ML_values %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "protein") %>%
  tidyr::pivot_longer(
    cols = -protein,
    names_to = "Cell",
    values_to = "ML"
  )

total_values_long <- total_values %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "protein") %>%
  tidyr::pivot_longer(
    cols = -protein,
    names_to = "Cell",
    values_to = "total"
  )

ML_MH_long <- merge(ML_values_long, MH_values_long, by = c("protein", "Cell"), all = TRUE)
ML_total_long <- merge(ML_values_long, total_values_long, by = c("protein", "Cell"), all = TRUE)
MH_total_long <- merge(MH_values_long, total_values_long, by = c("protein", "Cell"), all = TRUE)

ML_MH_long_plot <- ML_MH_long %>% filter(protein %in% selected_markers)
ML_total_long_plot <- ML_total_long %>% filter(protein %in% selected_markers)
MH_total_long_plot <- MH_total_long %>% filter(protein %in% selected_markers)

clusters <- assay_object@meta.data
clusters <- clusters[,c(4,5)]

plots_data <- list(
  ML = ML_total_long_plot,
  MH = MH_total_long_plot
)

for (label in names(plots_data)) {
  data_plot <- merge(plots_data[[label]], clusters, by.x = "Cell", by.y = "Well")

   stats_df <- data_plot %>%
    group_by(protein) %>%
    group_modify(~{
      model <- lm(.x[[label]] ~ .x$total)

      tibble(
        intercept = coef(model)[1],
        slope = coef(model)[2],
        resid_sd = sd(residuals(model))
      )
    })


  data_plot$clusters_total <- factor(data_plot$clusters_total,levels = c(0, 1, 2))

  p <- ggplot(data_plot, aes(x = total, y = .data[[label]])) +
    geom_point(aes(color = clusters_total)) +
    scale_color_manual(values = c("0" = "#00AEEF","1" = "#007169","2" = "#EB008B")) +
    #geom_pointdensity(adjust = 0.1) + 
    #scale_color_viridis_c() +
    geom_smooth(method = "lm", color = "black", se = FALSE, size = 1) + # Main Regression Line
    geom_abline(data = stats_df,                                        # Add the +2 SD line
                aes(intercept = intercept + (2 * resid_sd), slope = slope), 
                color = "red", linetype = "dashed") +
    geom_abline(data = stats_df,                                        # Add the -2 SD line
                aes(intercept = intercept - (2 * resid_sd), slope = slope), 
                color = "red", linetype = "dashed") +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
    theme_minimal() +
    labs(
      title = paste("log10 total vs log10", label, "with Outlier Boundaries"),
      subtitle = "Red dashed lines represent ±2 Standard Deviations of the residuals",
      x = "Log10 total",
      y = paste("Log10", label),
      color = "Cluster"
    ) +
    facet_grid(~protein)

  pdf(paste0(path, "Figures/Fig5f_Correlation_log10", label, "_log10total_SelectedMarkers_SD2_colClusters.pdf", sep = ""))
  print(p)
  dev.off()
}


cat(
  "\nScript completed successfully: ",
  format(Sys.time()),
  "\n",
  file = log_file,
  append = TRUE,
  sep = ""
)


# ---- Methodological note ----
# proDA was developed for log-intensity proteomics data with intensity-dependent
# dropout. Its assumptions are most direct for total-intensity data. Missingness
# in M/L values can additionally depend on whether either constituent channel
# was quantified; ratio-derived results should therefore be interpreted with care.

#########################################################################################################################################################################
# END OF SCRIPT
#########################################################################################################################################################################