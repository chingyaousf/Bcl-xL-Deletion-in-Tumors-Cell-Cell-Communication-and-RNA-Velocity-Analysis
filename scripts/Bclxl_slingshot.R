# ====================================================
# Slingshot Trajectory Analysis on Bcl-xL Seurat Object
# ====================================================

BiocManager::install("slingshot")

# Load libraries
library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(scater)


# ---------------------------------------------
# Step 1: Load the annotated Seurat object
# ---------------------------------------------
seurat_obj <- readRDS("/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/bclxl_annotated_noRBC.rds")

head(seurat_obj)

# ---------------------------------------------
# Step 2: Convert Seurat to SingleCellExperiment
# ---------------------------------------------
sce <- as.SingleCellExperiment(seurat_obj)

# ---------------------------------------------
# Step 3: Run Slingshot using UMAP and clustering
# ---------------------------------------------
# If you want to use another label (e.g., cell type), change `seurat_clusters` below
sce <- slingshot(
  sce,
  clusterLabels = sce$seurat_clusters,
  reducedDim = "UMAP"
)

# ---------------------------------------------
# Step 4: Plot UMAP with inferred trajectories
# ---------------------------------------------
# Save trajectory plot as PNG
png(
  filename = "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/slingshot/slingshot_trajectory_umap.png",
  width = 7, height = 7, units = "in", res = 300
)

plot(
  reducedDims(sce)$UMAP,
  col = sce$seurat_clusters,
  pch = 16,
  cex = 0.3,  # <<-- THIS makes points smaller
  main = "Slingshot Trajectory on UMAP"
)
lines(SlingshotDataSet(sce), lwd = 2)

dev.off()


# ---------------------------------------------
# Step 5: Extract and save pseudotime values
# ---------------------------------------------
pseudotime <- slingPseudotime(sce)
head(pseudotime)

# Save pseudotime to Seurat object (first lineage)
seurat_obj$pseudotime <- pseudotime[, 1]

# Save the updated Seurat object with pseudotime added
saveRDS(
  seurat_obj,
  file = "/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/slingshot/bclxl_with_pseudotime.rds"
)



## Plot Pseudotime for individual Lineage
# Extract UMAP coordinates and pseudotime
umap_df <- as.data.frame(reducedDims(sce)$UMAP)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$pseudotime <- slingPseudotime(sce)[, 1]  # First lineage

# Plot with pseudotime gradient
library(ggplot2)

png("/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/slingshot/umap_pseudotime_gradient_first_lineage.png",
    width = 6, height = 5, units = "in", res = 300)

ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c(option = "plasma", na.value = "grey90") +
  labs(title = "UMAP Colored by Slingshot Pseudotime First lineage") +
  theme_classic()

dev.off()


## Plot Pseudotime for combined Lineage
# ------------------------------------------------
# Extract UMAP coordinates and combined pseudotime
# ------------------------------------------------
umap_df <- as.data.frame(reducedDims(sce)$UMAP)
colnames(umap_df) <- c("UMAP_1", "UMAP_2")

# Combine pseudotime across all lineages (minimum pseudotime per cell)
pseudotime_all <- slingPseudotime(sce)
umap_df$pseudotime <- apply(pseudotime_all, 1, function(x) min(x, na.rm = TRUE))
umap_df$pseudotime[!is.finite(umap_df$pseudotime)] <- NA  # Clean up

# ------------------------------------------------
# Plot UMAP with pseudotime gradient
# ------------------------------------------------
library(ggplot2)

png("/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/slingshot/umap_combined_pseudotime_gradient.png",
    width = 6, height = 5, units = "in", res = 300)

ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c(option = "plasma", na.value = "grey90") +
  labs(title = "UMAP Colored by Combined Pseudotime (All Lineages)") +
  theme_classic()

dev.off()


## Code to Plot individual Genes Across Pseudotime
## Loop to Plot Pseudotime for Each Lineage
# Load required libraries
library(ggplot2)
library(viridis)

# Get UMAP coordinates
umap_df_base <- as.data.frame(reducedDims(sce)$UMAP)
colnames(umap_df_base) <- c("UMAP_1", "UMAP_2")

# Extract all pseudotime values
pseudotime_all <- slingPseudotime(sce)

# Loop through each lineage
for (i in 1:ncol(pseudotime_all)) {
  
  # Clone base UMAP data and add pseudotime for this lineage
  umap_df <- umap_df_base
  umap_df$pseudotime <- pseudotime_all[, i]
  
  # Clean up NAs (those not assigned to this lineage)
  umap_df$pseudotime[!is.finite(umap_df$pseudotime)] <- NA
  
  # Set output filename
  filename <- paste0("/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/slingshot/umap_pseudotime_gradient_lineage_", i, ".png")
  
  # Save plot
  png(filename, width = 6, height = 5, units = "in", res = 300)
  
  print(
    ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
      geom_point(size = 0.5) +
      scale_color_viridis_c(option = "plasma", na.value = "grey90") +
      labs(title = paste("UMAP Colored by Slingshot Pseudotime - Lineage", i)) +
      theme_classic()
  )
  
  dev.off()
}


## Loop Code to Plot Multiple Genes Across Pseudotime
# Add gene expression to the UMAP data
#umap_df$gene_expr <- seurat_obj[["RNA"]]@data["Mki67", ]
#umap_df$gene_expr <- seurat_obj[["RNA"]]@data["App", ]
#umap_df$gene_expr <- seurat_obj[["RNA"]]@data["Nrcam", ]
umap_df$gene_expr <- seurat_obj[["RNA"]]@data["Gja1", ]

# Plot gene expression across pseudotime
#png("/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/slingshot/gene_Mki67_vs_pseudotime.png",
#png("/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/slingshot/gene_App_vs_pseudotime.png",
#png("/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/slingshot/gene_Nrcam_vs_pseudotime.png",
png("/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/slingshot/gene_Gja1_vs_pseudotime.png",
    width = 6, height = 4, units = "in", res = 300)

ggplot(umap_df, aes(x = pseudotime, y = gene_expr)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_smooth(method = "loess", color = "firebrick") +
  #labs(title = "Mki67 Expression over Pseudotime",
  #labs(title = "App Expression over Pseudotime",
  #labs(title = "Nrcam Expression over Pseudotime",
  labs(title = "Gja1 Expression over Pseudotime",
       x = "Pseudotime", y = "Log-normalized Expression") +
  theme_minimal()

dev.off()



# Full list of genes to include (from previous and new)
genes_to_plot <- c("Mki67", "App", "Nrcam", "Gja1", "Ncam1", "Pcdhga9", "Psap", "Cadm1", "Bsg")

# Loop through each gene and generate a pseudotime plot
for (gene in genes_to_plot) {
  
  # Add expression to UMAP dataframe
  umap_df$gene_expr <- seurat_obj[["RNA"]]@data[gene, ]
  
  # Define output file path
  output_file <- paste0("/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/slingshot/gene_", gene, "_vs_pseudotime.png")
  
  # Save the plot
  png(output_file, width = 6, height = 4, units = "in", res = 300)
  
  # Plot expression vs pseudotime
  print(
    ggplot(umap_df, aes(x = pseudotime, y = gene_expr)) +
      geom_point(size = 0.5, alpha = 0.5) +
      geom_smooth(method = "loess", color = "firebrick") +
      labs(title = paste(gene, "Expression over Pseudotime"),
           x = "Pseudotime", y = "Log-normalized Expression") +
      theme_minimal()
  )
  
  dev.off()
}




# ============================================================
# Identify Branching Points and Diverging Fates
# ============================================================

# ----- Step 1: Check number of inferred lineages -----
n_lineages <- length(slingCurves(sce))
print(paste("Number of inferred lineages:", n_lineages))

# ----- Step 2: Plot UMAP with Slingshot trajectories -----
png("/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/slingshot/umap_with_trajectories.png",
    width = 6, height = 5, units = "in", res = 300)

plot(reducedDims(sce)$UMAP, col = as.factor(seurat_obj$cell_type), pch = 16, cex = 0.3,
     main = "UMAP with Slingshot Trajectories")
lines(SlingshotDataSet(sce), lwd = 2)

dev.off()

# ----- Step 3: Identify branching cells from lineage weights -----
sling_lineage_weights <- slingCurveWeights(sce)

# Cells with significant weight in more than one lineage = likely branching points
branching_cells <- which(rowSums(sling_lineage_weights > 0.25) > 1)

# Mark cells as branching or not
umap_df$branching <- ifelse(rownames(umap_df) %in% names(branching_cells), "Yes", "No")

# ----- Step 4: Plot UMAP highlighting branching points -----
png("/home/cyang40/chingyao/BclXL_PHGDH_project/BclXL/slingshot/umap_branching_cells.png",
    width = 6, height = 5, units = "in", res = 300)

ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = branching)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("Yes" = "purple", "No" = "lightgray")) +
  labs(title = "Branching Cells in UMAP") +
  theme_classic() +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
      ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

dev.off()



# ================================================================
# tradeSeq Lineage-Specific Gene Expression Analysis with Slingshot
# ================================================================

# Install tradeSeq from Bioconductor
BiocManager::install("tradeSeq")

# ---------------------------------------------
# Step 1: Load required packages
# ---------------------------------------------
library(SingleCellExperiment)
library(slingshot)
library(tradeSeq)
library(mgcv)

# =======================================================
# Step 2: Prepare Input Matrices (Remove NA Pseudotime)
# =======================================================

# Extract matrices from Slingshot-fitted SingleCellExperiment
counts <- as.matrix(assays(sce)$counts)
pseudotime <- slingPseudotime(sce)
cellWeights <- slingCurveWeights(sce)

# Step 2.1: Keep only cells with pseudotime defined for all lineages (no NAs)
valid_cells <- rowSums(is.na(pseudotime)) == 0

# Step 2.2: Subset matrices to valid cells
counts <- counts[, valid_cells]
pseudotime <- pseudotime[valid_cells, ]
cellWeights <- cellWeights[valid_cells, ]

# Step 2.3: Sanity check
stopifnot(all(colnames(counts) == rownames(pseudotime)))

# ---------------------------------------------
# Step 3: Fit GAM models across pseudotime
# ---------------------------------------------
set.seed(123)
gam_fit <- fitGAM(counts = counts,
                  pseudotime = pseudotime,
                  cellWeights = cellWeights,
                  nknots = 6,
                  verbose = TRUE)

# ---------------------------------------------
# Step 4: Association Test — Genes changing along lineages
# ---------------------------------------------
assoRes <- associationTest(gam_fit)

# Extract significant genes
sig_genes <- rownames(assoRes)[which(assoRes$padj < 0.05)]
write.csv(assoRes, "associationTest_results.csv")

# ---------------------------------------------
# Step 5: Fate-Specific Expression — End-point comparison
# ---------------------------------------------
diffEnd <- diffEndTest(gam_fit)

# Genes diverging between lineages
fate_genes <- rownames(diffEnd)[which(diffEnd$padj < 0.05)]
write.csv(diffEnd, "diffEndTest_fateGenes.csv")

# ---------------------------------------------
# Step 6: Plot smoothed expression for top genes
# ---------------------------------------------
top_genes <- head(rownames(diffEnd[order(diffEnd$padj), ]), 10)

# Create plots directory
dir.create("slingshot/tradeSeq_genePlots", showWarnings = FALSE, recursive = TRUE)

for (gene in top_genes) {
  png(paste0("slingshot/tradeSeq_genePlots/", gene, "_smoother.png"),
      width = 6, height = 4, units = "in", res = 300)
  plotSmoothers(gam_fit, counts, gene = gene)
  dev.off()
}
