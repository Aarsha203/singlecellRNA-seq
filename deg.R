# Load necessary packages
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("MAST", quietly = TRUE)) BiocManager::install("MAST")

library(Seurat)
library(MAST)

# Load Seurat object
seurat_obj <- readRDS("/media/user/8882652682651A48/scRNA/seurat_obj.rds")

# Preprocessing (if not already done)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Plot and save UMAP
png("/media/user/8882652682651A48/scRNA/umap_clusters.png", width = 1000, height = 800)
DimPlot(seurat_obj, reduction = "umap", label = TRUE)
dev.off()

# Differential Expression with MAST (cluster 0 vs all)
markers_cluster0 <- FindMarkers(seurat_obj, ident.1 = 0, test.use = "MAST")

# Save marker results
write.csv(markers_cluster0, "/media/user/8882652682651A48/scRNA/cluster0_markers_MAST.csv")

# Visualize top 5 marker genes
top5 <- head(rownames(markers_cluster0), 5)

# Violin plot
png("/media/user/8882652682651A48/scRNA/violin_top5.png", width = 1200, height = 800)
VlnPlot(seurat_obj, features = top5, pt.size = 0.1)
dev.off()

# FeaturePlot (on UMAP)
png("/media/user/8882652682651A48/scRNA/feature_top5.png", width = 1200, height = 800)
FeaturePlot(seurat_obj, features = top5)
dev.off()

# Heatmap of top 10 genes
top10 <- head(rownames(markers_cluster0), 10)
png("/media/user/8882652682651A48/scRNA/heatmap_top10.png", width = 1200, height = 800)
DoHeatmap(seurat_obj, features = top10)
dev.off()

cat("✅ Analysis complete. Results and plots saved to /media/user/8882652682651A48/scRNA/\n")
