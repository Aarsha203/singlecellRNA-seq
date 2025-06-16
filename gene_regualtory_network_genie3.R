library(GENIE3)
library(matrixStats)

# Read the matrix (genes as columns, cells as rows)
exprMatrix <- read.csv("/media/user/8882652682651A48/scRNA/expression_matrix_final.csv", row.names = 1)

# Transpose so that rows are genes
exprMatrix <- t(as.matrix(exprMatrix))

# Filter for top 2000 most variable genes
gene_vars <- rowVars(exprMatrix)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:2000]
exprMatrix <- exprMatrix[top_genes, ]
cat("✅ Filtered top 2000 variable genes. Matrix now has", nrow(exprMatrix), "genes.\n")

# Load TF list
TFs <- scan("/media/user/8882652682651A48/scRNA/tfs.txt", what = "character")

# Filter TFs to those present in the expression matrix
common_TFs <- intersect(TFs, rownames(exprMatrix))
cat("✅ Found", length(common_TFs), "TFs in expression matrix.\n")

# Run GENIE3 if sufficient TFs are found
if (length(common_TFs) >= 2) {
  set.seed(123)
  weightMatrix <- GENIE3(exprMatrix, regulators = common_TFs, nCores = 10, verbose = TRUE)
  linkList <- getLinkList(weightMatrix)
  write.csv(linkList, file = "/media/user/8882652682651A48/scRNA/genie3_output_links.csv", row.names = FALSE)
  cat("✅ GENIE3 network inference complete. Output saved to genie3_output_links.csv\n")
} else {
  stop("❌ Less than 2 TFs found in expression matrix. Please check gene names.")
}

#plot image
library(cowplot)

# Make your plot without title
p <- ggraph(tg, layout = "fr") +
  geom_edge_link(aes(edge_alpha = weight), 
                 arrow = arrow(length = unit(4, 'mm')), 
                 end_cap = circle(3, 'mm'),
                 edge_colour = "gray50") +
  geom_node_point(size = 5, aes(color = name %in% top_regulators)) +
  geom_node_text(aes(label = name), size = 4, family = "serif", repel = TRUE, max.overlaps = Inf) +
  scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#0072B2"),
                     labels = c("Target", "Regulator"),
                     name = "Node Type") +
  scale_edge_alpha(range = c(0.3, 1)) +
  theme_void(base_family = "serif") +
  theme(
    legend.position = "right",
    plot.margin = margin(20, 40, 20, 20)
  )

# Add title and underline using cowplot
title <- ggdraw() +
  draw_label("Gene Regulatory Network", fontface = 'bold', size = 18, x = 0.5, hjust = 0.5) +
  geom_hline(yintercept = 0.1, size = 0.8)

plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1))
# Save the plot
output_path <- "/media/user/8882652682651A48/scRNA/genie3_regulatory_network.png"

ggsave(output_path, plot = plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1)),
       width = 10, height = 8, dpi = 300)

cat("✅ PNG saved to:", output_path, "\n")

