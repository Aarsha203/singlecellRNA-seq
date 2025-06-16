#!/bin/bash

set -e  # Exit on error

echo "🔹 Step 1: Running Cell Ranger..."

if [ -d "P12_GEX/outs" ]; then
    echo "✅ Cell Ranger output found, skipping..."
else
    cellranger count \
      --id=P12_GEX \
      --transcriptome=/mnt/60975228-7325-434c-8414-1d8464f71392/scRNa/ensembl_reference/refdata-cellranger-GRCh38-3.0.0 \
      --fastqs=/mnt/60975228-7325-434c-8414-1d8464f71392/scRNa \
      --sample=Sample_P12 \
      --localcores=28 \
      --localmem=64 \
      --create-bam=true
    echo "✅ Cell Ranger completed"
fi

echo "🔹 Step 2: Running Velocyto..."

LOOM_OUTPUT="/mnt/60975228-7325-434c-8414-1d8464f71392/scRNa/velocyto_results/possorted_genome_P12_GEX.loom"

if [ -f "$LOOM_OUTPUT" ]; then
    echo "✅ Loom file found, skipping Velocyto..."
else
   velocyto run10x \
  P12_GEX \
  /mnt/60975228-7325-434c-8414-1d8464f71392/scRNa/ensembl_reference/Homo_sapiens.GRCh38.109.gtf

# Move loom file to your expected path
mkdir -p /mnt/60975228-7325-434c-8414-1d8464f71392/scRNa/velocyto_results
mv P12_GEX/velocyto/*.loom /mnt/60975228-7325-434c-8414-1d8464f71392/scRNa/velocyto_results/possorted_genome_P12_GEX.loom

    echo "✅ Velocyto completed, loom file created"
fi

echo "🔹 Step 3: Running scVelo in Python..."

SCVELO_OUT="/mnt/60975228-7325-434c-8414-1d8464f71392/scRNa/cache/mnt-60975228-7325-434c-8414-1d8464f71392-scRNa-velocyto_results-possorted_genome_P12_GEX.h5ad"

if [ -f "$SCVELO_OUT" ]; then
    echo "✅ scVelo output found, skipping Python step..."
else
    # Ensure conda is ready to activate envs
    eval "$(conda shell.bash hook)"
    conda activate scvelo_env

    python3 <<EOF
import scvelo as scv
import scanpy as sc

scv.settings.set_figure_params(dpi=100)
scv.logging.print_version()

adata = scv.read("$LOOM_OUTPUT", cache=True)

scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata)
scv.pp.log1p(adata)
scv.pp.moments(adata)

scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)

sc.pp.neighbors(adata)
sc.tl.leiden(adata)

if 'umap' not in adata.obsm:
    sc.tl.umap(adata)

scv.tl.velocity_confidence(adata)
scv.pl.scatter(adata, color='velocity_confidence', cmap='coolwarm', save='velocity_confidence.png')
scv.pl.velocity_embedding_stream(adata, basis='umap', color='leiden', cmap='viridis', save='velocity_plot.png')

adata.write("$SCVELO_OUT")
print("✅ scVelo analysis completed")
EOF
fi

echo "🔹 Step 4: Running CS-CORE co-expression analysis..."

CSCORE_DIR="/mnt/60975228-7325-434c-8414-1d8464f71392/scRNa/CS_core_coexpression"
CSCORE_OUT="$CSCORE_DIR/cs_core_network.csv"
EXPR_MATRIX="$CSCORE_DIR/expression_matrix_final.csv"

if [ -f "$CSCORE_OUT" ]; then
    echo "✅ CS-CORE output found, skipping CS-CORE step..."
else
    mkdir -p "$CSCORE_DIR"
    python3 <<EOF
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
from CSCORE import CSCORE
import time
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load expression matrix
expr = pd.read_csv("$EXPR_MATRIX", index_col=0)

# Create AnnData object
adata = sc.AnnData(X=csr_matrix(expr.values))
adata.var_names = expr.columns
adata.obs_names = expr.index
adata.raw = adata.copy()

# Run CS-CORE
gene_index = list(range(expr.shape[1]))
start = time.time()
res = CSCORE(adata, gene_index)
end = time.time()

# Save co-expression matrix
adj_df = pd.DataFrame(res[0], index=expr.columns, columns=expr.columns)
adj_df.to_csv("$CSCORE_OUT")
print(f"✅ CS-CORE completed in {end - start:.2f} seconds. Network saved to cs_core_network.csv")

# Top N network visualization
top_n = 100
genes = list(expr.columns[:top_n])
df = adj_df.iloc[:top_n, :top_n]

# Extract top edges
edges = []
for i in range(top_n):
    for j in range(i+1, top_n):
        weight = df.iloc[i, j]
        if pd.notnull(weight) and abs(weight) > 0:
            edges.append((df.index[i], df.columns[j], abs(weight)))
top_edges = sorted(edges, key=lambda x: x[2], reverse=True)[:10]

# Build graph
G = nx.Graph()
G.add_weighted_edges_from(top_edges)

# Plot network
plt.figure(figsize=(10, 8))
pos = nx.spring_layout(G, seed=42)
nx.draw(G, pos, with_labels=True, node_color='lightblue',
        edge_color='gray', node_size=2000, font_size=10, width=2)
plt.title("Top 10 Gene Co-expression Edges")
plt.tight_layout()
plt.savefig("$CSCORE_DIR/top10_network.png")
plt.close()

# Plot heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(df.iloc[:30, :30], cmap="vlag", center=0)
plt.title("Heatmap of Top 30 Highly Co-expressed Genes")
plt.tight_layout()
plt.savefig("$CSCORE_DIR/heatmap_top30.png")
plt.close()

# Plot histogram of co-expression strengths
values = [df.iloc[i, j] for i in range(top_n) for j in range(i+1, top_n)]
plt.figure(figsize=(8, 5))
plt.hist(values, bins=40, color='steelblue', edgecolor='black')
plt.title("Distribution of Co-expression Values (Top 100 Genes)")
plt.xlabel("Co-expression Estimate")
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig("$CSCORE_DIR/coexpression_histogram.png")
plt.close()
EOF
fi

echo "🎉 All steps checked and executed if needed!"

