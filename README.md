# singlecellRNA-seq
This repository provides **code for analyzing single-cell RNA sequencing (scRNA-seq) data**. It outlines the steps to begin with a **FASTQ file** and perform **fundamental scRNA-seq analysis**, offering a structured workflow for processing, quality control, and downstream interpretation.

1️⃣ Preprocessing with Cell Ranger – This step processes raw sequencing data, performing alignment, barcode assignment, and UMI counting to generate an expression matrix. (barcodes.tsv, genes.tsv and matrix.mtx)
2️⃣ Velocity Analysis – The output from Cell Ranger (.bam) is used for RNA velocity estimation with Velocyte and (.loop) scVelo, producing two velocity images that illustrate directionality and confidence in cellular transitions. 
3️⃣ Co-expression Analysis – Gene co-expression patterns are analyzed using CS_CORE, providing insights into regulatory interactions among genes. gives an output of .csv file and 3 plots.

This pipeline streamlines scRNA-seq workflows by integrating essential tools for data processing, velocity inference, and co-expression analysis in an automated manner. 

Please change the directories before using according to yours. the outputs will be saved to the same input directories. 

After running the main pipeline, if you've processed your data using Seurat, you'll have an .rds file containing the integrated single-cell dataset. This serves as the input for the next analysis: Differential Expression (DEG) using MAST. deg.R is the script name. 
🔧 DEG Analysis Workflow

1️⃣ Load the .rds File – Import your Seurat object to extract relevant metadata and gene expression data. 2️⃣ Run MAST for DEG – Identify differentially expressed genes using MAST, which will generate a .csv file containing gene details and statistical significance. 3️⃣ Visualize Results – MAST produces four key visualizations:

    UMAP plot – Displays cell clusters and expression patterns.

    Violin plot – Shows expression distributions for significant genes.

    Heatmap – Highlights differentially expressed genes across conditions.

    Feature plot – Maps gene expression onto clusters for interpretation.

This analysis enables a deep dive into cell-type-specific expression patterns, offering critical insights into biological mechanisms.

The final step in this pipeline is Gene Regulatory Network (GRN) analysis, performed using GENIE3, which models transcriptional regulatory relationships.
🔧 GRN Analysis Workflow

1️⃣ Input Preparation

    The expression matrix (.csv) serves as the input, containing gene expression values across single cells.

2️⃣ Running GENIE3 (Computationally Intensive)

    GENIE3 identifies potential gene-gene regulatory interactions using random forest regression.

    This process is time-consuming, averaging 24 hours on a 32-core CPU.

3️⃣ Output Files

    A .csv file containing inferred regulatory links between transcription factors and target genes.

    A visualization (network image) representing the GRN structure, highlighting key regulators.

This final step provides deep insights into transcriptional control mechanisms, uncovering regulatory relationships that drive biological processes.
