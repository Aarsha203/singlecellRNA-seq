# singlecellRNA-seq
This repository provides **code for analyzing single-cell RNA sequencing (scRNA-seq) data**. It outlines the steps to begin with a **FASTQ file** and perform **fundamental scRNA-seq analysis**, offering a structured workflow for processing, quality control, and downstream interpretation.

1️⃣ Preprocessing with Cell Ranger – This step processes raw sequencing data, performing alignment, barcode assignment, and UMI counting to generate an expression matrix. (barcodes.tsv, genes.tsv and matrix.mtx)
2️⃣ Velocity Analysis – The output from Cell Ranger (.bam) is used for RNA velocity estimation with Velocyte and (.loop) scVelo, producing two velocity images that illustrate directionality and confidence in cellular transitions. 
3️⃣ Co-expression Analysis – Gene co-expression patterns are analyzed using CS_CORE, providing insights into regulatory interactions among genes. gives an output of .csv file and 3 plots.

This pipeline streamlines scRNA-seq workflows by integrating essential tools for data processing, velocity inference, and co-expression analysis in an automated manner. 

Please change the directories before using according to yours. the outputs will be saved to the same input directories. 
