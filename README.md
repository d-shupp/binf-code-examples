# binf-code-examples
Bash, R, and Python script examples from RNA-seq and machine learning pipeline.

## download_sras.sh
For loop using sra toolkit's fastq-dump to download list of samples, list provided as argument.

## salmon_quant_PE.sh
Abundance quantification for paired-ended RNA-seq data. Loops over all fastq files, records sample names to use when generating quant output.

## batch_correction_deseq2.R
Uses ComBat_seq to perform batch correction on RNA-seq gene abundance files. Runs DESeq2 on corrected and non-corrected counts to generate PCA plot, Volcano plot, and Distance Heatmap.

## ml.ipynb
Machine learning script. Takes batch-correced, variance-stablized gene abudance matrix and runs two feature selection algorithms (logistic regression with L1 penatly, random forest) to select important gene features. These feature sets are then used to train three models (Logistic regression with L1 penalty, random forest, and support vector machine) and tune parameters. Model accuracy is displayed.
