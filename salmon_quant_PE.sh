#!/bin/sh

# Abundance quantification for paired-end RNA-seq

for r1 in PE_trimmed/*1_val_1.fq.gz; do
	r2=${r1/1_val_1.fq.gz/2_val_2.fq.gz}
	sample=$(basename "$r1" _1_val_1.fq.gz)

	salmon quant -i salmon_hg38_index -l A -1 "$r1" -2 "$r2" --validateMappings -o salmon_quant/"$sample"
done