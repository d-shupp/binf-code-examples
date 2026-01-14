#!/usr/bin/env bash

SRAS="$1"
SRABIN="/Users/dshupp/Documents/Bioinformatics/sratoolkit.3.2.1-mac-arm64/bin"

while read -r i; do
    echo "Processing $i"
    ${SRABIN}/fastq-dump --split-3 --gzip $i
done < "$SRAS"