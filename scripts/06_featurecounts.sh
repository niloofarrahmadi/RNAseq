#!/usr/bin/env bash

set -euo pipefail

PROJECT="/mnt/d/RNAseq"

ANNOTATION="${PROJECT}/reference/BDGP5.25.62/Drosophila_melanogaster.BDGP5.25.62.gtf"
ALIGN_DIR="${PROJECT}/alignment"
OUTPUT_DIR="${PROJECT}/results/counts"

mkdir -p "${OUTPUT_DIR}"

featureCounts \
    -T 4 \
    -t exon \
    -g gene_id \
    -a "${ANNOTATION}" \
    -o "${OUTPUT_DIR}/all_samples_gene_counts.txt" \
    "${ALIGN_DIR}/SRR031708.sorted.bam" \
    "${ALIGN_DIR}/SRR031712.sorted.bam" \
    "${ALIGN_DIR}/SRR031713.sorted.bam" \
    "${ALIGN_DIR}/SRR031720.sorted.bam" \
    "${ALIGN_DIR}/SRR031721.sorted.bam" \
    "${ALIGN_DIR}/SRR031722.sorted.bam"

echo "featureCounts completed."
echo "Output: ${OUTPUT_DIR}/all_samples_gene_counts.txt"
echo "Summary: ${OUTPUT_DIR}/all_samples_gene_counts.txt.summary"
