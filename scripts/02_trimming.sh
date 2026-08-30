#!/usr/bin/env bash

set -euo pipefail

PROJECT_DIR="/mnt/d/RNAseq"
FASTQ_DIR="${PROJECT_DIR}/fastq"
TRIMMED_DIR="${PROJECT_DIR}/trimmed"

mkdir -p "${TRIMMED_DIR}"

ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"

SAMPLES=(
    SRR031708
    SRR031712
    SRR031713
    SRR031720
    SRR031721
    SRR031722
)

echo "=========================================="
echo "Step 02: Adapter trimming with cutadapt"
echo "=========================================="

for srr in "${SAMPLES[@]}"; do

    echo
    echo "========== ${srr} =========="

    cutadapt \
        -a "${ADAPTER}" \
        -m 20 \
        -o "${TRIMMED_DIR}/${srr}.trimmed.fastq.gz" \
        "${FASTQ_DIR}/${srr}.fastq.gz"

done

echo
echo "Trimming completed successfully."
echo "Output directory: ${TRIMMED_DIR}"
