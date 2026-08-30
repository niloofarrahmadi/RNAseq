#!/usr/bin/env bash

set -euo pipefail

PROJECT_DIR="/mnt/d/RNAseq"
FASTQ_DIR="${PROJECT_DIR}/fastq"
QC_DIR="${PROJECT_DIR}/qc"

mkdir -p "${QC_DIR}"

echo "=========================================="
echo "Step 01: FastQC"
echo "=========================================="

for fq in "${FASTQ_DIR}"/*.fastq.gz; do
    echo "Running FastQC on: ${fq}"

    fastqc \
        "${fq}" \
        --outdir "${QC_DIR}"
done

echo
echo "FastQC completed successfully."
echo "Output directory: ${QC_DIR}"
