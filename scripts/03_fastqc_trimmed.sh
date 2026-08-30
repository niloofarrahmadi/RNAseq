#!/usr/bin/env bash

set -euo pipefail

PROJECT_DIR="/mnt/d/RNAseq"
TRIMMED_DIR="${PROJECT_DIR}/trimmed"
QC_DIR="${PROJECT_DIR}/qc_trimmed"

mkdir -p "${QC_DIR}"

echo "=========================================="
echo "Step 03: FastQC after trimming"
echo "=========================================="

for fq in "${TRIMMED_DIR}"/*.trimmed.fastq.gz; do

    echo "Running FastQC on: ${fq}"

    fastqc \
        "${fq}" \
        --outdir "${QC_DIR}"

done

echo
echo "Post-trimming FastQC completed successfully."
echo "Output directory: ${QC_DIR}"
