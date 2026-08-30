#!/usr/bin/env bash

set -euo pipefail

PROJECT="/mnt/d/RNAseq"

TRIMMED_DIR="${PROJECT}/trimmed"
ALIGN_DIR="${PROJECT}/alignment"
LOG_DIR="${ALIGN_DIR}/logs"

INDEX="${PROJECT}/reference/BDGP5.25.62/Drosophila_melanogaster.BDGP5.25.62"

mkdir -p "${ALIGN_DIR}" "${LOG_DIR}"

SAMPLES=(
    SRR031708
    SRR031712
    SRR031713
    SRR031720
    SRR031721
    SRR031722
)

echo "=========================================="
echo "Step 04: HISAT2 alignment"
echo "Reference: Drosophila melanogaster BDGP5.25.62"
echo "=========================================="

for sample in "${SAMPLES[@]}"; do

    echo
    echo "========== ${sample} =========="

    hisat2 \
        -x "${INDEX}" \
        -U "${TRIMMED_DIR}/${sample}.trimmed.fastq.gz" \
        -S "${ALIGN_DIR}/${sample}.sam" \
        2> "${LOG_DIR}/${sample}.hisat2.log"

    echo "${sample} alignment completed."

done

echo
echo "All HISAT2 alignments completed."
