#!/usr/bin/env bash

set -euo pipefail

PROJECT="/mnt/d/RNAseq"
ALIGN_DIR="${PROJECT}/alignment"

for sample in \
    SRR031708 \
    SRR031712 \
    SRR031713 \
    SRR031720 \
    SRR031721 \
    SRR031722
do
    echo "=========================================="
    echo "Processing ${sample}"
    echo "=========================================="

    samtools sort \
        -o "${ALIGN_DIR}/${sample}.sorted.bam" \
        "${ALIGN_DIR}/${sample}.sam"

    samtools index \
        "${ALIGN_DIR}/${sample}.sorted.bam"

    samtools quickcheck -v \
        "${ALIGN_DIR}/${sample}.sorted.bam"

    samtools flagstat \
        "${ALIGN_DIR}/${sample}.sorted.bam" \
        > "${ALIGN_DIR}/logs/${sample}.flagstat.txt"

    echo "${sample} BAM processing completed."
done

echo "All BAM files processed and checked."
