#!/usr/bin/env bash
set -euo pipefail

###############################################################################
#  Viral contig abundance pipeline (coverage-based)
#
#  Author: Hayde Saracho
#  Purpose:
#    1. Build Bowtie2 indexes for viral contigs from each mat type
#    2. Map each replicate independently against the corresponding viral contigs
#    3. Generate sorted/indexed BAM files
#    4. Compute per-base depth including zero-coverage positions
#    5. Compute mean depth per contig
#
#  Input FASTQ files (interleaved paired-end):
#    HMAb_1 -> H3.fastq
#    HMAb_2 -> HMAb.fastq
#    HMNa_1 -> H4.fastq
#    HMNa_2 -> HMNa.fastq
#    HMR_1  -> H8.fastq
#    HMR_2  -> HMR.fastq
#
#  Input viral contig FASTA files:
#    HMAb -> Green_Mat_Viruses.fasta
#    HMNa -> Black_Mat_Viruses.fasta
#    HMR  -> Pink_Mat_Viruses.fasta
#
#  Output:
#    bowtie2_index/
#    bam/
#    depth/
#    coverage/
###############################################################################

# ----------------------------- User-configurable ----------------------------- #

THREADS=8

# Viral contig FASTA files
HMAB_FASTA="Green_Mat_Viruses.fasta"
HMNA_FASTA="Black_Mat_Viruses.fasta"
HMR_FASTA="Pink_Mat_Viruses.fasta"

# Interleaved FASTQ files
HMAB_R1_FASTQ="H3.fastq"
HMAB_R2_FASTQ="HMAb.fastq"

HMNA_R1_FASTQ="H4.fastq"
HMNA_R2_FASTQ="HMNa.fastq"

HMR_R1_FASTQ="H8.fastq"
HMR_R2_FASTQ="HMR.fastq"

# Output directories
INDEX_DIR="bowtie2_index"
BAM_DIR="bam"
DEPTH_DIR="depth"
COV_DIR="coverage"

# --------------------------------- Checks ---------------------------------- #

echo "Checking required software..."
command -v bowtie2 >/dev/null 2>&1 || { echo "Error: bowtie2 not found in PATH"; exit 1; }
command -v bowtie2-build >/dev/null 2>&1 || { echo "Error: bowtie2-build not found in PATH"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "Error: samtools not found in PATH"; exit 1; }
command -v awk >/dev/null 2>&1 || { echo "Error: awk not found in PATH"; exit 1; }

echo "Checking input files..."
for f in \
    "$HMAB_FASTA" "$HMNA_FASTA" "$HMR_FASTA" \
    "$HMAB_R1_FASTQ" "$HMAB_R2_FASTQ" \
    "$HMNA_R1_FASTQ" "$HMNA_R2_FASTQ" \
    "$HMR_R1_FASTQ" "$HMR_R2_FASTQ"
do
    [[ -f "$f" ]] || { echo "Error: input file not found -> $f"; exit 1; }
done

mkdir -p "$INDEX_DIR" "$BAM_DIR" "$DEPTH_DIR" "$COV_DIR"

# ------------------------------ Build indexes ------------------------------- #

echo "Building Bowtie2 indexes..."
bowtie2-build "$HMAB_FASTA" "${INDEX_DIR}/HMAb"
bowtie2-build "$HMNA_FASTA" "${INDEX_DIR}/HMNa"
bowtie2-build "$HMR_FASTA"  "${INDEX_DIR}/HMR"

# ----------------------------- Sample metadata ------------------------------ #

declare -A INDEX_PREFIX
declare -A FASTQ_FILE

INDEX_PREFIX["HMAb_1"]="${INDEX_DIR}/HMAb"
INDEX_PREFIX["HMAb_2"]="${INDEX_DIR}/HMAb"
INDEX_PREFIX["HMNa_1"]="${INDEX_DIR}/HMNa"
INDEX_PREFIX["HMNa_2"]="${INDEX_DIR}/HMNa"
INDEX_PREFIX["HMR_1"]="${INDEX_DIR}/HMR"
INDEX_PREFIX["HMR_2"]="${INDEX_DIR}/HMR"

FASTQ_FILE["HMAb_1"]="$HMAB_R1_FASTQ"
FASTQ_FILE["HMAb_2"]="$HMAB_R2_FASTQ"
FASTQ_FILE["HMNa_1"]="$HMNA_R1_FASTQ"
FASTQ_FILE["HMNa_2"]="$HMNA_R2_FASTQ"
FASTQ_FILE["HMR_1"]="$HMR_R1_FASTQ"
FASTQ_FILE["HMR_2"]="$HMR_R2_FASTQ"

SAMPLES=("HMAb_1" "HMAb_2" "HMNa_1" "HMNa_2" "HMR_1" "HMR_2")

# ------------------------------- Main loop --------------------------------- #

for sample in "${SAMPLES[@]}"; do
    echo
    echo "=============================="
    echo "Processing sample: $sample"
    echo "FASTQ:  ${FASTQ_FILE[$sample]}"
    echo "Index:  ${INDEX_PREFIX[$sample]}"
    echo "=============================="

    BAM_OUT="${BAM_DIR}/${sample}.sorted.bam"
    DEPTH_OUT="${DEPTH_DIR}/${sample}.depth.txt"
    MEAN_OUT="${COV_DIR}/${sample}.mean_depth.txt"

    # 1. Map interleaved paired-end reads and sort BAM
    echo "[1/4] Mapping with Bowtie2..."
    bowtie2 \
        --interleaved "${FASTQ_FILE[$sample]}" \
        -x "${INDEX_PREFIX[$sample]}" \
        -p "$THREADS" | \
    samtools view -bS - | \
    samtools sort -@ "$THREADS" -o "$BAM_OUT"

    # 2. Index BAM
    echo "[2/4] Indexing BAM..."
    samtools index "$BAM_OUT"

    # 3. Compute depth including zero-coverage positions
    echo "[3/4] Calculating per-base depth with zero-covered positions included..."
    samtools depth -a "$BAM_OUT" > "$DEPTH_OUT"

    # 4. Compute mean depth per contig
    echo "[4/4] Calculating mean depth per contig..."
    awk '
        {
            sum[$1] += $3;
            count[$1] += 1;
        }
        END {
            for (contig in sum) {
                mean = sum[contig] / count[contig];
                print contig "\t" mean;
            }
        }
    ' "$DEPTH_OUT" | sort -k1,1 > "$MEAN_OUT"

    echo "Done: $sample"
    echo "  BAM:       $BAM_OUT"
    echo "  DEPTH:     $DEPTH_OUT"
    echo "  MEAN DEPTH:$MEAN_OUT"
done

echo
echo "Pipeline finished successfully."
echo "Mean depth files are in: ${COV_DIR}/"