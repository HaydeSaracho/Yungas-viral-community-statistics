#!/usr/bin/env bash
set -euo pipefail

###############################################################################
#  Viral contig abundance pipeline (coverage-based)
#
#  Author: Hayde Saracho
#  Purpose:
#    1. Build Bowtie2 indexes for viral contigs from each mat type
#    2. Map each replicate independently against the corresponding viral contigs
#    3. Generate sorted/indexed BAM files (mapped reads only)
#    4. Compute per-base depth including zero-coverage positions
#    5. Compute mean depth per contig
#
#  Context:
#    Three metagenomes from microbial mats (hot springs, Yungas rainforest,
#    Jujuy - Argentina). Each metagenome was sequenced twice (two replicates).
#    Reads are interleaved paired-end FASTQ files.
#
#  Input FASTQ files (interleaved paired-end):
#    HMAb replicate 1 -> H3.fastq
#    HMAb replicate 2 -> HMAb.fastq
#    HMNa replicate 1 -> H4.fastq
#    HMNa replicate 2 -> HMNa.fastq
#    HMR  replicate 1 -> H8.fastq
#    HMR  replicate 2 -> HMR.fastq
#
#  Input viral contig FASTA files (one per mat type):
#    HMAb -> Green_Mat_Viruses.fasta
#    HMNa -> Black_Mat_Viruses.fasta
#    HMR  -> Pink_Mat_Viruses.fasta
#
#  Output:
#    bowtie2_index/   Bowtie2 indexes
#    bam/             Sorted, indexed BAM files (mapped reads only)
#    depth/           Per-base depth files
#    data/coverage/   Mean depth per contig (TSV, one file per replicate)
###############################################################################

# ----------------------------- User-configurable ----------------------------- #

THREADS=8

# Viral contig FASTA files (one per mat type)
HMAB_FASTA="Green_Mat_Viruses.fasta"
HMNA_FASTA="Black_Mat_Viruses.fasta"
HMR_FASTA="Pink_Mat_Viruses.fasta"

# Interleaved paired-end FASTQ files (replicate 1 and replicate 2 per mat)
HMAB_REP1_FASTQ="H3.fastq"
HMAB_REP2_FASTQ="HMAb.fastq"

HMNA_REP1_FASTQ="H4.fastq"
HMNA_REP2_FASTQ="HMNa.fastq"

HMR_REP1_FASTQ="H8.fastq"
HMR_REP2_FASTQ="HMR.fastq"

# Output directories
INDEX_DIR="bowtie2_index"
BAM_DIR="bam"
DEPTH_DIR="depth"
COV_DIR="data/coverage"   # matches repository structure

# --------------------------------- Checks ---------------------------------- #

echo "Checking required software..."
for cmd in bowtie2 bowtie2-build samtools awk; do
    command -v "$cmd" >/dev/null 2>&1 || { echo "Error: $cmd not found in PATH"; exit 1; }
done

echo "Checking input files..."
for f in \
    "$HMAB_FASTA" "$HMNA_FASTA" "$HMR_FASTA" \
    "$HMAB_REP1_FASTQ" "$HMAB_REP2_FASTQ" \
    "$HMNA_REP1_FASTQ" "$HMNA_REP2_FASTQ" \
    "$HMR_REP1_FASTQ"  "$HMR_REP2_FASTQ"
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

FASTQ_FILE["HMAb_1"]="$HMAB_REP1_FASTQ"
FASTQ_FILE["HMAb_2"]="$HMAB_REP2_FASTQ"
FASTQ_FILE["HMNa_1"]="$HMNA_REP1_FASTQ"
FASTQ_FILE["HMNa_2"]="$HMNA_REP2_FASTQ"
FASTQ_FILE["HMR_1"]="$HMR_REP1_FASTQ"
FASTQ_FILE["HMR_2"]="$HMR_REP2_FASTQ"

SAMPLES=("HMAb_1" "HMAb_2" "HMNa_1" "HMNa_2" "HMR_1" "HMR_2")

# ------------------------------- Main loop ---------------------------------- #

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

    # 1. Map interleaved paired-end reads, retain only mapped reads, sort BAM
    #    -F 4: exclude unmapped reads to avoid inflating depth calculations
    echo "[1/4] Mapping with Bowtie2 (mapped reads only)..."
    bowtie2 \
        --interleaved "${FASTQ_FILE[$sample]}" \
        -x "${INDEX_PREFIX[$sample]}" \
        -p "$THREADS" \
        2>"${BAM_DIR}/${sample}.bowtie2.log" | \
    samtools view -bS -F 4 - | \
    samtools sort -@ "$THREADS" -o "$BAM_OUT"

    # 2. Index BAM
    echo "[2/4] Indexing BAM..."
    samtools index "$BAM_OUT"

    # Log mapping summary
    echo "  Mapping stats:"
    samtools flagstat "$BAM_OUT" | grep -E "mapped|paired"

    # 3. Compute per-base depth including zero-coverage positions
    echo "[3/4] Calculating per-base depth (zero-coverage positions included)..."
    samtools depth -a "$BAM_OUT" > "$DEPTH_OUT"

    # 4. Compute mean depth per contig, sorted by contig name
    echo "[4/4] Calculating mean depth per contig..."
    awk '
        {
            sum[$1]   += $3
            count[$1] += 1
        }
        END {
            for (contig in sum) {
                mean = sum[contig] / count[contig]
                print contig "\t" mean
            }
        }
    ' "$DEPTH_OUT" | sort -k1,1 > "$MEAN_OUT"

    echo "Done: $sample"
    echo "  BAM:        $BAM_OUT"
    echo "  Bowtie2 log:${BAM_DIR}/${sample}.bowtie2.log"
    echo "  DEPTH:      $DEPTH_OUT"
    echo "  MEAN DEPTH: $MEAN_OUT"
done

echo
echo "Pipeline finished successfully."
echo "Mean depth files are in: ${COV_DIR}/"
