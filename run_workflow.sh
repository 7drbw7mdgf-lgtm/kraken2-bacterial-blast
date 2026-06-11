#!/bin/bash
set -euo pipefail

# ── Paths ────────────────────────────────────────────────────────────────────
WORKDIR="/home/lshpk18/unmapped reads"
BAM="$WORKDIR/24h_2_Run_2_S38_L001.host.unmapped.bam"
KRAKEN_DB="$WORKDIR/kraken2_db"
FASTQ_DIR="$WORKDIR/fastq"
KRAKEN_OUT="$WORKDIR/kraken2_output"
ASSEMBLY_OUT="$WORKDIR/assembly"
THREADS=8

mkdir -p "$KRAKEN_DB" "$FASTQ_DIR" "$KRAKEN_OUT" "$ASSEMBLY_OUT"

# ── Step 1: Download Kraken2 standard 16 GB database ─────────────────────────
echo "[1/4] Downloading Kraken2 standard-16 database..."
if [ ! -f "$KRAKEN_DB/hash.k2d" ]; then
    wget -q --show-progress \
        https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20240605.tar.gz \
        -O "$KRAKEN_DB/k2_db.tar.gz"
    echo "Extracting database..."
    tar -xzf "$KRAKEN_DB/k2_db.tar.gz" -C "$KRAKEN_DB"
    rm "$KRAKEN_DB/k2_db.tar.gz"
    echo "Database ready."
else
    echo "Database already present, skipping download."
fi

# ── Step 2: BAM → paired FASTQ ───────────────────────────────────────────────
echo "[2/4] Converting unmapped BAM to paired FASTQ..."
R1="$FASTQ_DIR/unmapped_R1.fastq.gz"
R2="$FASTQ_DIR/unmapped_R2.fastq.gz"

if [ ! -f "$R1" ]; then
    # Sort by name first so samtools fastq pairs reads correctly
    samtools sort -n -@ "$THREADS" "$BAM" \
    | samtools fastq -@ "$THREADS" \
        -1 "$R1" \
        -2 "$R2" \
        -0 /dev/null -s /dev/null
    echo "FASTQ written: $R1, $R2"
else
    echo "FASTQ files already exist, skipping conversion."
fi

# ── Step 3: Kraken2 classification ───────────────────────────────────────────
echo "[3/4] Running Kraken2 classification..."
kraken2 \
    --db "$KRAKEN_DB" \
    --threads "$THREADS" \
    --paired \
    --gzip-compressed \
    --output "$KRAKEN_OUT/kraken2.output" \
    --report "$KRAKEN_OUT/kraken2.report" \
    --report-minimizer-data \
    "$R1" "$R2"
echo "Kraken2 report: $KRAKEN_OUT/kraken2.report"

# ── Step 4: Unicycler de novo assembly ───────────────────────────────────────
echo "[4/4] Running Unicycler assembly on all unmapped reads..."
conda run -n unicycler unicycler \
    -1 "$R1" \
    -2 "$R2" \
    -o "$ASSEMBLY_OUT" \
    -t "$THREADS"
echo "Assembly complete: $ASSEMBLY_OUT/assembly.fasta"

echo ""
echo "=== Workflow complete ==="
echo "  Kraken2 report : $KRAKEN_OUT/kraken2.report"
echo "  Assembly        : $ASSEMBLY_OUT/assembly.fasta"
