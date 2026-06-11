#!/bin/bash
set -euo pipefail

BAM_DIR="/home/lshpk18/OzanGundogdu_SOUK011275/Analysis/Alignment/bam/host"
KRAKEN_DB="/home/lshpk18/unmapped reads/kraken2_db"
OUT_DIR="/home/lshpk18/unmapped reads/kraken2_output"
TMP_DIR="/home/lshpk18/unmapped reads/tmp_fastq"
THREADS=8

mkdir -p "$OUT_DIR" "$TMP_DIR"

for BAM in "$BAM_DIR"/*.bam; do
    SAMPLE=$(basename "$BAM" .host.Aligned.sortedByCoord.out.bam)
    REPORT="$OUT_DIR/${SAMPLE}.kraken2.report"

    if [ -f "$REPORT" ]; then
        echo "[SKIP] $SAMPLE — report already exists"
        continue
    fi

    echo "[$(date +%H:%M:%S)] Processing: $SAMPLE"

    R1="$TMP_DIR/${SAMPLE}_R1.fastq.gz"
    R2="$TMP_DIR/${SAMPLE}_R2.fastq.gz"

    # Extract both-mates-unmapped reads and convert to FASTQ
    samtools view -f 12 -b -@ "$THREADS" "$BAM" \
    | samtools sort -n -@ "$THREADS" \
    | samtools fastq -@ "$THREADS" -1 "$R1" -2 "$R2" -0 /dev/null -s /dev/null

    # Run Kraken2
    kraken2 \
        --db "$KRAKEN_DB" \
        --threads "$THREADS" \
        --paired \
        --gzip-compressed \
        --output /dev/null \
        --report "$REPORT" \
        "$R1" "$R2"

    # Clean up FASTQ to save space
    rm -f "$R1" "$R2"
    echo "[$(date +%H:%M:%S)] Done: $SAMPLE"
done

rmdir "$TMP_DIR" 2>/dev/null || true

# ── Build summary matrix ──────────────────────────────────────────────────────
echo ""
echo "Building summary matrix..."
MATRIX="$OUT_DIR/kraken2_summary_matrix.txt"

# Header
printf "Sample\tUnclassified_pct\tBacteria_pct\tViruses_pct\tEukaryota_pct\tHomo_sapiens_pct\tCampylobacter_pct\tTop_Bacteria\n" > "$MATRIX"

for REPORT in "$OUT_DIR"/*.kraken2.report; do
    SAMPLE=$(basename "$REPORT" .kraken2.report)

    unclass=$(awk '$4=="U" && $5==0 {printf "%.2f", $1}' "$REPORT")
    bacteria=$(awk '$4=="D" && $6==2 {printf "%.2f", $1}' "$REPORT")
    viruses=$(awk '$4=="D" && $6==10239 {printf "%.2f", $1}' "$REPORT")
    eukaryota=$(awk '$4=="D" && $6==2759 {printf "%.2f", $1}' "$REPORT")
    human=$(awk '$4=="S" && $6==9606 {printf "%.2f", $1}' "$REPORT")
    campy=$(awk '$4=="G" && $6==194 {printf "%.2f", $1}' "$REPORT")

    # Top bacterial species (highest % at species level within Bacteria)
    top_bacteria=$(awk '$4=="S" {print $1"\t"$7}' "$REPORT" | sort -k1,1rn | head -3 | \
        awk '{printf "%s(%s%%) ", $2, $1}' | sed 's/ $//')

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$SAMPLE" "${unclass:-0}" "${bacteria:-0}" "${viruses:-0}" \
        "${eukaryota:-0}" "${human:-0}" "${campy:-0}" "${top_bacteria:-NA}" >> "$MATRIX"
done

echo "Matrix written to: $MATRIX"
echo "=== All done ==="
