#!/bin/bash
set -euo pipefail

BAM_DIR="/home/lshpk18/OzanGundogdu_SOUK011275/Analysis/Alignment/bam/host"
OUT_DIR="/home/lshpk18/unmapped reads (kraken)/sequences_by_timepoint"
SAMTOOLS="/home/lshpk18/miniconda3/envs/dualrnaseq_deseq2/bin/samtools"
THREADS=8

# Time point prefixes
TIMEPOINTS=("15m" "30m" "1h" "3h" "24h" "uninfected_0h")

for TP in "${TIMEPOINTS[@]}"; do
    mkdir -p "$OUT_DIR/$TP"
done

echo "[$(date)] Starting extraction of unmapped reads from all BAMs"
echo "Output directory: $OUT_DIR"
echo ""

TOTAL=$(ls "$BAM_DIR"/*.bam | wc -l)
COUNT=0

for BAM in "$BAM_DIR"/*.bam; do
    SAMPLE=$(basename "$BAM" .host.Aligned.sortedByCoord.out.bam)
    COUNT=$((COUNT + 1))

    # Determine time point from sample name
    TIMEPOINT=""
    for TP in "${TIMEPOINTS[@]}"; do
        if [[ "$SAMPLE" == ${TP}_* ]]; then
            TIMEPOINT="$TP"
            break
        fi
    done

    if [ -z "$TIMEPOINT" ]; then
        echo "[WARN] Could not determine time point for: $SAMPLE — skipping"
        continue
    fi

    OUT_FASTA="$OUT_DIR/$TIMEPOINT/${SAMPLE}.unmapped.fasta"

    if [ -f "$OUT_FASTA" ]; then
        echo "[SKIP $COUNT/$TOTAL] $SAMPLE — already exists"
        continue
    fi

    echo "[$(date +%H:%M:%S) $COUNT/$TOTAL] $SAMPLE → $TIMEPOINT/"

    # Extract both-mates-unmapped reads, sort by name, output as FASTA
    # -f 12: both mates unmapped; /1 /2 suffixes added for R1/R2
    $SAMTOOLS view -f 12 -b -@ "$THREADS" "$BAM" \
    | $SAMTOOLS sort -n -@ "$THREADS" \
    | $SAMTOOLS fasta -@ "$THREADS" \
        -1 >(sed 's/^>/>\/1 /' >> "$OUT_FASTA") \
        -2 >(sed 's/^>/>\/2 /' >> "$OUT_FASTA") \
        -0 /dev/null -s /dev/null

    NSEQ=$(grep -c "^>" "$OUT_FASTA")
    echo "  -> $NSEQ sequences written"
done

echo ""
echo "[$(date)] All done. Files written to: $OUT_DIR"
echo ""
echo "Summary:"
for TP in "${TIMEPOINTS[@]}"; do
    N=$(ls "$OUT_DIR/$TP/"*.fasta 2>/dev/null | wc -l)
    echo "  $TP: $N files"
done
