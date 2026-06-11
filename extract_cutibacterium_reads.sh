#!/bin/bash
# Re-runs Kraken2 per sample keeping the per-read output, extracts reads
# classified as Cutibacterium acnes (taxid 1743), and writes a combined FASTA
# for BLAST.
set -euo pipefail

BAM_DIR="/home/lshpk18/OzanGundogdu_SOUK011275/Analysis/Alignment/bam/host"
KRAKEN_DB="/home/lshpk18/unmapped reads/kraken2_db"
OUT_DIR="/home/lshpk18/unmapped reads (kraken)/cutibacterium_reads"
TMP_DIR="/home/lshpk18/unmapped reads (kraken)/tmp_fastq_cutib"
TAXID=1743          # Cutibacterium acnes
THREADS=8

mkdir -p "$OUT_DIR" "$TMP_DIR"

COMBINED_FASTA="$OUT_DIR/cutibacterium_acnes_all_samples.fasta"
> "$COMBINED_FASTA"   # empty/create

TOTAL_READS=0

for BAM in "$BAM_DIR"/*.bam; do
    [[ "$BAM" == *.bai ]] && continue
    SAMPLE=$(basename "$BAM" .host.Aligned.sortedByCoord.out.bam)

    echo "[$(date +%H:%M:%S)] $SAMPLE"

    R1="$TMP_DIR/${SAMPLE}_R1.fastq.gz"
    R2="$TMP_DIR/${SAMPLE}_R2.fastq.gz"
    K2_OUTPUT="$TMP_DIR/${SAMPLE}.kraken2.output"

    # Extract both-mates-unmapped reads → paired FASTQ
    samtools view -f 12 -b -@ "$THREADS" "$BAM" \
    | samtools sort -n -@ "$THREADS" \
    | samtools fastq -@ "$THREADS" -1 "$R1" -2 "$R2" -0 /dev/null -s /dev/null

    # Run Kraken2 — keep the per-read output this time
    kraken2 \
        --db "$KRAKEN_DB" \
        --threads "$THREADS" \
        --paired \
        --gzip-compressed \
        --output "$K2_OUTPUT" \
        --report /dev/null \
        "$R1" "$R2"

    # Extract read IDs assigned to taxid 1743
    READIDS="$TMP_DIR/${SAMPLE}_cutib_readids.txt"
    awk '$1=="C" && $3=='"$TAXID"' {print $2}' "$K2_OUTPUT" > "$READIDS"
    N=$(wc -l < "$READIDS")

    if [ "$N" -gt 0 ]; then
        echo "  → $N read-pairs assigned to taxid $TAXID"
        TOTAL_READS=$((TOTAL_READS + N))

        # Extract those reads from R1 and R2 → append to combined FASTA
        python3 - "$READIDS" "$R1" "$R2" "$SAMPLE" >> "$COMBINED_FASTA" <<'PYEOF'
import sys, gzip

readid_file, r1_path, r2_path, sample = sys.argv[1:]

with open(readid_file) as fh:
    wanted = set(line.strip() for line in fh)

def iter_fastq(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        while True:
            header = fh.readline().strip()
            if not header:
                break
            seq  = fh.readline().strip()
            fh.readline()   # +
            fh.readline()   # qual
            # Read name is everything after '@' up to the first space
            name = header[1:].split()[0]
            yield name, seq

for path, mate in [(r1_path, "R1"), (r2_path, "R2")]:
    for name, seq in iter_fastq(path):
        if name in wanted:
            print(f">{sample}_{mate}_{name}")
            print(seq)
PYEOF
    else
        echo "  → 0 reads for taxid $TAXID, skipping"
    fi

    # Clean up temp files for this sample
    rm -f "$R1" "$R2" "$K2_OUTPUT" "$READIDS"
done

rmdir "$TMP_DIR" 2>/dev/null || true

echo ""
echo "=== Done ==="
echo "Total read-pairs extracted : $TOTAL_READS"
echo "Combined FASTA             : $COMBINED_FASTA"
echo ""
echo "To BLAST online: https://blast.ncbi.nlm.nih.gov/Blast.cgi"
echo "  Program: blastn, Database: nt, check 'Somewhat similar sequences'"
echo "  Or nucleotide BLAST against nr/nt"
