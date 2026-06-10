#!/usr/bin/env python3
"""
Pipeline: BLAST bacterial reads (Kraken2-classified) against C. jejuni reference.

Steps:
  1. Parse all Kraken2 reports to collect bacterial taxon IDs and species names.
  2. Extract read IDs from kraken2.output that were classified as bacteria.
  3. Extract those sequences from the unmapped FASTQ files.
  4. Build a BLAST database from the C. jejuni reference genome.
  5. Run blastn.
  6. Generate a per-species summary report.
"""

import os
import gzip
import subprocess
import sys
from collections import defaultdict

# ── Paths ───────────────────────────────────────────────────────────────────
BASE = "/home/lshpk18/unmapped reads (kraken)"
REPORT_DIR = os.path.join(BASE, "kraken2_output")
KRAKEN_OUTPUT = os.path.join(REPORT_DIR, "kraken2.output")
FASTQ_R1 = os.path.join(BASE, "fastq", "unmapped_R1.fastq.gz")
FASTQ_R2 = os.path.join(BASE, "fastq", "unmapped_R2.fastq.gz")
CJEJUNI_REF = "/home/lshpk18/OzanGundogdu_SOUK011275/Analysis/Alignment/references/AL111168.1.fasta"

OUT_DIR = BASE  # parent of kraken2_output
# Use /tmp for BLAST db to avoid issues with spaces in path
BLAST_DB_DIR = "/tmp/cjejuni_blastdb"
BLAST_DB = os.path.join(BLAST_DB_DIR, "AL111168")
BACTERIAL_FASTA = os.path.join(OUT_DIR, "bacterial_reads.fasta")
BLAST_OUT = os.path.join(OUT_DIR, "bacterial_vs_cjejuni.blast6")
REPORT_FILE = os.path.join(OUT_DIR, "bacterial_vs_cjejuni_report.txt")
REPORT_TSV = os.path.join(OUT_DIR, "bacterial_vs_cjejuni_report.tsv")

os.makedirs(BLAST_DB_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)


# ── Step 1: Parse Kraken2 reports ───────────────────────────────────────────
print("=" * 70)
print("STEP 1: Parsing Kraken2 reports for bacterial taxa")
print("=" * 70)

taxid_to_name = {}   # taxon_id (int) -> canonical name
taxid_to_rank = {}   # taxon_id (int) -> rank string
bacterial_taxids = set()  # all taxon IDs under Bacteria domain

report_files = sorted(
    f for f in os.listdir(REPORT_DIR) if f.endswith(".kraken2.report")
)
print(f"  Found {len(report_files)} report files.")

for rfile in report_files:
    path = os.path.join(REPORT_DIR, rfile)
    in_bacteria = False
    bacteria_indent = -1

    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            rank = parts[3]
            try:
                taxid = int(parts[4])
            except ValueError:
                continue
            name_field = parts[5]
            indent = len(name_field) - len(name_field.lstrip())
            name = name_field.strip()

            # Store name/rank for every taxon we see
            if taxid not in taxid_to_name:
                taxid_to_name[taxid] = name
                taxid_to_rank[taxid] = rank

            if taxid == 2:  # Bacteria domain
                in_bacteria = True
                bacteria_indent = indent
                bacterial_taxids.add(taxid)
                continue

            if in_bacteria:
                if indent <= bacteria_indent:
                    # Left the Bacteria subtree
                    in_bacteria = False
                else:
                    bacterial_taxids.add(taxid)

print(f"  Collected {len(bacterial_taxids)} bacterial taxon IDs across all reports.")

# Species-level bacterial taxa (rank == S)
bacterial_species = {
    tid: taxid_to_name[tid]
    for tid in bacterial_taxids
    if taxid_to_rank.get(tid, "").startswith("S")
}
print(f"  Of which {len(bacterial_species)} are species-level (rank S).")


# ── Step 2: Extract bacterial read IDs ──────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 2: Extracting bacterial read IDs from kraken2.output")
print("=" * 70)

read_to_taxid = {}   # read_id -> taxid (for all bacterial reads)

with open(KRAKEN_OUTPUT) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        status = parts[0]
        if status != "C":
            continue
        read_id = parts[1]
        try:
            taxid = int(parts[2])
        except ValueError:
            continue
        if taxid in bacterial_taxids:
            read_to_taxid[read_id] = taxid

print(f"  Found {len(read_to_taxid):,} bacterial reads.")

# Summary by taxon ID
taxid_counts = defaultdict(int)
for tid in read_to_taxid.values():
    taxid_counts[tid] += 1

print("\n  Top 20 bacterial taxa by read count:")
print(f"  {'Taxon ID':<12} {'Reads':>8}  {'Rank':<6}  Name")
print("  " + "-" * 60)
for tid, cnt in sorted(taxid_counts.items(), key=lambda x: -x[1])[:20]:
    rank = taxid_to_rank.get(tid, "?")
    name = taxid_to_name.get(tid, "unknown")
    print(f"  {tid:<12} {cnt:>8}  {rank:<6}  {name}")


# ── Step 3: Extract sequences from FASTQ ────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 3: Extracting bacterial sequences from FASTQ")
print("=" * 70)

bacterial_read_ids = set(read_to_taxid.keys())
written = 0

if os.path.exists(BACTERIAL_FASTA) and os.path.getsize(BACTERIAL_FASTA) > 0:
    print(f"  {BACTERIAL_FASTA} already exists — counting existing sequences ...")
    with open(BACTERIAL_FASTA) as fh:
        written = sum(1 for l in fh if l.startswith(">"))
    print(f"  Reusing {written:,} existing sequences.")
else:
    with open(BACTERIAL_FASTA, "w") as out_fh:
        for fastq_path, suffix in [(FASTQ_R1, "/1"), (FASTQ_R2, "/2")]:
            print(f"  Processing {os.path.basename(fastq_path)} ...")
            try:
                fh = gzip.open(fastq_path, "rt")
            except Exception as e:
                print(f"  ERROR opening {fastq_path}: {e}")
                continue

            while True:
                header = fh.readline()
                if not header:
                    break
                seq = fh.readline().rstrip("\n")
                fh.readline()  # +
                fh.readline()  # quality

                if not header.startswith("@"):
                    continue

                read_id = header[1:].split()[0]
                if read_id in bacterial_read_ids:
                    taxid = read_to_taxid[read_id]
                    species = taxid_to_name.get(taxid, f"taxid_{taxid}").replace(" ", "_")
                    out_fh.write(f">{read_id}{suffix} taxid={taxid} species={species}\n{seq}\n")
                    written += 1

            fh.close()

    print(f"  Wrote {written:,} sequences to {BACTERIAL_FASTA}")


# ── Step 4: Build BLAST database ────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 4: Building BLAST database from C. jejuni reference")
print("=" * 70)

cmd = [
    "makeblastdb",
    "-in", CJEJUNI_REF,
    "-dbtype", "nucl",
    "-out", BLAST_DB,
    "-title", "Campylobacter_jejuni_NCTC11168",
]
print(f"  Running: {' '.join(cmd)}")
result = subprocess.run(cmd, capture_output=True, text=True)
if result.returncode != 0:
    print(f"  ERROR: {result.stderr}")
    sys.exit(1)
print(f"  {result.stdout.strip()}")
print("  BLAST database built successfully.")


# ── Step 5: Run BLASTn ──────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 5: Running BLASTn (bacterial reads vs C. jejuni reference)")
print("=" * 70)

blast_fmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"

cmd = [
    "blastn",
    "-query", BACTERIAL_FASTA,
    "-db", BLAST_DB,
    "-out", BLAST_OUT,
    "-outfmt", blast_fmt,
    "-evalue", "1e-5",
    "-perc_identity", "70",
    "-word_size", "11",
    "-num_threads", "4",
    "-max_target_seqs", "1",
]
print(f"  Running blastn with evalue=1e-5, perc_identity>=70%, word_size=11 ...")
print(f"  Query: {written:,} bacterial reads  |  DB: C. jejuni NCTC 11168")
result = subprocess.run(cmd, capture_output=True, text=True)
if result.returncode != 0:
    print(f"  ERROR: {result.stderr}")
    sys.exit(1)
print(f"  BLASTn completed. Results in {BLAST_OUT}")

# Count hits
with open(BLAST_OUT) as fh:
    blast_lines = [l for l in fh if l.strip()]
print(f"  Total BLAST hits (≥70% identity, e-value ≤1e-5): {len(blast_lines):,}")


# ── Step 6: Generate report ──────────────────────────────────────────────────
print("\n" + "=" * 70)
print("STEP 6: Generating report")
print("=" * 70)

# Parse BLAST output
# Fields: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen
hits_by_taxid = defaultdict(list)

with open(BLAST_OUT) as fh:
    for line in fh:
        if not line.strip():
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 13:
            continue
        qseqid = cols[0]
        pident = float(cols[2])
        aln_len = int(cols[3])
        evalue = float(cols[10])
        bitscore = float(cols[11])
        qlen = int(cols[12])
        coverage = aln_len / qlen * 100 if qlen > 0 else 0

        # Extract taxid from qseqid (format: readID taxid=X species=Y)
        taxid = None
        read_base = qseqid.split("/")[0]  # strip /1 or /2
        if read_base in read_to_taxid:
            taxid = read_to_taxid[read_base]

        if taxid is not None:
            hits_by_taxid[taxid].append({
                "pident": pident,
                "aln_len": aln_len,
                "evalue": evalue,
                "bitscore": bitscore,
                "qlen": qlen,
                "coverage": coverage,
            })

# Compile per-taxon statistics
species_stats = []

all_taxids_queried = set(read_to_taxid.values())

for taxid in all_taxids_queried:
    name = taxid_to_name.get(taxid, f"taxid_{taxid}")
    rank = taxid_to_rank.get(taxid, "?")
    total_reads = taxid_counts[taxid]
    hits = hits_by_taxid.get(taxid, [])
    n_hits = len(hits)
    hit_rate = n_hits / total_reads * 100 if total_reads > 0 else 0

    if hits:
        avg_pident = sum(h["pident"] for h in hits) / n_hits
        max_pident = max(h["pident"] for h in hits)
        avg_cov = sum(h["coverage"] for h in hits) / n_hits
        avg_bitscore = sum(h["bitscore"] for h in hits) / n_hits
        min_evalue = min(h["evalue"] for h in hits)
    else:
        avg_pident = max_pident = avg_cov = avg_bitscore = 0.0
        min_evalue = float("inf")

    species_stats.append({
        "taxid": taxid,
        "name": name,
        "rank": rank,
        "total_reads": total_reads,
        "n_hits": n_hits,
        "hit_rate": hit_rate,
        "avg_pident": avg_pident,
        "max_pident": max_pident,
        "avg_coverage": avg_cov,
        "avg_bitscore": avg_bitscore,
        "min_evalue": min_evalue,
    })

# Sort by number of hits descending
species_stats.sort(key=lambda x: (-x["n_hits"], -x["avg_pident"]))

# Write TSV
with open(REPORT_TSV, "w") as fh:
    fh.write(
        "taxon_id\tspecies_name\trank\ttotal_reads\thits_vs_cjejuni\thit_rate_%\t"
        "avg_pident_%\tmax_pident_%\tavg_query_coverage_%\tavg_bitscore\tmin_evalue\n"
    )
    for s in species_stats:
        ev = f"{s['min_evalue']:.2e}" if s["min_evalue"] < float("inf") else "N/A"
        fh.write(
            f"{s['taxid']}\t{s['name']}\t{s['rank']}\t{s['total_reads']}\t"
            f"{s['n_hits']}\t{s['hit_rate']:.1f}\t{s['avg_pident']:.1f}\t"
            f"{s['max_pident']:.1f}\t{s['avg_coverage']:.1f}\t"
            f"{s['avg_bitscore']:.1f}\t{ev}\n"
        )

# Write human-readable report
total_hits_all = sum(s["n_hits"] for s in species_stats)
taxa_with_hits = [s for s in species_stats if s["n_hits"] > 0]
taxa_no_hits = [s for s in species_stats if s["n_hits"] == 0]

with open(REPORT_FILE, "w") as fh:
    fh.write("=" * 78 + "\n")
    fh.write("BLAST HOMOLOGY REPORT: Bacterial Reads vs Campylobacter jejuni NCTC 11168\n")
    fh.write("=" * 78 + "\n\n")
    fh.write(f"Reference genome : AL111168.1 (C. jejuni NCTC 11168 complete genome)\n")
    fh.write(f"BLAST parameters : evalue ≤ 1e-5, identity ≥ 70%, word_size = 11\n")
    fh.write(f"Kraken2 reports  : {len(report_files)} samples\n")
    fh.write(f"FASTQ files      : unmapped_R1.fastq.gz + unmapped_R2.fastq.gz\n\n")
    fh.write(f"Bacterial taxa identified across all Kraken2 reports : {len(bacterial_taxids):,}\n")
    fh.write(f"Bacterial species (S-rank) identified                : {len(bacterial_species):,}\n")
    fh.write(f"Total bacterial reads queried                        : {len(bacterial_read_ids):,} "
             f"(from {written:,} FASTA entries)\n")
    fh.write(f"Total BLAST hits to C. jejuni                        : {total_hits_all:,}\n")
    fh.write(f"Taxa with ≥1 BLAST hit                               : {len(taxa_with_hits):,}\n")
    fh.write(f"Taxa with 0 BLAST hits                               : {len(taxa_no_hits):,}\n\n")

    fh.write("─" * 78 + "\n")
    fh.write("TAXA WITH BLAST HITS TO C. jejuni (sorted by number of hits)\n")
    fh.write("─" * 78 + "\n\n")
    fh.write(
        f"{'Taxon Name':<45} {'Rank':<5} {'Reads':>7} {'Hits':>6} {'Hit%':>6} "
        f"{'AvgID%':>7} {'MaxID%':>7} {'AvgCov%':>8} {'min-E':>9}\n"
    )
    fh.write("-" * 108 + "\n")

    for s in taxa_with_hits:
        ev = f"{s['min_evalue']:.1e}" if s["min_evalue"] < float("inf") else "N/A"
        fh.write(
            f"{s['name'][:45]:<45} {s['rank']:<5} {s['total_reads']:>7} "
            f"{s['n_hits']:>6} {s['hit_rate']:>5.1f}% "
            f"{s['avg_pident']:>6.1f}% {s['max_pident']:>6.1f}% "
            f"{s['avg_coverage']:>7.1f}% {ev:>9}\n"
        )

    if taxa_no_hits:
        fh.write("\n" + "─" * 78 + "\n")
        fh.write("TAXA WITH NO BLAST HITS TO C. jejuni\n")
        fh.write("─" * 78 + "\n\n")
        fh.write(f"{'Taxon Name':<50} {'Rank':<5} {'Reads':>7}\n")
        fh.write("-" * 65 + "\n")
        for s in sorted(taxa_no_hits, key=lambda x: -x["total_reads"]):
            fh.write(f"{s['name'][:50]:<50} {s['rank']:<5} {s['total_reads']:>7}\n")

    fh.write("\n" + "─" * 78 + "\n")
    fh.write("NOTES\n")
    fh.write("─" * 78 + "\n")
    fh.write(
        "* Reads classified at genus, family, or higher levels are included.\n"
        "* 'Hit%' = percentage of reads for that taxon with ≥1 BLAST hit to C. jejuni.\n"
        "* 'AvgID%' = mean percent identity across all hits for that taxon.\n"
        "* 'AvgCov%' = mean query coverage (alignment length / read length × 100).\n"
        "* High hit rates with low identity (<85%) may indicate conserved regions\n"
        "  shared across Epsilonproteobacteria or other related organisms.\n"
        "* Results are for the single sample in fastq/unmapped_R1/R2.fastq.gz.\n"
        "  The 100 Kraken2 reports provided the species catalogue used as the\n"
        "  taxon reference list.\n"
    )

print(f"\n  Written reports:")
print(f"    {REPORT_FILE}")
print(f"    {REPORT_TSV}")

print("\n" + "=" * 70)
print("PIPELINE COMPLETE")
print("=" * 70)
