#!/usr/bin/env python3
"""Extract sequences for all reads in kraken2.output from paired fastq files."""

import gzip
import sys
from pathlib import Path

base = Path("/home/lshpk18/unmapped reads (kraken)")
kraken_output = base / "kraken2_output/kraken2.output"
r1_fastq = base / "fastq/unmapped_R1.fastq.gz"
r2_fastq = base / "fastq/unmapped_R2.fastq.gz"
out_fasta = base / "all_kraken2_reads.fasta"

print("Loading Kraken2 classifications...", flush=True)
classifications = {}
with open(kraken_output) as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        status = parts[0]   # C or U
        read_id = parts[1]
        taxid = parts[2]
        classifications[read_id] = (status, taxid)

print(f"  Loaded {len(classifications):,} classifications", flush=True)

def fastq_gz_records(path):
    with gzip.open(path, "rt") as f:
        while True:
            header = f.readline().rstrip("\n")
            if not header:
                break
            seq = f.readline().rstrip("\n")
            f.readline()  # +
            f.readline()  # quality
            read_id = header[1:].split()[0]  # strip @ and any trailing description
            yield read_id, seq

print("Writing FASTA output...", flush=True)
written = 0
missing = 0

with open(out_fasta, "w") as out:
    for mate, fastq_path in [("1", r1_fastq), ("2", r2_fastq)]:
        print(f"  Processing R{mate}...", flush=True)
        for read_id, seq in fastq_gz_records(fastq_path):
            if read_id in classifications:
                status, taxid = classifications[read_id]
                out.write(f">{read_id}/{mate} status={status} taxid={taxid}\n{seq}\n")
                written += 1
            else:
                missing += 1

print(f"\nDone. Written: {written:,} sequences | Missing classifications: {missing:,}")
print(f"Output: {out_fasta}")
