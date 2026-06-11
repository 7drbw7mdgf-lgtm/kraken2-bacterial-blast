# Unmapped Reads — Kraken2 Metagenomics Workflow

This document records every step taken on the host-unmapped reads from the dual RNA-seq time-course experiment (OzanGundogdu_SOUK011275). The goal is to characterise the microbial content of reads that did not align to the human genome.

---

## Dataset

- **Project:** OzanGundogdu_SOUK011275 — dual RNA-seq Campylobacter jejuni / human intestinal cell time-course
- **Input BAMs:** 96 host-aligned BAMs in `Analysis/Alignment/bam/host/`
- **Timepoints:** uninfected_0h, 15m, 30m, 1h, 3h, 24h (4 biological × 4 technical replicates = 16 files per timepoint)
- **Sample naming:** `{timepoint}_{bio_rep}_{tech_rep}` (e.g. `15m_1_1`, `uninfected_0h_4_4`)

---

## Phase 1 — Tool Installation

**Script:** `install_tools.sh`

Installs Kraken2 into the base conda environment and creates a dedicated `unicycler` conda environment (Python 3.12) for assembly.

```bash
bash install_tools.sh
```

---

## Phase 2 — Pilot: Single Sample Run

**Script:** `run_workflow.sh`

Prototype run on a single BAM (`24h_2_Run_2_S38_L001.host.unmapped.bam`) to validate the full pipeline:

1. Download Kraken2 standard-16 GB database (k2_standard_16gb_20240605)
2. Convert unmapped BAM → paired FASTQ (`samtools sort -n | samtools fastq`)
3. Run Kraken2 classification
4. Run Unicycler de novo assembly on all unmapped reads

**Outputs:**
- `fastq/unmapped_R1.fastq.gz`, `fastq/unmapped_R2.fastq.gz`
- `kraken2_output/kraken2.output`, `kraken2_output/kraken2.report`
- `assembly/` (Unicycler/SPAdes output)

**Assembly details (unicycler.log):**
- Unicycler v0.5.1, SPAdes 4.2.0
- k-mer range: 25, 49, 69, 83, 97, 107, 115, 123
- Started: 2026-05-19 13:52

---

## Phase 3 — Batch Kraken2 Across All 96 Samples

**Scripts:** `batch_kraken2.sh`, `submit_kraken2.sh`

Runs Kraken2 on every host BAM. For each sample:
1. Extracts both-mates-unmapped reads (`samtools view -f 12`)
2. Name-sorts and converts to paired FASTQ (temp files, deleted after)
3. Runs Kraken2 (`--report-minimizer-data` mode, `--paired`)
4. Deletes temp FASTQs to save disk

After all samples, builds a summary matrix (`kraken2_output/kraken2_summary_matrix.txt`) with per-sample columns: unclassified %, Bacteria %, Viruses %, Eukaryota %, Homo sapiens %, Campylobacter %, top bacterial species.

**SLURM submission:**
```bash
sbatch submit_kraken2.sh   # 8 CPUs, 32 GB, 12 h
```

**Outputs:**
- `kraken2_output/{sample}.kraken2.report` — 96 files (one per sample)
- `kraken2_output/kraken2_summary_matrix.txt`

---

## Phase 4 — Species Matrix

**Script:** `build_species_matrix.py`

Parses all 96 Kraken2 reports and builds a species × sample matrix of bacterial species percentages. Used to identify which species appear consistently across timepoints.

**Output:** `kraken2_output/bacterial_species_matrix.txt`

---

## Phase 5 — BLASTn Homology: Bacterial Reads vs *C. jejuni*

**Script:** `blast_vs_cjejuni_pipeline.py`

End-to-end pipeline to test whether bacteria identified by Kraken2 share sequence homology with *C. jejuni* NCTC 11168:

1. Parse all 96 Kraken2 reports → collect bacterial taxon IDs and species names
2. Extract read IDs classified as bacteria from `kraken2_output/kraken2.output`
3. Extract those sequences from the pilot-sample FASTQ (`fastq/unmapped_R1/R2.fastq.gz`)
4. Build a BLAST nucleotide database from the *C. jejuni* reference genome (AL111168.1)
5. Run BLASTn (evalue ≤ 1e-5, identity ≥ 70%, word_size 11)
6. Generate per-species report

**Outputs:**
- `bacterial_reads.fasta` — bacterial reads extracted for querying
- `bacterial_vs_cjejuni.blast6` — tabular BLAST output (fmt 6)
- `bacterial_vs_cjejuni_report.txt` — human-readable report with hit rates, avg identity, coverage per taxon
- `bacterial_vs_cjejuni_report.tsv` — machine-readable equivalent

---

## Phase 6 — Interactive HTML Homology Report

**Script:** `generate_homology_html.py`

Generates an interactive HTML table (species × sample) combining:
- **Cell colour (blue gradient):** read count abundance per species per sample (log scale)
- **Cell value:** % of reads for that species with a BLAST hit to *C. jejuni*

Features: sortable columns, search box, min-reads and min-samples filters. File size ~865 KB.

**Output:** `kraken2_cjejuni_homology.html`

**Committed versions:**
- Initial table (2026-06-10, commit `ddeba77`): 5,302 taxa shown
- Rebuild with BLAST % colouring (commit `65056d6`)
- Hybrid dual-encoding: value = BLAST %, colour = read count (commit `8c1df6f`)

---

## Phase 7 — Sample Naming Fixes

Two naming corrections applied to Kraken2 report filenames and `bacterial_species_matrix.txt`:

1. **Commit `a2a5e32`:** Remove `H_` prefix from non-0h timepoints (15m, 30m, 1h, 3h, 24h). Host vs bacteria labelling was not needed for post-alignment unmapped reads.
2. **Commit `6dfdf84`:** Rename `B_0h` → `uninfected_0h`. The 0h timepoint represents uninfected controls, so naming by infection status is more informative.

Final naming convention: `{timepoint}_{bio_rep}_{tech_rep}` (e.g. `uninfected_0h_1_1`, `24h_2_3`).

---

## Phase 8 — Extract Unmapped Reads by Timepoint

**Scripts:** `extract_by_timepoint.sh`, `submit_extract_timepoints.sh`

Extracts both-mates-unmapped reads from all 96 host BAMs and writes them as FASTA files organised by timepoint. Intended for downstream per-timepoint metagenomics or assembly.

**Method:** `samtools view -f 12` (both mates unmapped) → name-sort → `samtools fasta` with `/1`/`/2` mate suffixes

**SLURM submission:**
```bash
sbatch submit_extract_timepoints.sh   # 8 CPUs, 16 GB, 8 h
```

**Job 6141636** ran on 2026-06-10 16:24–16:51 (27 min). All 96 samples completed without errors.

**Read counts per timepoint (total sequences written):**

| Timepoint | Files | Notes |
|-----------|-------|-------|
| 15m | 16 | 585 K – 6.6 M reads per file |
| 30m | 16 | 443 K – 5.4 M reads per file |
| 1h | 16 | 409 K – 4.7 M reads per file |
| 3h | 16 | 560 K – 3.5 M reads per file |
| 24h | 16 | 558 K – 3.5 M reads per file |
| uninfected_0h | 16 | 585 K – 9.5 M reads per file |

**Output:** `sequences_by_timepoint/{timepoint}/{sample}.unmapped.fasta` — 96 files

---

## Phase 9 — Cutibacterium acnes Targeted Extraction

**Scripts:** `extract_cutibacterium_reads.sh`, `submit_extract_cutibacterium.sh`

*Cutibacterium acnes* (taxid 1743) emerged as a consistently present species in the Kraken2 reports. This script re-runs Kraken2 per sample with per-read output enabled, extracts only reads assigned to taxid 1743, and concatenates them into a single FASTA for BLAST verification.

```bash
sbatch submit_extract_cutibacterium.sh
```

**Output:** `cutibacterium_reads/cutibacterium_acnes_all_samples.fasta`

---

## Phase 10 — Extract All Classified Reads for Inspection

**Script:** `extract_kraken2_sequences.py`

Extracts sequences for all reads present in `kraken2.output` (both classified and unclassified) from the pilot FASTQ files, writing status and taxid into FASTA headers. Useful for inspecting the full classification landscape of the single-sample run.

**Output:** `all_kraken2_reads.fasta`

---

## Script Index

| Script | Purpose |
|--------|---------|
| `install_tools.sh` | Install Kraken2 and Unicycler conda environments |
| `run_workflow.sh` | Pilot: single BAM → FASTQ → Kraken2 → Unicycler assembly |
| `batch_kraken2.sh` | Batch Kraken2 across all 96 host BAMs |
| `submit_kraken2.sh` | SLURM wrapper for `batch_kraken2.sh` |
| `build_species_matrix.py` | Build species × sample matrix from Kraken2 reports |
| `blast_vs_cjejuni_pipeline.py` | Extract bacterial reads and BLAST vs *C. jejuni* NCTC 11168 |
| `generate_homology_html.py` | Interactive HTML: species × sample with BLAST homology |
| `extract_by_timepoint.sh` | Extract unmapped reads per sample, organised by timepoint |
| `submit_extract_timepoints.sh` | SLURM wrapper for `extract_by_timepoint.sh` |
| `extract_cutibacterium_reads.sh` | Re-run Kraken2 per sample and extract *C. acnes* reads |
| `submit_extract_cutibacterium.sh` | SLURM wrapper for `extract_cutibacterium_reads.sh` |
| `extract_kraken2_sequences.py` | Extract all reads from kraken2.output for inspection |

---

## Key Outputs

| Output | Path |
|--------|------|
| Per-sample Kraken2 reports (96) | `kraken2_output/{sample}.kraken2.report` |
| Summary matrix (all samples) | `kraken2_output/kraken2_summary_matrix.txt` |
| Species × sample matrix | `kraken2_output/bacterial_species_matrix.txt` |
| BLAST results vs *C. jejuni* | `bacterial_vs_cjejuni.blast6` |
| BLAST report (text) | `bacterial_vs_cjejuni_report.txt` |
| BLAST report (TSV) | `bacterial_vs_cjejuni_report.tsv` |
| Interactive HTML report | `kraken2_cjejuni_homology.html` |
| Unmapped reads by timepoint | `sequences_by_timepoint/{timepoint}/{sample}.unmapped.fasta` |
| *C. acnes* combined FASTA | `cutibacterium_reads/cutibacterium_acnes_all_samples.fasta` |

---

## Environment

- **Kraken2 database:** k2_standard_16gb_20240605 (16 GB standard database)
- **Conda environments:** `base` (Kraken2), `unicycler` (Unicycler v0.5.1, SPAdes 4.2.0), `dualrnaseq_deseq2` (samtools)
- **Cluster:** SLURM on LSHTM HPC
- **BLAST version:** 2.17.0+
- **Reference:** *C. jejuni* NCTC 11168 complete genome (AL111168.1)
