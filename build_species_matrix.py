import os
import glob
from collections import defaultdict

REPORT_DIR = "/home/lshpk18/unmapped reads/kraken2_output"
OUT_FILE = os.path.join(REPORT_DIR, "bacterial_species_matrix.txt")

def parse_report(filepath):
    """Return dict of {species_name: pct} for bacterial species only."""
    species = {}
    in_bacteria = False
    bacteria_depth = None

    with open(filepath) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            pct = float(parts[0].strip())
            rank = parts[3].strip()
            taxid = parts[4].strip()
            name = parts[5].strip()

            # Detect entry into Bacteria domain
            if rank == "D" and taxid == "2":
                in_bacteria = True
                bacteria_depth = len(parts[5]) - len(parts[5].lstrip())
                continue

            if in_bacteria:
                depth = len(parts[5]) - len(parts[5].lstrip())
                # Left of bacteria domain means we've exited it
                if rank == "D" and taxid != "2":
                    in_bacteria = False
                    continue
                # Collect species-level entries
                if rank == "S" and pct > 0:
                    species[name] = pct

    return species

# Parse all reports
report_files = sorted(glob.glob(os.path.join(REPORT_DIR, "*.kraken2.report")))
print(f"Found {len(report_files)} reports")

sample_data = {}
all_species = set()

for rf in report_files:
    sample = os.path.basename(rf).replace(".kraken2.report", "")
    data = parse_report(rf)
    sample_data[sample] = data
    all_species.update(data.keys())

# Sort species by mean abundance across samples (descending)
species_means = {}
for sp in all_species:
    vals = [sample_data[s].get(sp, 0.0) for s in sample_data]
    species_means[sp] = sum(vals) / len(vals)

sorted_species = sorted(all_species, key=lambda s: species_means[s], reverse=True)

print(f"Total bacterial species found: {len(sorted_species)}")

# Write matrix
with open(OUT_FILE, "w") as out:
    header = "Sample\t" + "\t".join(sorted_species)
    out.write(header + "\n")

    for sample in sorted(sample_data.keys()):
        row = [sample]
        for sp in sorted_species:
            row.append(f"{sample_data[sample].get(sp, 0.0):.4f}")
        out.write("\t".join(row) + "\n")

print(f"Matrix written to: {OUT_FILE}")
print(f"Dimensions: {len(sample_data)} samples x {len(sorted_species)} species")
