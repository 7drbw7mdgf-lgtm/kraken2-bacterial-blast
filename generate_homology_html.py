#!/usr/bin/env python3
"""
Generate an interactive HTML report showing, for every taxon found across all
96 Kraken2 reports, its per-sample read counts and its % homology to the
C. jejuni NCTC 11168 reference genome (AL111168.1) derived from BLASTn.
"""

import os, json, math
from collections import defaultdict

# ── Paths ────────────────────────────────────────────────────────────────────
REPORT_DIR   = "/home/lshpk18/unmapped reads (kraken)/kraken2_output"
KRAKEN_OUT   = os.path.join(REPORT_DIR, "kraken2.output")
BLAST6       = "/home/lshpk18/unmapped reads (kraken)/bacterial_vs_cjejuni.blast6"
HTML_OUT     = "/home/lshpk18/unmapped reads (kraken)/kraken2_cjejuni_homology.html"
REF          = "AL111168.1 — C. jejuni NCTC 11168 (complete genome)"

CONDITION_ORDER = ["uninfected_0h", "15m", "30m", "1h", "3h", "24h"]

# ── Step 1: Parse Kraken2 reports ────────────────────────────────────────────
print("Parsing Kraken2 reports ...")

report_files = sorted(f for f in os.listdir(REPORT_DIR) if f.endswith(".kraken2.report"))

# sample name = filename without extension
samples = [f.replace(".kraken2.report", "") for f in report_files]

# taxid → {name, rank}
taxid_meta = {}
# sample → taxid → direct_reads
sample_taxid_reads = {}

for rfile, sname in zip(report_files, samples):
    path = os.path.join(REPORT_DIR, rfile)
    counts = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            try:
                direct = int(parts[2])
                taxid  = int(parts[4])
            except ValueError:
                continue
            rank = parts[3]
            name = parts[5].strip()
            if taxid not in taxid_meta:
                taxid_meta[taxid] = {"name": name, "rank": rank}
            if direct > 0:
                counts[taxid] = counts.get(taxid, 0) + direct
    sample_taxid_reads[sname] = counts

print(f"  {len(samples)} samples, {len(taxid_meta)} unique taxa.")

# ── Step 2: Build per-taxon stats across samples ─────────────────────────────
# For each taxon: total reads across all samples, # samples present
taxid_total  = defaultdict(int)
taxid_nsamps = defaultdict(int)

for sname, counts in sample_taxid_reads.items():
    for tid, cnt in counts.items():
        taxid_total[tid]  += cnt
        taxid_nsamps[tid] += 1

# ── Step 3: Load kraken2.output → read_id → taxid ───────────────────────────
print("Loading kraken2.output read-taxon map ...")
read_to_taxid = {}
with open(KRAKEN_OUT) as fh:
    for line in fh:
        parts = line.split("\t")
        if len(parts) < 3 or parts[0] != "C":
            continue
        try:
            read_to_taxid[parts[1]] = int(parts[2])
        except ValueError:
            pass

# reads per taxid (from the single sample)
taxid_read_count = defaultdict(int)
for tid in read_to_taxid.values():
    taxid_read_count[tid] += 1

print(f"  {len(read_to_taxid):,} classified reads loaded.")

# ── Step 4: Parse BLAST results → per-taxid homology ────────────────────────
print("Computing per-taxon BLAST homology ...")

# For each read: best hit pident, bitscore
# blast6 fields: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen
read_best_hit = {}  # read_id (without /1 /2) → max pident

with open(BLAST6) as fh:
    for line in fh:
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 12:
            continue
        qid   = cols[0].rsplit("/", 1)[0]   # strip /1 or /2
        pident = float(cols[2])
        if qid not in read_best_hit or pident > read_best_hit[qid]:
            read_best_hit[qid] = pident

# Aggregate per taxid
taxid_hits      = defaultdict(int)    # reads with any blast hit
taxid_pident_sum = defaultdict(float) # sum of best pidents for hit reads
taxid_pident_max = defaultdict(float) # max pident

for rid, pident in read_best_hit.items():
    tid = read_to_taxid.get(rid)
    if tid is None:
        continue
    taxid_hits[tid]       += 1
    taxid_pident_sum[tid] += pident
    if pident > taxid_pident_max[tid]:
        taxid_pident_max[tid] = pident

# Per-taxid homology stats
def taxid_homology(tid):
    n_reads = taxid_read_count.get(tid, 0)
    n_hits  = taxid_hits.get(tid, 0)
    if n_reads == 0:
        return None   # not in the BLAST sample
    hit_pct   = n_hits / n_reads * 100
    avg_pident = taxid_pident_sum[tid] / n_hits if n_hits > 0 else 0
    max_pident = taxid_pident_max.get(tid, 0)
    return {"reads_in_blast_sample": n_reads,
            "blast_hits": n_hits,
            "hit_pct": round(hit_pct, 1),
            "avg_pident": round(avg_pident, 1),
            "max_pident": round(max_pident, 1)}

# ── Step 5: Sort conditions and samples ──────────────────────────────────────
def sample_condition(sname):
    for cond in CONDITION_ORDER:
        if sname.startswith(cond):
            return cond
    return "other"

samples_sorted = []
for cond in CONDITION_ORDER:
    grp = sorted(s for s in samples if sample_condition(s) == cond)
    samples_sorted.extend(grp)
other = sorted(s for s in samples if sample_condition(s) == "other")
samples_sorted.extend(other)

# Condition boundaries for column grouping
cond_spans = []
prev_cond = None
for i, s in enumerate(samples_sorted):
    c = sample_condition(s)
    if c != prev_cond:
        cond_spans.append({"cond": c, "start": i, "count": 1})
        prev_cond = c
    else:
        cond_spans[-1]["count"] += 1

# ── Step 6: Filter taxa for display ─────────────────────────────────────────
# Keep taxa with ≥5 total reads across all samples, or with any BLAST hit
display_taxa = []
for tid, total in sorted(taxid_total.items(), key=lambda x: -x[1]):
    if total < 2:
        continue
    meta = taxid_meta.get(tid, {"name": f"taxid_{tid}", "rank": "?"})
    hom  = taxid_homology(tid)
    per_sample = {s: sample_taxid_reads[s].get(tid, 0) for s in samples_sorted}
    display_taxa.append({
        "taxid":       tid,
        "name":        meta["name"],
        "rank":        meta["rank"],
        "total_reads": total,
        "n_samples":   taxid_nsamps[tid],
        "homology":    hom,
        "per_sample":  per_sample,
    })

print(f"  {len(display_taxa)} taxa for display (≥2 total reads).")

# ── Step 7: Build JavaScript data ────────────────────────────────────────────
max_reads_per_cell = max(
    (v for t in display_taxa for v in t["per_sample"].values()), default=1
)

# Serialise to JSON for embedding
js_samples     = json.dumps(samples_sorted)
js_cond_spans  = json.dumps(cond_spans)
js_taxa = json.dumps([{
    "taxid":       t["taxid"],
    "name":        t["name"],
    "rank":        t["rank"],
    "total":       t["total_reads"],
    "n_samp":      t["n_samples"],
    "hit_pct":     t["homology"]["hit_pct"]    if t["homology"] else None,
    "avg_pid":     t["homology"]["avg_pident"] if t["homology"] else None,
    "max_pid":     t["homology"]["max_pident"] if t["homology"] else None,
    "blast_reads": t["homology"]["reads_in_blast_sample"] if t["homology"] else 0,
    "blast_hits":  t["homology"]["blast_hits"] if t["homology"] else 0,
    "counts":      [t["per_sample"][s] for s in samples_sorted],
} for t in display_taxa])

js_max = max_reads_per_cell

# ── Step 8: Render HTML ───────────────────────────────────────────────────────
print("Rendering HTML ...")

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Kraken2 Taxa vs C. jejuni Homology</title>
<style>
*{{box-sizing:border-box;margin:0;padding:0;font-family:'Segoe UI',Arial,sans-serif}}
body{{background:#f4f6f9;color:#222;font-size:13px}}
header{{background:#1a3a5c;color:#fff;padding:18px 24px}}
header h1{{font-size:20px;font-weight:600}}
header p{{opacity:.8;font-size:12px;margin-top:4px}}
.controls{{background:#fff;padding:14px 24px;border-bottom:1px solid #dde;display:flex;flex-wrap:wrap;gap:16px;align-items:center}}
.controls label{{font-weight:600;margin-right:4px}}
.controls input[type=range]{{width:160px;vertical-align:middle}}
.controls input[type=text]{{border:1px solid #ccc;border-radius:4px;padding:4px 8px;width:200px}}
.controls select{{border:1px solid #ccc;border-radius:4px;padding:4px 8px}}
.stat-val{{font-weight:700;color:#1a3a5c}}
.summary{{display:flex;flex-wrap:wrap;gap:12px;padding:14px 24px;background:#fff;border-bottom:1px solid #dde}}
.card{{background:#f0f4ff;border:1px solid #c5d0ea;border-radius:6px;padding:10px 18px;min-width:140px}}
.card .val{{font-size:22px;font-weight:700;color:#1a3a5c}}
.card .lbl{{font-size:11px;color:#555;margin-top:2px}}
.wrap{{overflow:auto;max-height:calc(100vh - 220px);padding:0 24px 24px}}
table{{border-collapse:collapse;width:max-content;min-width:100%}}
thead th{{position:sticky;top:0;background:#1a3a5c;color:#fff;padding:6px 8px;white-space:nowrap;cursor:pointer;user-select:none;z-index:2;font-size:12px}}
thead th:hover{{background:#2a5a8c}}
thead .cond-hdr{{background:#0e2a45;text-align:center;font-size:11px;letter-spacing:.5px}}
tbody tr:nth-child(even){{background:#f9fafb}}
tbody tr:hover{{background:#e8eef8}}
td{{padding:4px 8px;white-space:nowrap;border-bottom:1px solid #eee;font-size:12px}}
td.num{{text-align:right;font-variant-numeric:tabular-nums}}
td.name{{max-width:260px;overflow:hidden;text-overflow:ellipsis}}
td.cell{{width:36px;text-align:center;font-size:10px;border:1px solid rgba(0,0,0,.04)}}
.rank{{display:inline-block;background:#dde;border-radius:3px;padding:1px 5px;font-size:10px;color:#444}}
.hom-bar{{display:inline-block;height:12px;border-radius:3px;vertical-align:middle;margin-right:4px}}
.sort-arrow{{font-size:10px;opacity:.6}}
.na{{color:#bbb;font-size:11px}}
.legend{{display:flex;align-items:center;gap:8px;font-size:11px}}
.legend-grad{{width:100px;height:12px;border-radius:3px;background:linear-gradient(to right,#f0f4ff,#d32f2f)}}
#rowcount{{font-size:12px;color:#555;padding:8px 24px;background:#f8f8f8;border-bottom:1px solid #eee}}
.filters-active{{color:#c00;font-weight:700}}
</style>
</head>
<body>
<header>
  <h1>Kraken2 Taxa — C. jejuni NCTC 11168 Homology Report</h1>
  <p>Reference: {REF} &nbsp;|&nbsp;
     BLAST: blastn, e-value ≤ 1e-5, identity ≥ 70%, word_size 11 &nbsp;|&nbsp;
     {len(samples)} samples &nbsp;|&nbsp; {len(display_taxa)} taxa shown</p>
</header>

<div class="summary">
  <div class="card"><div class="val">{len(samples)}</div><div class="lbl">Samples</div></div>
  <div class="card"><div class="val">{len(display_taxa)}</div><div class="lbl">Taxa (≥2 reads)</div></div>
  <div class="card"><div class="val">{sum(t["total_reads"] for t in display_taxa):,}</div><div class="lbl">Total classified reads</div></div>
  <div class="card"><div class="val">{sum(1 for t in display_taxa if t["homology"] and t["homology"]["hit_pct"]>0)}</div><div class="lbl">Taxa with BLAST hits</div></div>
  <div class="card"><div class="val">{sum(1 for t in display_taxa if t["homology"] and t["homology"]["hit_pct"]>=50)}</div><div class="lbl">Taxa ≥50% homology</div></div>
</div>

<div class="controls">
  <span><label>Search taxon:</label><input type="text" id="searchBox" placeholder="e.g. Campylobacter" oninput="applyFilters()"></span>
  <span><label>Min total reads:</label><input type="range" id="minReads" min="0" max="1000" step="1" value="2" oninput="document.getElementById('minReadsVal').textContent=this.value;applyFilters()"><span id="minReadsVal" class="stat-val">2</span></span>
  <span><label>Min samples:</label><input type="range" id="minSamps" min="1" max="{len(samples)}" step="1" value="1" oninput="document.getElementById('minSampsVal').textContent=this.value;applyFilters()"><span id="minSampsVal" class="stat-val">1</span></span>
  <span><label>Rank:</label>
    <select id="rankFilter" onchange="applyFilters()">
      <option value="">All</option>
      <option value="S">Species (S)</option>
      <option value="G">Genus (G)</option>
      <option value="F">Family (F)</option>
      <option value="D">Domain (D)</option>
    </select>
  </span>
  <span><label>BLAST tested only:</label><input type="checkbox" id="blastOnly" onchange="applyFilters()"></span>
  <span class="legend"><span>Homology %:</span><span class="legend-grad"></span><span>0 → 100%</span></span>
</div>
<div id="rowcount"></div>

<div class="wrap">
<table id="mainTable">
<thead>
<tr id="condRow">
  <th colspan="7" style="background:#0e2a45"></th>
"""

# Condition header row
for span in cond_spans:
    html += f'  <th colspan="{span["count"]}" class="cond-hdr">{span["cond"].replace("_"," ")}</th>\n'

html += """</tr>
<tr id="headerRow">
  <th onclick="sortBy('name')">Taxon <span class="sort-arrow">▲▼</span></th>
  <th onclick="sortBy('rank')">Rank <span class="sort-arrow">▲▼</span></th>
  <th onclick="sortBy('total')">Total reads <span class="sort-arrow">▲▼</span></th>
  <th onclick="sortBy('n_samp')"># Samples <span class="sort-arrow">▲▼</span></th>
  <th onclick="sortBy('hit_pct')">% Reads hit C. jejuni <span class="sort-arrow">▲▼</span></th>
  <th onclick="sortBy('avg_pid')">Avg identity % <span class="sort-arrow">▲▼</span></th>
  <th onclick="sortBy('max_pid')">Max identity % <span class="sort-arrow">▲▼</span></th>
"""

for s in samples_sorted:
    html += f'  <th title="{s}">{s.replace("uninfected_0h","U0h").replace("_"," ")}</th>\n'

html += """</tr></thead>
<tbody id="tableBody"></tbody>
</table>
</div>

<script>
const SAMPLES = """ + js_samples + """;
const COND_SPANS = """ + js_cond_spans + """;
const TAXA = """ + js_taxa + """;
const MAX_READS = """ + str(js_max) + """;

let sortKey = 'total';
let sortDir = -1;
let filtered = [...TAXA];

function heatColor(val, maxVal) {
  if (val === 0) return '#f9fafb';
  const t = Math.log1p(val) / Math.log1p(maxVal);
  const r = Math.round(255 * t);
  const g = Math.round(255 * (1 - t * 0.6));
  const b = Math.round(255 * (1 - t));
  return `rgb(${r},${g},${b})`;
}

function homColor(pct) {
  if (pct === null) return '#f0f0f0';
  const t = pct / 100;
  const r = Math.round(220 * t + 240 * (1-t));
  const g = Math.round(47 * t + 244 * (1-t));
  const b = Math.round(47 * t + 255 * (1-t));
  return `rgb(${r},${g},${b})`;
}

function renderTable() {
  const tbody = document.getElementById('tableBody');
  const rows = [];
  for (const t of filtered) {
    const hp = t.hit_pct;
    const hpStr = hp !== null ? hp.toFixed(1)+'%' : '<span class="na">—</span>';
    const barW = hp !== null ? Math.max(1, Math.round(hp)) : 0;
    const barStyle = hp !== null ? `background:${homColor(hp)};width:${barW}px;` : '';
    const avgStr = t.avg_pid !== null ? t.avg_pid.toFixed(1)+'%' : '<span class="na">—</span>';
    const maxStr = t.max_pid !== null ? t.max_pid.toFixed(1)+'%' : '<span class="na">—</span>';

    let cells = '';
    for (let i = 0; i < SAMPLES.length; i++) {
      const v = t.counts[i];
      const bg = heatColor(v, MAX_READS);
      cells += `<td class="cell" style="background:${bg}" title="${SAMPLES[i]}: ${v.toLocaleString()} reads">${v > 0 ? v.toLocaleString() : ''}</td>`;
    }

    rows.push(`<tr>
      <td class="name" title="${t.name} (taxid ${t.taxid})"><span class="rank">${t.rank}</span> ${t.name}</td>
      <td class="num">${t.rank}</td>
      <td class="num">${t.total.toLocaleString()}</td>
      <td class="num">${t.n_samp}</td>
      <td class="num"><span class="hom-bar" style="${barStyle}"></span>${hpStr}</td>
      <td class="num">${avgStr}</td>
      <td class="num">${maxStr}</td>
      ${cells}
    </tr>`);
  }
  tbody.innerHTML = rows.join('');
  document.getElementById('rowcount').textContent =
    `Showing ${filtered.length.toLocaleString()} of ${TAXA.length.toLocaleString()} taxa`;
}

function applyFilters() {
  const search   = document.getElementById('searchBox').value.toLowerCase();
  const minReads = parseInt(document.getElementById('minReads').value);
  const minSamps = parseInt(document.getElementById('minSamps').value);
  const rank     = document.getElementById('rankFilter').value;
  const blastOnly= document.getElementById('blastOnly').checked;

  filtered = TAXA.filter(t =>
    t.total >= minReads &&
    t.n_samp >= minSamps &&
    (!rank || t.rank.startsWith(rank)) &&
    (!blastOnly || t.hit_pct !== null) &&
    (search === '' || t.name.toLowerCase().includes(search))
  );
  filtered.sort((a, b) => {
    let av = a[sortKey], bv = b[sortKey];
    if (av === null) av = -Infinity;
    if (bv === null) bv = -Infinity;
    return sortDir * (bv - av || String(bv).localeCompare(String(av)));
  });
  renderTable();
}

function sortBy(key) {
  if (sortKey === key) sortDir *= -1;
  else { sortKey = key; sortDir = -1; }
  applyFilters();
}

// Initial render
applyFilters();
</script>
</body>
</html>
"""

with open(HTML_OUT, "w") as fh:
    fh.write(html)

print(f"\nHTML written to: {HTML_OUT}")
print("Done.")
