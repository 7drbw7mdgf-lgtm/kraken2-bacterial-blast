#!/usr/bin/env python3
"""
Generates a compact interactive HTML table:
  Rows    = bacterial species (from Kraken2 reports)
  Columns = sample names (96 samples)
  Cells   = % of reads for that species that matched C. jejuni in BLASTn
            Coloured white→orange→red by %; grey = species absent from sample
Data is embedded as compact JSON; rows are rendered by JS (virtual scroll).
"""

import os, json
from collections import defaultdict

REPORT_DIR = "/home/lshpk18/unmapped reads (kraken)/kraken2_output"
KRAKEN_OUT = os.path.join(REPORT_DIR, "kraken2.output")
BLAST6     = "/home/lshpk18/unmapped reads (kraken)/bacterial_vs_cjejuni.blast6"
HTML_OUT   = "/home/lshpk18/unmapped reads (kraken)/kraken2_cjejuni_homology.html"
CONDITION_ORDER = ["uninfected_0h", "15m", "30m", "1h", "3h", "24h"]

# ── 1. Parse reports ──────────────────────────────────────────────────────────
print("Parsing Kraken2 reports ...")
report_files = sorted(f for f in os.listdir(REPORT_DIR) if f.endswith(".kraken2.report"))
samples      = [f.replace(".kraken2.report", "") for f in report_files]

taxid_meta   = {}        # taxid → {name, rank}
sample_counts = {}       # sample → {taxid: direct_reads}

for rfile, sname in zip(report_files, samples):
    counts = {}
    with open(os.path.join(REPORT_DIR, rfile)) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6: continue
            try: d, tid = int(p[2]), int(p[4])
            except ValueError: continue
            if tid not in taxid_meta:
                taxid_meta[tid] = {"name": p[5].strip(), "rank": p[3]}
            if d > 0:
                counts[tid] = d
    sample_counts[sname] = counts

# ── 2. Load kraken2.output ────────────────────────────────────────────────────
print("Loading kraken2.output ...")
read_to_taxid = {}
with open(KRAKEN_OUT) as fh:
    for line in fh:
        p = line.split("\t")
        if len(p) < 3 or p[0] != "C": continue
        try: read_to_taxid[p[1]] = int(p[2])
        except ValueError: pass

taxid_n_reads = defaultdict(int)
for tid in read_to_taxid.values():
    taxid_n_reads[tid] += 1

# ── 3. Parse BLAST6 ───────────────────────────────────────────────────────────
print("Computing BLAST homology ...")
read_best = {}
with open(BLAST6) as fh:
    for line in fh:
        c = line.rstrip("\n").split("\t")
        if len(c) < 12: continue
        qb = c[0].rsplit("/", 1)[0]
        pi = float(c[2])
        if qb not in read_best or pi > read_best[qb]:
            read_best[qb] = pi

taxid_hits = defaultdict(int)
taxid_pidsum = defaultdict(float)
for rid, pi in read_best.items():
    tid = read_to_taxid.get(rid)
    if tid is not None:
        taxid_hits[tid]   += 1
        taxid_pidsum[tid] += pi

def get_homology(tid):
    n = taxid_n_reads.get(tid, 0)
    h = taxid_hits.get(tid, 0)
    if n == 0:
        return None, None
    hit_pct = round(h / n * 100, 1)
    avg_pid = round(taxid_pidsum[tid] / h, 1) if h else 0.0
    return hit_pct, avg_pid

# ── 4. Filter to species-level with ≥5 total reads ───────────────────────────
taxid_total = defaultdict(int)
taxid_nsamp = defaultdict(int)
for sc in sample_counts.values():
    for tid, cnt in sc.items():
        taxid_total[tid] += cnt
        taxid_nsamp[tid] += 1

rows = []
for tid, meta in taxid_meta.items():
    if not meta["rank"].startswith("S"): continue
    total = taxid_total[tid]
    if total < 5: continue
    hp, avg = get_homology(tid)
    rows.append({
        "taxid":  tid,
        "name":   meta["name"],
        "total":  total,
        "nsamp":  taxid_nsamp[tid],
        "hp":     hp,     # % reads hitting C. jejuni (None if not in blast sample)
        "avg":    avg,
    })

rows.sort(key=lambda x: (-(x["hp"] or -1), -x["total"]))
print(f"  {len(rows)} species for display")

# ── 5. Sort samples by condition ──────────────────────────────────────────────
def cond(s):
    for c in CONDITION_ORDER:
        if s.startswith(c): return c
    return "other"

samples_sorted = []
cond_groups = {}
for c in CONDITION_ORDER:
    grp = sorted(s for s in samples if cond(s) == c)
    cond_groups[c] = grp
    samples_sorted.extend(grp)

# ── 6. Build compact JSON payload ────────────────────────────────────────────
# species_data: list of [name, total, nsamp, hp_or_null, avg_or_null, [counts per sample]]
species_data = []
for sp in rows:
    tid = sp["taxid"]
    counts = [sample_counts[s].get(tid, 0) for s in samples_sorted]
    species_data.append([
        sp["name"],
        sp["total"],
        sp["nsamp"],
        sp["hp"],
        sp["avg"],
        counts,
    ])

cond_info = []
for c in CONDITION_ORDER:
    grp = cond_groups[c]
    if grp:
        label = c.replace("uninfected_0h","Uninfected 0h").replace("_"," ")
        cond_info.append({"label": label, "n": len(grp)})

sample_labels = [s.replace("uninfected_0h_","U0h ").replace("_"," ") for s in samples_sorted]

js_species  = json.dumps(species_data, separators=(',',':'))
js_samples  = json.dumps(sample_labels,  separators=(',',':'))
js_conds    = json.dumps(cond_info,      separators=(',',':'))

# ── 7. Write HTML ─────────────────────────────────────────────────────────────
print("Writing HTML ...")

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Bacterial Species vs C. jejuni Homology</title>
<style>
*{{box-sizing:border-box;margin:0;padding:0;font-family:'Segoe UI',Arial,sans-serif;font-size:12px}}
body{{background:#f4f6f9}}
header{{background:#1a3a5c;color:#fff;padding:14px 20px}}
header h1{{font-size:17px;font-weight:600}}
header p{{font-size:11px;opacity:.75;margin-top:3px;line-height:1.5}}
.bar{{display:flex;flex-wrap:wrap;gap:12px;align-items:center;background:#fff;padding:9px 20px;border-bottom:1px solid #dde}}
.bar label{{font-weight:600}}
.bar input[type=text]{{border:1px solid #ccc;border-radius:4px;padding:3px 8px;width:200px}}
.bar input[type=range]{{width:130px;vertical-align:middle}}
.sv{{font-weight:700;color:#1a3a5c;min-width:34px;display:inline-block}}
.leg{{display:flex;align-items:center;gap:5px;font-size:11px}}
.lg{{width:110px;height:11px;border-radius:3px;background:linear-gradient(to right,#fff,#ff8800,#c00);border:1px solid #ccc}}
.lg2{{width:18px;height:11px;background:#ddd;border-radius:2px;border:1px solid #ccc}}
#info{{padding:5px 20px;background:#f0f0f0;border-bottom:1px solid #dde;color:#555;font-size:11px}}
.wrap{{overflow:auto;height:calc(100vh - 175px)}}
table{{border-collapse:collapse;table-layout:fixed}}
thead th{{position:sticky;top:0;background:#1a3a5c;color:#fff;padding:5px 7px;
          white-space:nowrap;cursor:pointer;user-select:none;z-index:3;
          border-right:1px solid #2a5a8c;font-weight:600}}
thead th:hover{{background:#2a5a8c}}
.ch{{background:#0e2a45 !important;text-align:center;font-size:10px;letter-spacing:.3px;
     border-right:2px solid #1a3a5c !important}}
.sf{{position:sticky;z-index:4}}
.sf1{{left:0;min-width:220px;max-width:220px}}
.sf2{{left:220px;min-width:75px;max-width:75px}}
.sf3{{left:295px;min-width:65px;max-width:65px}}
.sf4{{left:360px;min-width:55px;max-width:55px;border-right:2px solid #2a5a8c !important}}
tbody td{{padding:2px 6px;border-bottom:1px solid #eee;border-right:1px solid #eee;white-space:nowrap}}
tbody tr:hover td{{filter:brightness(.92)}}
.tn{{overflow:hidden;text-overflow:ellipsis;background:#fff}}
.tn.sf1{{border-right:1px solid #ccc}}
.tc{{text-align:right;background:#f7f7f7;font-variant-numeric:tabular-nums}}
.tc.sf4{{border-right:2px solid #ccc}}
.cell{{text-align:center;min-width:52px;max-width:52px;font-weight:600}}
.absent{{background:#e0e0e0 !important;color:#aaa;font-weight:400}}
.nt{{color:#888;font-weight:400}}
</style>
</head>
<body>
<header>
  <h1>Bacterial Species &mdash; <em>C. jejuni</em> NCTC 11168 BLAST % Homology</h1>
  <p>Rows = species from 96 Kraken2 reports &nbsp;|&nbsp;
     Columns = samples grouped by timepoint &nbsp;|&nbsp;
     Cell = % reads matching <em>C. jejuni</em> (BLASTn, e&#8209;value&nbsp;&le;&nbsp;1e&#8209;5, id&nbsp;&ge;&nbsp;70%) &nbsp;|&nbsp;
     <span style="background:#ddd;padding:0 4px;border-radius:2px">grey</span> = species absent &nbsp;|&nbsp;
     N/T = present but not in BLAST sample
  </p>
</header>
<div class="bar">
  <span><label>Search: </label>
    <input type="text" id="q" placeholder="species name…" oninput="applyFilter()"></span>
  <span><label>Min total reads: </label>
    <input type="range" id="mr" min="5" max="5000" step="5" value="5"
      oninput="document.getElementById('mrv').textContent=this.value;applyFilter()">
    <span class="sv" id="mrv">5</span></span>
  <span><label>Min samples: </label>
    <input type="range" id="ms" min="1" max="96" value="1"
      oninput="document.getElementById('msv').textContent=this.value;applyFilter()">
    <span class="sv" id="msv">1</span></span>
  <span class="leg">0%<div class="lg"></div>100%&nbsp;
    <div class="lg2"></div>absent</span>
</div>
<div id="info"></div>
<div class="wrap" id="wrap">
<table id="tbl">
<thead>
<tr id="condRow"></tr>
<tr id="hdRow"></tr>
</thead>
<tbody id="tb"></tbody>
</table>
</div>

<script>
const SP   = {js_species};
const SAMP = {js_samples};
const COND = {js_conds};

// Build header
(function(){{
  const cr = document.getElementById('condRow');
  const hr = document.getElementById('hdRow');

  // Fixed col headers
  const fixedConds = ['','','',''];
  ['Species','Total','#Samp','Avg ID%'].forEach((t,i)=>{{
    const th = document.createElement('th');
    th.rowSpan=2; th.textContent=t;
    th.className='sf sf'+(i+1);
    if(i===0) th.onclick=()=>sortBy(0);
    if(i===1) th.onclick=()=>sortBy(1);
    if(i===2) th.onclick=()=>sortBy(2);
    if(i===3) th.onclick=()=>sortBy(4);
    cr.appendChild(th);
  }});

  // Condition spans
  COND.forEach(c=>{{
    const th=document.createElement('th');
    th.colSpan=c.n; th.textContent=c.label; th.className='ch';
    cr.appendChild(th);
  }});

  // Sample headers
  SAMP.forEach((s,i)=>{{
    const th=document.createElement('th');
    th.textContent=s; th.title=s;
    th.style.fontSize='10px'; th.style.fontWeight='400';
    th.style.minWidth='52px'; th.style.maxWidth='52px';
    hr.appendChild(th);
  }});
}})();

function pctColor(p){{
  if(p===null) return null; // not tested
  if(p===0) return '#ffffff';
  const t=p/100;
  const r=255, g=Math.round(255*(1-t)), b=Math.round(220*(1-t));
  return `rgb(${{r}},${{g}},${{b}})`;
}}
function textCol(p){{ return (p!==null&&p>=55)?'#111':'#333'; }}

let sortCol=1, sortDir=-1;
let visible=[...SP];

function sortBy(col){{
  if(sortCol===col) sortDir*=-1; else{{sortCol=col;sortDir=-1;}}
  applyFilter();
}}

function applyFilter(){{
  const q  =document.getElementById('q').value.toLowerCase();
  const mr =parseInt(document.getElementById('mr').value);
  const ms =parseInt(document.getElementById('ms').value);
  visible=SP.filter(r=>
    r[1]>=mr && r[2]>=ms && (q===''||r[0].toLowerCase().includes(q))
  );
  visible.sort((a,b)=>{{
    let av=a[sortCol], bv=b[sortCol];
    if(av===null) av=-999; if(bv===null) bv=-999;
    return sortDir*(bv>av?1:bv<av?-1:0);
  }});
  render();
}}

function render(){{
  const tb=document.getElementById('tb');
  const frags=[];
  for(const r of visible){{
    const[name,total,nsamp,hp,avg,counts]=r;
    const avgStr = avg!==null ? avg.toFixed(1)+'%' : '<span class="nt">N/T</span>';
    let cells='';
    for(let i=0;i<SAMP.length;i++){{
      const cnt=counts[i];
      if(cnt===0){{
        cells+=`<td class="cell absent" title="${{SAMP[i]}}: absent">—</td>`;
      }}else{{
        const bg=pctColor(hp);
        const valStr=hp!==null ? hp.toFixed(1)+'%' : '<span class="nt">N/T</span>';
        const tip=`${{SAMP[i]}}: ${{cnt.toLocaleString()}} reads | ${{hp!==null?hp.toFixed(1)+'% hit C.jejuni':'not in BLAST sample'}}`;
        if(bg){{
          cells+=`<td class="cell" style="background:${{bg}};color:${{textCol(hp)}}" title="${{tip}}">${{valStr}}</td>`;
        }}else{{
          cells+=`<td class="cell" title="${{tip}}">${{valStr}}</td>`;
        }}
      }}
    }}
    frags.push(`<tr>
      <td class="tn sf sf1" title="${{name}}">${{name}}</td>
      <td class="tc sf sf2">${{total.toLocaleString()}}</td>
      <td class="tc sf sf3">${{nsamp}}</td>
      <td class="tc sf sf4">${{avgStr}}</td>
      ${{cells}}
    </tr>`);
  }}
  tb.innerHTML=frags.join('');
  document.getElementById('info').textContent=
    `Showing ${{visible.length.toLocaleString()}} of ${{SP.length.toLocaleString()}} species · ${{SAMP.length}} samples`;
}}

applyFilter();
</script>
</body></html>
"""

with open(HTML_OUT, "w") as fh:
    fh.write(html)
print(f"Written: {HTML_OUT}  ({os.path.getsize(HTML_OUT)/1024:.0f} KB)")
