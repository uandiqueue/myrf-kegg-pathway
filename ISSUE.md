# MYRF Pathway — Open Issues & Verification Checklist

Generated: 2026-05-26 | Last updated: 2026-05-28  
Source audit: `raw_master_gene_list.csv` → `raw_edges.csv` → `final_edges.csv` → SVG

---

## Issue 1 — ~~CRITICAL~~ RESOLVED: FYN direct-regulation edge silently dropped

### What happened
Rule 10 in the classifier treated all FYN+transcription edges as `indirect effect`, causing both
FYN rows to collapse to the same subtype. Post-classification deduplication then dropped one.

### Fix applied (2026-05-28)
Removed `KINASE_INDIRECT` rule from `02_build_pathway_csv.py`. Classification is now driven
purely by the curated fields. Result: FYN now has **two distinct edges** in `final_edges.csv`:

| Edge | Subtype | Basis |
|---|---|---|
| FYN → MYRF_gene | `expression` (-->) | "Regulates MYRF gene", transcription, positive |
| FYN → MYRF_gene | `indirect effect` (..>) | "Indirectly regulates MYRF gene" |

### Still to verify (optional)
- [ ] Confirm both FYN mechanisms are genuinely distinct in the literature.
  Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC9976406/#Sec2

---

## Issue 2 — MODERATE: FBXW7 — classification fixed, mechanism still needs verification

### What happened
Rule 1 hard-overrode the curated `regulation_level=transcription` and classified FBXW7 as
`ubiquitination (+u, PPrel)` based on its known biology as an E3 ligase.

### Fix applied (2026-05-28)
Removed the gene-specific Rule 1 override. FBXW7 now classifies from curated fields:

| Field | Value | Result |
|---|---|---|
| regulation_level | transcription | GErel |
| direction | negative | `repression` (--|) |
| needs_manual_review | False | no flag |

### Still to verify (required)
The mechanism question is **not resolved** — only the classification inconsistency is fixed.

- [ ] **Does FBXW7 act at the protein level or transcription level?**
  - Protein degradation (ubiquitinates MYRF protein) → correct node should be `MYRF_protein`,
    subtype `inhibition`, type `PPrel`. Update `regulation_level=protein` in raw Excel.
  - Transcriptional repression (represses MYRF gene expression) → current classification is correct.
  - Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC12370921/

---

## Issue 3 — MODERATE (open): Direction of indirect effects invisible in SVG

### What happened
`03_preview_kegg_graph.py` renders all `indirect effect` edges identically — grey dashed arrow —
regardless of `direction`. Inhibitory indirect effects cannot be distinguished visually.

**Affected edges:**

| Edge | direction | SVG appearance |
|---|---|---|
| GPR37 → MYRF_gene | **negative** | grey dashed → (same as positive) |
| HES5 → MYRF_gene | **negative** | grey dashed → (same as positive) |
| MAG → MYRF_gene | positive | grey dashed → |
| ANGPTL2 → MYRF_gene | positive | grey dashed → |
| FYN → MYRF_gene (indirect) | positive | grey dashed → |
| MAPK1/3 → MYRF_gene | positive | grey dashed → |

### What you need to verify

- [ ] **GPR37 direction**: The note says "absence of GPR37 → enhanced ERK1/2 → decreased Hes5 →
  increased Myrf." This is a double-negative — GPR37 suppresses a pathway that promotes MYRF.
  Confirm net effect: is GPR37 a negative regulator of MYRF (direction=negative correct)?
  Source: https://www.nature.com/articles/ncomms10884
- [ ] **HES5 direction**: Confirm HES5 is a negative indirect regulator (sequesters SOX10).
  Source: https://www-nature-com.libproxy1.nus.edu.sg/articles/ncomms10884#Sec1

### Fix (code side — no data change needed once directions confirmed)
Update `03_preview_kegg_graph.py`: for `indirect effect` edges where `direction=negative`,
use a red/dark colour or `tee` arrowhead to distinguish from activating indirect edges.

---

## Issue 4 — MINOR (open): PLP1 missing second source URL (Wmn1 enhancer reference)

### What happened
`raw_master_gene_list.csv` has **two PLP1 rows** with different sources:

| Row | Source bucket | URL |
|---|---|---|
| Row 24 | Table S3 (MYRF ChIP peaks; rat) | plos.org/plosbiology + pubmed/34230963 |
| Row 50 | Various literatures (Wmn1 enhancer) | **pmc.ncbi.nlm.nih.gov/articles/PMC3814293/** |

Deduplication in `01_excel_to_raw.py` kept only the first row's URLs. The Wmn1 reference
(PMC3814293) is absent from all downstream files.

### What you need to verify

- [ ] **Is PMC3814293 a distinct, important reference for PLP1 regulation by MYRF?** If yes,
  the URL should be merged into the edge's `source_urls`.
- [ ] **Is the Wmn1 enhancer mechanistically distinct** from the Table S3 ChIP peak?

### Fix (once verified)
Update `01_excel_to_raw.py` to concatenate `source_urls` from all rows sharing the same
(source, target) instead of keeping only the first. Re-run pipeline from step 1.

---

## Issue 5 — INFO (open): `primary_relationship` not carried to final_edges.csv

### What happened
`raw_edges.csv` carries a human-readable `primary_relationship` column (e.g. "Indirectly regulates
MYRF gene", "Gene regulated by MYRF") that is dropped by `02_build_pathway_csv.py`. Not preserved
in the KGML, SVG tooltips, or evidence CSV.

### Action (optional)
No biological verification needed. If you want the description preserved downstream, ask to add a
`relationship_description` column in `finalise_edges()` in `02_build_pathway_csv.py`.

---

## Quick-reference: verification sources

| Gene | Key source to check |
|---|---|
| FYN | https://pmc.ncbi.nlm.nih.gov/articles/PMC9976406/#Sec2 |
| FBXW7 | https://pmc.ncbi.nlm.nih.gov/articles/PMC12370921/ |
| GPR37 | https://www.nature.com/articles/ncomms10884 |
| HES5 | https://www-nature-com.libproxy1.nus.edu.sg/articles/ncomms10884#Sec1 |
| PLP1 (Wmn1) | https://pmc.ncbi.nlm.nih.gov/articles/PMC3814293/ |

---

## Status tracker

| Issue | Priority | Your action | Code fix needed | Status |
|---|---|---|---|---|
| 1 — FYN direct edge | ~~CRITICAL~~ | Optionally confirm 2 mechanisms | Done | **Resolved** |
| 2 — FBXW7 mechanism | MODERATE | Verify protein vs. transcription level | Done (if data changes) | Classification fixed; **mechanism open** |
| 3 — Indirect direction in SVG | MODERATE | Confirm GPR37 + HES5 direction | Yes (script 03 colours) | Open |
| 4 — PLP1 Wmn1 URL | MINOR | Confirm if URL is needed | Yes (script 01 dedup) | Open |
| 5 — primary_relationship dropped | INFO | Decide if needed | Yes (script 02 output) | Open |
