# MYRF pathway curation workflow

Your proposed workflow is broadly correct. The main corrections are:

- Define the pathway scope and evidence rules before searching, so inclusion decisions are consistent.
- Record experimental context and uncertainty, not only key sentences and the inferred relationship.
- Treat the Excel workbook as the source of truth; generated CSV, SVG, DOT, and KGML files should not be edited manually.
- Generate and inspect the graph before exporting KGML/XML. In this repository, the KGML exporter can reuse the DOT layout.
- Audit twice: first for biological interpretation, then after regeneration for literature and artifact consistency.

## End-to-end workflow

### 1. Define the scope and inclusion rules

Before collecting papers, state the biological question and what qualifies as a pathway edge. For this project, an edge should be relevant to MYRF-centered oligodendrocyte differentiation or myelination and supported by traceable evidence.

Prefer primary research. Reviews and database pages can help discover literature, but they should not be the sole support for a retained edge when a primary study is available.

### 2. Find and screen relevant literature

Search for studies involving MYRF regulators, MYRF protein processing or stability, MYRF binding targets, transcriptional effects, and interacting transcription factors or signaling pathways.

Screen each paper for:

- organism, tissue, and cell type;
- developmental or disease context;
- perturbation used, such as knockout, knockdown, overexpression, inhibitor, or mutation;
- assay type, such as ChIP-seq, reporter assay, RNA expression, protein abundance, binding, or ubiquitination;
- whether the result supports direct regulation, an indirect effect, association only, or an unresolved mechanism.

### 3. Extract the evidence

Record the key sentence or result together with enough context to interpret it correctly. Capture:

- citation and stable URL or DOI;
- relevant sentence, figure, table, or supplementary-data location;
- source and target entities;
- measured outcome;
- organism and experimental system;
- assay and perturbation;
- limitations or alternative interpretations.

Do not turn correlation, co-expression, or pathway membership into a directional regulatory edge unless the experiment supports that inference.

### 4. Derive the pathway relationship

Translate the evidence into explicit edge fields:

- source and target;
- positive, negative, or neutral direction;
- gene/transcription or protein regulation level;
- direct, indirect, binding/association, inhibition, repression, ubiquitination, or unknown mechanism;
- Strong, Moderate, or Weak confidence;
- whether manual review is still needed.

Keep gene-locus regulation separate from protein abundance, processing, binding, and degradation. In particular, this project uses `MYRF_gene` for transcriptional regulation and `MYRF_protein` for active protein-level relationships.

When a study tests a complex or pathway jointly, do not claim independently resolved component effects. Preserve this limitation in the notes and confidence level.

### 5. Log the curation in the workbook

Enter the result in `raw_data/myrf-final-consolidation-documentation.xlsx`, which is the source of truth.

- `Master_Gene_List`: curated source rows and relationship fields.
- `References`: portable citation list.
- `Revision_Log`: material corrections and their rationale.
- `Literature_Audit`: one row for every retained graph edge and its primary evidence.

Retain rejected or superseded rows as `Excluded` with a reason instead of silently deleting them. This preserves the audit trail.

### 6. Perform the first biological audit

Before generating artifacts, check the workbook for:

- unsupported directness or direction;
- gene-versus-protein routing errors;
- duplicated observations represented as multiple edges;
- citations that do not support the stated mechanism;
- confidence levels that overstate the experimental evidence;
- identifier errors and institution-only proxy URLs.

Correct the workbook first. Do not patch generated CSV or XML files as substitutes for fixing the source data.

### 7. Convert the workbook to normalized raw tables

Run:

```powershell
python scripts/01_excel_to_raw.py
```

This reads the workbook, normalizes the master sheet, merges true duplicates without losing citations, splits MYRF gene/protein routing, and writes:

- `data/raw_master_gene_list.csv`
- `data/raw_edges.csv`

### 8. Classify and validate pathway nodes and edges

Run:

```powershell
python scripts/02_build_pathway_csv.py
```

This assigns KEGG-style node and relationship semantics, validates the model, and writes:

- `data/final_nodes.csv`
- `data/final_edges.csv`
- `output/kegg_style_validation_report.txt`

Stop and resolve any validation or manual-review finding before continuing.

### 9. Generate and inspect the graph

Run:

```powershell
python scripts/03_preview_kegg_graph.py
```

This creates:

- `output/myrf_final_pathway.dot`
- `output/myrf_final_pathway.svg`

Inspect the SVG for missing nodes, wrong directions, incorrect inhibitory arrowheads, duplicate edges, misleading labels, overlap, and unreadable layout. The SVG is the authoritative rendered pathway for this custom model.

### 10. Export KGML/XML and the evidence table

After the graph layout is available, run:

```powershell
python scripts/04_export_kgml.py
```

This creates:

- `output/myrf_pathway.kgml`
- `output/myrf_pathway_evidence.csv`
- `output/myrf_pathway_entry_id_map.csv`
- `custom_kegg/hsa99999.xml`

The exporter automatically keeps `custom_kegg/hsa99999.xml` byte-identical to the current KGML output.

This is a custom KGML-like pathway, not an official KEGG pathway. Pathview compatibility is therefore optional and separate from biological correctness.

### 11. Perform the second-round audit

Run:

```powershell
python scripts/07_audit_literature.py
```

The second audit checks that:

- every retained edge exactly matches the literature-audit register;
- every edge has a supporting primary citation;
- excluded overstatements remain absent;
- relationship descriptions and review states survive export;
- node, edge, subtype, SVG, and KGML counts agree;
- KGML relation references resolve;
- graph styles preserve biological direction;
- both KGML copies are synchronized;
- workbook audit sheets and identifiers are valid.

The result is written to `output/literature_audit_round2.md`.

### 12. Complete manual verification and release

Even after the automated audit passes, manually verify:

- each literature interpretation against the cited passage, figure, or table;
- the final SVG visually;
- the workbook audit and revision sheets;
- the Git diff, ensuring generated artifacts changed only because their source data or generator changed.

If any problem is found, return to the workbook, correct the source row, regenerate all downstream artifacts, and rerun the second audit.

Only commit or publish when the biological review, generated-artifact validation, and manual visual check all pass.

## Condensed workflow

```text
Define scope and evidence rules
  -> find and screen primary literature
  -> extract evidence plus experimental context
  -> derive a qualified relationship
  -> log it in the source workbook
  -> first biological audit and correction
  -> workbook to raw CSV
  -> classify and validate nodes/edges
  -> generate and visually inspect DOT/SVG
  -> export and synchronize KGML/XML plus evidence
  -> second literature/artifact audit
  -> final manual verification
  -> commit or publish
```
