"""
07_build_pathway_csv.py  —  Step 2 of 4
Apply KEGG-style classification to the raw edges from script 06 and build
the final node/edge CSVs used by all downstream scripts.

Inputs  (from script 06):
    data/raw_edges.csv
    data/raw_master_gene_list.csv

Outputs  (verify before running scripts 08 / 09):
    data/final_nodes.csv
    data/final_edges.csv
    output/kegg_style_validation_report.txt

Intermediate checks to make:
  - final_edges.csv: kegg_edge_subtype correct for each edge?
                     repression for FBXW7 (transcription, negative)?
                     indirect effect for GPR37 / HES5 / ANGPTL2 / FYN indirect row / MAPK1 / MAPK3?
                     expression for FYN direct row (transcription, positive)?
                     repression for downstream negative edges?
                     activation for SOX10 protein-level?
  - final_nodes.csv: kegg_node_type, biological_entity_type, role columns present?
  - validation report: all 10 checks passing?
"""

import os
import sys

import csv
import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────────

DATA_DIR     = "data"
OUT_DIR      = "output"
IN_EDGES     = os.path.join(DATA_DIR, "raw_edges.csv")
IN_MASTER    = os.path.join(DATA_DIR, "raw_master_gene_list.csv")
OUT_NODES    = os.path.join(DATA_DIR, "final_nodes.csv")
OUT_EDGES    = os.path.join(DATA_DIR, "final_edges.csv")
OUT_REPORT   = os.path.join(OUT_DIR,  "kegg_style_validation_report.txt")

# ── KEGG subtype → value mapping ──────────────────────────────────────────────

SUBTYPE_VALUES: dict[str, str] = {
    "expression":          "-->",
    "repression":          "--|",
    "activation":          "-->",
    "inhibition":          "--|",
    "indirect effect":     "..>",
    "state change":        "...",
    "binding/association": "---",
    "dissociation":        "-+-",
    "phosphorylation":     "+p",
    "dephosphorylation":   "-p",
    "glycosylation":       "+g",
    "ubiquitination":      "+u",
    "methylation":         "+m",
    "others/unknown":      "?",
}

VALID_KGML_TYPES = {"GErel", "PPrel", "PCrel", "maplink"}

# ── Biological-entity-type lookup ─────────────────────────────────────────────

BIO_TYPES: dict[str, str] = {
    # Transcription factors
    "SOX10":        "transcription_factor",
    "SOX8":         "transcription_factor",
    "SOX17":        "transcription_factor",
    "HES5":         "transcription_factor",
    "EGR2":         "transcription_factor",
    "ID4":          "transcription_factor",
    "ZEB2":         "transcription_factor",
    "NKX2-2":       "transcription_factor",
    "ASCL1":        "transcription_factor",
    "MYRF_protein": "transcription_factor",
    "MYRF":         "transcription_factor",
    "OLIG1":        "transcription_factor",
    "OLIG2":        "transcription_factor",
    # Enzymes (kinases, ligases, phosphatases)
    "FYN":          "enzyme",
    "GSK3B":        "enzyme",
    "MAPK3":        "enzyme",
    "MAPK1":        "enzyme",
    "AATK":         "enzyme",
    "CDK18":        "enzyme",
    "DUSP15":       "enzyme",
    "FBXW7":        "enzyme",
    "RFFL":         "enzyme",
    "RNF220":       "enzyme",
    "PPP1R14A":     "enzyme",
    # Structural / signalling proteins
    "TMEM98":       "protein",
    "MAG":          "protein",
    "MBP":          "protein",
    "MOG":          "protein",
    "PLP1":         "protein",
    "MPZ":          "protein",
    "NFASC":        "protein",
    "GJC2":         "protein",
    "GJB1":         "protein",
    "CNTN2":        "protein",
    "TF":           "protein",
    "ANGPTL2":      "protein",
    "GPR37":        "protein",
    "LPAR1":        "protein",
    "CSPG4":        "protein",
    "PDGFRA":       "protein",
    "TGFB2":        "protein",
    "WNT7A":        "protein",
    "MOB3B":        "protein",
    "SEPT4":        "protein",
    "DBNDD2":       "protein",
    "KIF21A":       "protein",
    # Special synthetic node
    "MYRF_gene":    "gene",
}
BIO_TYPE_DEFAULT = "gene_or_protein"

POSITION_ROLE: dict[str, str] = {
    "Center":     "core_hub",
    "Upstream":   "upstream_regulator",
    "Downstream": "downstream_target",
}


# ── Helpers ───────────────────────────────────────────────────────────────────

def _s(val) -> str:
    if val is None:
        return ""
    try:
        if pd.isna(val):
            return ""
    except (TypeError, ValueError):
        pass
    return str(val).strip()


def _lower(val) -> str:
    return _s(val).lower()


# ── Edge classifier ───────────────────────────────────────────────────────────

def classify_edge(
    primary_relationship: str,
    direction: str,
    regulation_level: str,
    evidence_strength: str = "",
) -> tuple[str, str, str, bool]:
    """Return (kegg_edge_subtype, kegg_edge_value, kgml_type, needs_manual_review).

    Classification is derived strictly from the curated fields — no gene-specific
    overrides. Subtypes produced: expression, repression, activation, inhibition,
    indirect effect, others/unknown.
    """

    pr   = _lower(primary_relationship)
    lvl  = _lower(regulation_level)
    dirn = _lower(direction)

    # Explicit molecular mechanisms take precedence over generic sign/level rules.
    if "ubiquitin" in pr:
        return "ubiquitination", "+u", "PPrel", False

    if any(term in pr for term in (
        "physically interact", "binding/association", "functionally cooperates",
    )):
        return "binding/association", "---", "PPrel", False

    # Curated indirect effects are review-complete when their sign is known.
    if "indirect" in pr:
        ktype = "PPrel" if lvl == "protein" else "GErel"
        return "indirect effect", "..>", ktype, dirn not in ("positive", "negative")

    # Unknown mechanism remains visibly flagged.
    if "unknown mechanism" in pr:
        ktype = "PPrel" if lvl == "protein" else "GErel"
        return "others/unknown", "?", ktype, True

    # Rule 3: transcription + positive → expression
    if lvl == "transcription" and dirn == "positive":
        return "expression", "-->", "GErel", False

    # Rule 4: transcription + negative → repression
    if lvl == "transcription" and dirn == "negative":
        return "repression", "--|", "GErel", False

    # Rule 5: protein + positive → activation
    if lvl == "protein" and dirn == "positive":
        return "activation", "-->", "PPrel", False

    # Rule 6: protein + negative → inhibition
    if lvl == "protein" and dirn == "negative":
        return "inhibition", "--|", "PPrel", False

    # Rule 7: synthetic gene→protein edge
    if lvl == "gene":
        return "expression", "-->", "GErel", False

    # Default: unknown level or unmatched
    return "others/unknown", "?", "GErel", True


def classify_all_edges(edges: pd.DataFrame) -> pd.DataFrame:
    subtypes, values, ktypes, reviews = [], [], [], []
    for _, r in edges.iterrows():
        pr = _s(r.get("primary_relationship", ""))
        # Synthetic MYRF_gene → MYRF_protein edge
        if pr == "gene_to_protein_expression":
            subtypes.append("expression"); values.append("-->")
            ktypes.append("GErel"); reviews.append(False)
            continue
        sub, val, kt, rev = classify_edge(
            primary_relationship = pr,
            direction            = _s(r.get("direction", "")),
            regulation_level     = _s(r.get("regulation_level", "")),
            evidence_strength    = _s(r.get("confidence", "")),
        )
        subtypes.append(sub); values.append(val)
        ktypes.append(kt);    reviews.append(rev)

    out = edges.copy()
    out["kegg_edge_subtype"]   = subtypes
    out["kegg_edge_value"]     = values
    out["kgml_type"]           = ktypes
    out["needs_manual_review"] = reviews
    return out


def post_classify_dedup(edges: pd.DataFrame) -> pd.DataFrame:
    """
    Drop edges where (source, target, kegg_edge_subtype) is identical —
    can happen when two MGL rows for the same gene both classify to the
    same subtype (e.g. FYN's two rows both become 'indirect effect').
    Keep the first occurrence (strongest reference).
    """
    before = len(edges)
    out = edges.drop_duplicates(
        subset=["source", "target", "kegg_edge_subtype"], keep="first"
    ).reset_index(drop=True)
    after = len(out)
    if after < before:
        print(f"  Post-classification dedup: {before} -> {after} edges")
    return out


def finalise_edges(edges: pd.DataFrame) -> pd.DataFrame:
    def _ev_type(bucket: str) -> str:
        b = _s(bucket).lower()
        if "primary literature" in b:
            return "Primary literature"
        if "biology" in b:
            return "Biological mechanism"
        for k in ("hTFtarget", "Table S3", "UniProt", "Literature"):
            if k.lower() in b:
                return k
        return _s(bucket)

    return pd.DataFrame({
        "source":               edges["source"],
        "target":               edges["target"],
        "primary_relationship": edges["primary_relationship"],
        "kegg_edge_subtype":    edges["kegg_edge_subtype"],
        "kegg_edge_value":      edges["kegg_edge_value"],
        "kgml_type":           edges["kgml_type"],
        "direction":           edges["direction"],
        "regulation_level":    edges["regulation_level"],
        "confidence":          edges["confidence"],
        "evidence_type":       edges["source_bucket"].apply(_ev_type),
        "source_urls":         edges["source_urls"],
        "notes":               edges["notes"],
        "needs_manual_review": edges["needs_manual_review"],
    }).reset_index(drop=True)


# ── Node builder ──────────────────────────────────────────────────────────────

def build_nodes(master: pd.DataFrame, edges: pd.DataFrame) -> pd.DataFrame:
    """Build final_nodes.csv from all genes that appear in edges."""

    all_ids: set[str] = set(edges["source"]) | set(edges["target"])
    all_ids.update({"MYRF_gene", "MYRF_protein"})

    # First-occurrence metadata per gene
    first: dict[str, pd.Series] = {}
    all_refs: dict[str, list[str]] = {}
    for _, r in master.iterrows():
        sym = _s(r.get("Gene_Symbol", ""))
        if not sym:
            continue
        if sym not in first:
            first[sym] = r
        ref = _s(r.get("References", ""))
        if ref:
            all_refs.setdefault(sym, []).append(ref)

    rows = []
    for nid in sorted(all_ids):
        base = "MYRF" if nid in ("MYRF_gene", "MYRF_protein") else nid
        r = first.get(base, pd.Series(dtype=object))

        if nid == "MYRF_gene":
            display_name = "MYRF (gene)"
            role, pos = "core_hub", "Center"
        elif nid == "MYRF_protein":
            display_name = "MYRF (protein)"
            role, pos = "core_hub", "Center"
        else:
            display_name = _s(r.get("Original_Label")) or nid
            pos  = _s(r.get("Pathway_Position"))
            role = POSITION_ROLE.get(pos, "candidate_context")

        source_urls = "; ".join(all_refs.get(base, []))
        bio_type    = BIO_TYPES.get(nid) or BIO_TYPES.get(base, BIO_TYPE_DEFAULT)

        rows.append({
            "node_id":                nid,
            "display_name":           display_name,
            "original_label":         _s(r.get("Original_Label")),
            "ensembl_id":             _s(r.get("Ensembl_ID")),
            "kegg_node_type":         "gene",
            "biological_entity_type": bio_type,
            "role":                   role,
            "compartment":            "",
            "stage":                  "",
            "confidence":             _s(r.get("Evidence_Strength")),
            "source_urls":            source_urls,
            "notes":                  _s(r.get("Notes")),
            "needs_manual_review":    False,
        })

    return pd.DataFrame(rows)


# ── Validation ────────────────────────────────────────────────────────────────

def validate(nodes: pd.DataFrame, edges: pd.DataFrame, lines: list[str]) -> int:
    errors, warnings = [], []
    nids     = set(nodes["node_id"])
    valid_st = set(SUBTYPE_VALUES.keys())

    # 1. Valid kegg_node_type
    bad = nodes[~nodes["kegg_node_type"].isin({"gene", "group", "compound", "map"})]
    for _, r in bad.iterrows():
        errors.append(f"[1] node '{r['node_id']}': invalid kegg_node_type='{r['kegg_node_type']}'")

    # 2. Valid kegg_edge_subtype
    bad = edges[~edges["kegg_edge_subtype"].isin(valid_st)]
    for _, r in bad.iterrows():
        errors.append(f"[2] edge {r['source']}->{r['target']}: invalid subtype='{r['kegg_edge_subtype']}'")

    # 3. kegg_edge_value matches subtype
    for _, r in edges.iterrows():
        expected = SUBTYPE_VALUES.get(_s(r["kegg_edge_subtype"]))
        if expected and _s(r["kegg_edge_value"]) != expected:
            errors.append(f"[3] edge {r['source']}->{r['target']}: value mismatch "
                          f"(subtype={r['kegg_edge_subtype']}, "
                          f"got={r['kegg_edge_value']}, expected={expected})")

    # 4. Valid kgml_type
    bad = edges[~edges["kgml_type"].isin(VALID_KGML_TYPES)]
    for _, r in bad.iterrows():
        errors.append(f"[4] edge {r['source']}->{r['target']}: invalid kgml_type='{r['kgml_type']}'")

    # 5. No direction=negative with subtype=expression
    bad = edges[(edges["direction"].str.lower() == "negative") &
                (edges["kegg_edge_subtype"] == "expression")]
    for _, r in bad.iterrows():
        errors.append(f"[5] edge {r['source']}->{r['target']}: negative+expression conflict")

    # 6. No direction=positive with inhibition/repression
    bad = edges[(edges["direction"].str.lower() == "positive") &
                (edges["kegg_edge_subtype"].isin({"inhibition", "repression"}))]
    for _, r in bad.iterrows():
        errors.append(f"[6] edge {r['source']}->{r['target']}: positive+{r['kegg_edge_subtype']} conflict")

    # 7. Strong edges have source_urls (excl. synthetic)
    strong = edges[edges["confidence"].str.lower() == "strong"]
    no_url = strong[
        strong["source_urls"].isna() | (strong["source_urls"].str.strip() == "")
    ]
    no_url = no_url[~((no_url["source"] == "MYRF_gene") & (no_url["target"] == "MYRF_protein"))]
    for _, r in no_url.iterrows():
        warnings.append(f"[7] Strong edge {r['source']}->{r['target']} has no source_urls")

    # 8. Edge endpoints in node table
    for col in ("source", "target"):
        missing = edges[~edges[col].isin(nids)][col].unique()
        for m in missing:
            errors.append(f"[8] edge {col} '{m}' not in node table")

    # 9. MYRF_gene and MYRF_protein present
    for req in ("MYRF_gene", "MYRF_protein"):
        if req not in nids:
            errors.append(f"[9] required node '{req}' missing")

    # 10. Synthetic edge present
    if not ((edges["source"] == "MYRF_gene") & (edges["target"] == "MYRF_protein")).any():
        errors.append("[10] synthetic MYRF_gene->MYRF_protein edge missing")

    all_issues = errors + warnings
    for line in all_issues:
        lines.append(line)
    if not all_issues:
        lines.append("All 10 validation checks PASSED.")
    return len(errors)


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> int:
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(OUT_DIR,  exist_ok=True)

    for path in (IN_EDGES, IN_MASTER):
        if not os.path.exists(path):
            sys.exit(f"ERROR: input file not found: {path}\n"
                     "Run scripts/06_excel_to_raw.py first.")

    print(f"Reading: {IN_EDGES}")
    raw_edges = pd.read_csv(IN_EDGES)
    print(f"  {len(raw_edges)} raw edges")

    print(f"Reading: {IN_MASTER}")
    master = pd.read_csv(IN_MASTER)
    print(f"  {len(master)} master gene list rows")

    # Classify
    classified = classify_all_edges(raw_edges)
    classified = post_classify_dedup(classified)

    # Build nodes from master metadata + edge gene list
    nodes_df   = build_nodes(master, classified)
    final_edges = finalise_edges(classified)
    print(f"  {len(nodes_df)} nodes built")
    print(f"  {len(final_edges)} edges after classification")

    # Validate
    report_lines: list[str] = [
        "KEGG-style Pathway Validation Report",
        "=" * 62, "",
    ]
    n_failures = validate(nodes_df, final_edges, report_lines)
    report_lines.append("")

    print("\n" + ("Validation PASSED." if n_failures == 0
                  else f"Validation FAILED — {n_failures} error(s)."))

    # Summary
    sep = "=" * 62
    summary = [sep, "SUMMARY", sep, ""]
    summary.append("Nodes by kegg_node_type:")
    for k, c in nodes_df["kegg_node_type"].value_counts().items():
        summary.append(f"  {k}: {c}")
    summary.append("\nEdges by kegg_edge_subtype:")
    for k, c in final_edges["kegg_edge_subtype"].value_counts().items():
        summary.append(f"  {k}: {c}")
    summary.append("\nEdges by kgml_type:")
    for k, c in final_edges["kgml_type"].value_counts().items():
        summary.append(f"  {k}: {c}")
    review_n = int(final_edges["needs_manual_review"].sum())
    summary.append(f"\nneeds_manual_review: {review_n} edges")
    if review_n:
        rev_edges = final_edges[final_edges["needs_manual_review"]]
        for _, r in rev_edges.iterrows():
            summary.append(f"  {r['source']} -> {r['target']}  [{r['kegg_edge_subtype']}]")
    summary.append(sep)
    for line in summary:
        print(line)
    report_lines.extend(summary)

    # Write
    nodes_df.to_csv(OUT_NODES,  index=False, encoding="utf-8", quoting=csv.QUOTE_NONNUMERIC)
    final_edges.to_csv(OUT_EDGES, index=False, encoding="utf-8", quoting=csv.QUOTE_NONNUMERIC)
    with open(OUT_REPORT, "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines) + "\n")

    print(f"\n[OK] {OUT_NODES}")
    print(f"[OK] {OUT_EDGES}")
    print(f"[OK] {OUT_REPORT}")
    print("\nVerify before running scripts/03_preview_kegg_graph.py")
    return 0 if n_failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
