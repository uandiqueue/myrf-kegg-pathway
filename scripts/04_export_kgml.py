"""
04_export_kgml.py
Converts curated pathway CSV files into a KGML-like XML file.

Inputs  (preferred)  : data/final_nodes.csv, data/final_edges.csv
        (fallback)   : data/nodes.csv, data/edges.csv
Layout  (optional)   : output/myrf_final_pathway.dot  (or myrf_pathway.dot)

Outputs:
  output/myrf_pathway.kgml
  output/myrf_pathway_evidence.csv
  output/myrf_pathway_entry_id_map.csv

MYRF is represented as two KGML entries (gene locus + active TF protein).
If final_nodes.csv already contains MYRF_gene and MYRF_protein (produced by
06_build_kegg_pathway.py), the split step is skipped.  If only MYRF exists,
the split is performed here.
Edges are routed to the appropriate MYRF entry based on regulation_level.
kegg_edge_subtype and kgml_type are read directly from the edges CSV when
available; otherwise they are inferred from interaction/level fields.
"""

import csv
import os
import shutil
import sys
import subprocess
from xml.dom import minidom
from xml.etree import ElementTree as ET

import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────────

OUT_DIR    = "output"
OUT_KGML   = os.path.join(OUT_DIR, "myrf_pathway.kgml")
OUT_EV     = os.path.join(OUT_DIR, "myrf_pathway_evidence.csv")
OUT_MAP    = os.path.join(OUT_DIR, "myrf_pathway_entry_id_map.csv")
OUT_CUSTOM_KGML = os.path.join("custom_kegg", "hsa99999.xml")

NODE_W, NODE_H = 80, 30
KGML_W, KGML_H, MARGIN = 1200, 900, 80   # KGML canvas dimensions in px

# Separation (px) between MYRF_gene and MYRF_protein centres in the final canvas
MYRF_SPLIT_OFFSET = 25

# ── Entrez ID map (confirmed via org.Hs.eg.db, 2026-05-21) ───────────────────
# Used to write  hsa:ENTREZID  entry names so Pathview can map gene data.
# MYRF_gene / MYRF_protein both map to MYRF (hsa:745).
# Genes not present here fall back to  custom:SYMBOL.
ENTREZ_MAP: dict[str, str] = {
    "AATK":    "9625",
    "CNTN2":   "6900",
    "CSPG4":   "1464",
    "DUSP15":  "128853",
    "EGR2":    "1959",
    "FBXW7":   "55294",
    "FYN":     "2534",
    "GJB1":    "2705",
    "GJC2":    "57165",
    "GSK3B":   "2932",
    "ID4":     "3400",
    "MAG":     "4099",
    "MAPK1":   "5594",
    "MAPK3":   "5595",
    "MBP":     "4155",
    "MOG":     "4340",
    "MPZ":     "4359",
    "MYRF":    "745",
    "MYRF_gene":    "745",   # synthetic split; maps back to MYRF Entrez
    "MYRF_protein": "745",   # synthetic split; maps back to MYRF Entrez
    "NFASC":   "23114",
    "PDGFRA":  "5156",
    "PLP1":    "5354",
    "RFFL":    "117584",
    "SOX10":   "6663",
    "SOX8":    "30812",
    "TF":      "7018",
    "TGFB2":   "7042",
    "TMEM98":  "26022",
    "WNT7A":   "7476",
}


# ── Helpers ───────────────────────────────────────────────────────────────────

def _s(val) -> str:
    return "" if pd.isna(val) else str(val).strip()


def _pick_files():
    pairs = [
        ("data/final_nodes.csv", "data/final_edges.csv"),
        ("data/nodes.csv",       "data/edges.csv"),
    ]
    for n, e in pairs:
        if os.path.exists(n) and os.path.exists(e):
            return n, e
    sys.exit("Error: no node/edge CSV files found. Run 02_build_final_pathway_files.py first.")


def _find_dot() -> str | None:
    for name in ("output/myrf_final_pathway.dot", "output/myrf_pathway.dot"):
        if os.path.exists(name):
            return name
    return None


# ── KGML subtype → value mapping (full KEGG list) ────────────────────────────

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


# ── KGML type / subtype resolution ───────────────────────────────────────────

def _kgml_type(level: str, row: "pd.Series | None" = None) -> str:
    """
    Return KEGG relation type.
    If row contains a populated 'kgml_type' column (written by script 06),
    use that directly.  Otherwise fall back to level-based inference.
    """
    if row is not None:
        col_val = _s(row.get("kgml_type", ""))
        if col_val in ("GErel", "PPrel", "PCrel", "maplink"):
            return col_val
    return "PPrel" if level == "protein" else "GErel"


def _kgml_subtype(interaction: str, row: "pd.Series | None" = None) -> tuple[str, str]:
    """
    Return (subtype_name, subtype_value).
    If row contains a populated 'kegg_edge_subtype' column (written by script 06),
    use that and look up the value from SUBTYPE_VALUES.
    Otherwise infer from the interaction string.
    """
    if row is not None:
        subtype = _s(row.get("kegg_edge_subtype", ""))
        if subtype in SUBTYPE_VALUES:
            return subtype, SUBTYPE_VALUES[subtype]

    ix = interaction.lower()
    if any(k in ix for k in ("activat", "expression", "upregulat", "gene_to_protein")):
        return "expression", SUBTYPE_VALUES["expression"]
    if any(k in ix for k in ("inhibit", "repress", "degradat", "downregulat")):
        return "inhibition", SUBTYPE_VALUES["inhibition"]
    if any(k in ix for k in ("bind", "associat", "coactiv")):
        return "binding/association", SUBTYPE_VALUES["binding/association"]
    if "phospho" in ix:
        return "phosphorylation", SUBTYPE_VALUES["phosphorylation"]
    if "ubiquit" in ix:
        return "ubiquitination", SUBTYPE_VALUES["ubiquitination"]
    if "indirect" in ix:
        return "indirect effect", SUBTYPE_VALUES["indirect effect"]
    if "repress" in ix:
        return "repression", SUBTYPE_VALUES["repression"]
    # Default
    return "expression", SUBTYPE_VALUES["expression"]


# ── Layout: parse dot -Tplain ─────────────────────────────────────────────────

def _dot_plain_positions(dot_file: str) -> dict[str, tuple[float, float]] | None:
    """
    Run  dot -Tplain  and extract node centres in inches.
    Returns {name: (x_in, y_in)} with graphviz bottom-left origin, or None.
    """
    try:
        r = subprocess.run(
            ["dot", "-Tplain", dot_file],
            capture_output=True, text=True, timeout=30,
        )
        if r.returncode != 0:
            return None
        pos = {}
        page_h = None
        for line in r.stdout.splitlines():
            parts = line.split()
            if not parts:
                continue
            if parts[0] == "graph" and len(parts) >= 4:
                page_h = float(parts[3])   # page height in inches
            elif parts[0] == "node" and len(parts) >= 4:
                name = parts[1].strip('"')
                pos[name] = (float(parts[2]), float(parts[3]))
        if pos and page_h is not None:
            # Flip Y: graphviz origin is bottom-left; KGML is top-left
            pos = {n: (x, page_h - y) for n, (x, y) in pos.items()}
        return pos or None
    except Exception:
        return None


def _scale_to_canvas(
    raw: dict[str, tuple[float, float]],
    canvas_w: int = KGML_W,
    canvas_h: int = KGML_H,
    margin: int = MARGIN,
) -> dict[str, tuple[int, int]]:
    """Scale raw inch positions (already Y-flipped) into KGML pixel canvas."""
    xs = [v[0] for v in raw.values()]
    ys = [v[1] for v in raw.values()]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    span_x = max_x - min_x or 1.0
    span_y = max_y - min_y or 1.0
    rx = (canvas_w - 2 * margin) / span_x
    ry = (canvas_h - 2 * margin) / span_y
    return {
        name: (
            int(margin + (x - min_x) * rx),
            int(margin + (y - min_y) * ry),
        )
        for name, (x, y) in raw.items()
    }


# ── Layout: algorithmic fallback ─────────────────────────────────────────────

def _spread_y(names: list[str], x: int, canvas_h: int, margin: int) -> dict[str, tuple[int, int]]:
    n = len(names)
    if n == 0:
        return {}
    step = (canvas_h - 2 * margin) / n
    return {name: (x, int(margin + (i + 0.5) * step)) for i, name in enumerate(names)}


def _algorithmic_layout(
    all_ids: set[str],
    upstream: set[str],
    downstream: set[str],
) -> dict[str, tuple[int, int]]:
    pos: dict[str, tuple[int, int]] = {}
    cx = KGML_W // 2

    # Centre
    pos["MYRF_gene"]    = (cx, KGML_H // 2 - MYRF_SPLIT_OFFSET)
    pos["MYRF_protein"] = (cx, KGML_H // 2 + MYRF_SPLIT_OFFSET)

    # Upstream – left column
    up_x = MARGIN + NODE_W // 2
    pos.update(_spread_y(sorted(upstream), up_x, KGML_H, MARGIN))

    # Downstream – one or two right columns
    dn = sorted(downstream)
    if len(dn) > 12:
        mid = len(dn) // 2
        pos.update(_spread_y(dn[:mid],       KGML_W - MARGIN - NODE_W - 10, KGML_H, MARGIN))
        pos.update(_spread_y(dn[mid:],       KGML_W - MARGIN + 10,          KGML_H, MARGIN))
    else:
        pos.update(_spread_y(dn, KGML_W - MARGIN - NODE_W // 2, KGML_H, MARGIN))

    # Anything else – bottom row
    others = sorted(
        n for n in all_ids
        if n not in pos and n not in {"MYRF_gene", "MYRF_protein"}
    )
    for i, name in enumerate(others):
        pos[name] = (MARGIN + i * (NODE_W + 10), KGML_H - MARGIN)

    return pos


# ── MYRF split helpers ────────────────────────────────────────────────────────

def _route_myrf_target(level: str) -> str:
    """Return 'MYRF_protein' for protein-level edges, else 'MYRF_gene'."""
    return "MYRF_protein" if level == "protein" else "MYRF_gene"


def _split_myrf(
    nodes_df: pd.DataFrame,
    edges_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Replace the single MYRF node with MYRF_gene + MYRF_protein entries.
    Route incoming edges by regulation_level; outgoing edges from MYRF → MYRF_protein.
    Append a synthetic MYRF_gene → MYRF_protein expression edge.

    If MYRF_gene and MYRF_protein already exist in the node table (because
    script 06 already performed the split), this function returns the data
    unchanged.
    """
    existing_ids = set(nodes_df["node_id"])
    if "MYRF_gene" in existing_ids and "MYRF_protein" in existing_ids:
        # Script 06 already did the split — nothing to do
        return nodes_df, edges_df

    myrf_rows = nodes_df[nodes_df["node_id"] == "MYRF"]
    if myrf_rows.empty:
        # Already split or absent — check and return as-is
        return nodes_df, edges_df

    base = myrf_rows.iloc[0].to_dict()

    gene_row = {**base,
                "node_id":        "MYRF_gene",
                "display_name":   "MYRF (gene)",
                "regulation_level": "transcription",
                "pathway_role":   "MYRF transcriptional locus"}
    prot_row = {**base,
                "node_id":        "MYRF_protein",
                "display_name":   "MYRF (protein)",
                "regulation_level": "protein",
                "pathway_role":   "MYRF N-terminal TF (autoproteolytic cleavage product)"}

    nodes_df = pd.concat(
        [nodes_df[nodes_df["node_id"] != "MYRF"],
         pd.DataFrame([gene_row, prot_row])],
        ignore_index=True,
    )

    def _reroute(row):
        r = row.copy()
        if r["target"] == "MYRF":
            # Support both old 'level' column and new 'regulation_level' column
            lvl = _s(r.get("regulation_level", "")) or _s(r.get("level", ""))
            r["target"] = _route_myrf_target(lvl)
        if r["source"] == "MYRF":
            r["source"] = "MYRF_protein"
        return r

    edges_df = edges_df.apply(_reroute, axis=1)

    synth = {
        "source":            "MYRF_gene",
        "source_ensembl_id": _s(base.get("ensembl_id", "")),
        "target":            "MYRF_protein",
        "target_ensembl_id": _s(base.get("ensembl_id", "")),
        "interaction":       "gene_to_protein_expression",
        "level":             "gene",
        "evidence":          "",   # no URL — this is a known biological fact, not a curated edge
        "remarks":           ("Biological fact: MYRF undergoes autoproteolytic cleavage; "
                              "N-terminal fragment translocates to nucleus as active TF."),
    }
    edges_df = pd.concat([edges_df, pd.DataFrame([synth])], ignore_index=True)
    return nodes_df, edges_df


# ── Validation ────────────────────────────────────────────────────────────────

def _validate(nodes_df: pd.DataFrame, edges_df: pd.DataFrame) -> list[str]:
    errors = []
    valid_ids = set(nodes_df["node_id"])

    # 1-2: every edge endpoint exists
    for col in ("source", "target"):
        missing = set(edges_df[~edges_df[col].isin(valid_ids)][col])
        for m in missing:
            errors.append(f"Edge {col} '{m}' not found in node table")

    # 3: no duplicate node_id
    dups = nodes_df[nodes_df.duplicated("node_id", keep=False)]["node_id"].unique()
    for d in dups:
        errors.append(f"Duplicate node_id: {d}")

    # 6: MYRF_gene and MYRF_protein present
    for req in ("MYRF_gene", "MYRF_protein"):
        if req not in valid_ids:
            errors.append(f"Required node '{req}' missing after split")

    # 7: MYRF_gene → MYRF_protein edge exists
    has_link = ((edges_df["source"] == "MYRF_gene") & (edges_df["target"] == "MYRF_protein")).any()
    if not has_link:
        errors.append("No MYRF_gene → MYRF_protein edge")

    # 8: at least one upstream edge into MYRF_gene
    if not (edges_df["target"] == "MYRF_gene").any():
        errors.append("No upstream edges targeting MYRF_gene")

    # 9: at least one downstream edge from MYRF_protein
    if not (edges_df["source"] == "MYRF_protein").any():
        errors.append("No downstream edges sourced from MYRF_protein")

    return errors


# ── XML pretty-printer ────────────────────────────────────────────────────────

def _prettify(element: ET.Element) -> str:
    rough = ET.tostring(element, encoding="unicode")
    dom = minidom.parseString(rough)
    pretty = dom.toprettyxml(indent="  ")
    # Remove blank lines that toprettyxml inserts
    lines = [ln for ln in pretty.splitlines() if ln.strip()]
    return "\n".join(lines) + "\n"


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> int:
    os.makedirs(OUT_DIR, exist_ok=True)

    nodes_path, edges_path = _pick_files()
    nodes_df = pd.read_csv(nodes_path)
    edges_df = pd.read_csv(edges_path)
    print(f"Loaded {len(nodes_df)} nodes  ({nodes_path})")
    print(f"Loaded {len(edges_df)} edges  ({edges_path})")

    # ── MYRF split ────────────────────────────────────────────────────────────
    nodes_df, edges_df = _split_myrf(nodes_df, edges_df)
    print(f"After MYRF split: {len(nodes_df)} nodes, {len(edges_df)} edges")

    # ── Validate ──────────────────────────────────────────────────────────────
    errors = _validate(nodes_df, edges_df)
    if errors:
        print("\nValidation FAILED:")
        for e in errors:
            print(f"  ERROR: {e}")
        return 1
    print("Validation passed (all 9 checks).")

    # ── Entry ID assignment ───────────────────────────────────────────────────
    # MYRF_gene=1, MYRF_protein=2, rest alphabetical
    priority = ["MYRF_gene", "MYRF_protein"]
    rest     = sorted(n for n in nodes_df["node_id"] if n not in priority)
    ordered  = priority + rest
    entry_id = {nid: i + 1 for i, nid in enumerate(ordered)}

    # 4-5: duplicate entry IDs impossible by construction; check 5 at relation time
    assert len(entry_id) == len(ordered), "BUG: duplicate entry ID"

    # ── Layout ────────────────────────────────────────────────────────────────
    dot_file = _find_dot()
    positions: dict[str, tuple[int, int]] = {}

    if dot_file:
        raw = _dot_plain_positions(dot_file)
        if raw:
            scaled = _scale_to_canvas(raw)
            # DOT file has single MYRF; split into gene/protein offset vertically
            if "MYRF" in scaled:
                mx, my = scaled.pop("MYRF")
                scaled["MYRF_gene"]    = (mx, max(MARGIN, my - MYRF_SPLIT_OFFSET))
                scaled["MYRF_protein"] = (mx, min(KGML_H - MARGIN, my + MYRF_SPLIT_OFFSET))
            positions = scaled
            print(f"Layout: scaled from DOT positions ({dot_file})")

    if not positions:
        valid_ids = set(nodes_df["node_id"])
        upstream_set   = (set(edges_df[edges_df["target"] == "MYRF_gene"]["source"])
                          | set(edges_df[edges_df["target"] == "MYRF_protein"]["source"])
                          ) - {"MYRF_gene", "MYRF_protein"}
        downstream_set = (set(edges_df[edges_df["source"] == "MYRF_protein"]["target"])
                          | set(edges_df[edges_df["source"] == "MYRF_gene"]["target"])
                          ) - {"MYRF_gene", "MYRF_protein"}
        positions = _algorithmic_layout(valid_ids, upstream_set, downstream_set)
        print("Layout: algorithmic (dot -Tplain unavailable)")

    # Fill fallback for any node still missing a position
    for nid in nodes_df["node_id"]:
        if nid not in positions:
            positions[nid] = (KGML_W // 2, KGML_H // 2)

    # ── Build KGML XML ────────────────────────────────────────────────────────
    pathway = ET.Element("pathway")
    pathway.set("name",   "path:custom_myrf0001")
    pathway.set("org",    "hsa")
    pathway.set("number", "00001")
    pathway.set("title",  "Custom MYRF-centered oligodendrocyte differentiation pathway")

    # Entries
    node_lookup = nodes_df.set_index("node_id")

    for nid in ordered:
        row  = node_lookup.loc[nid]
        eid  = entry_id[nid]
        x, y = positions[nid]

        entrez       = ENTREZ_MAP.get(nid)
        kgml_name    = f"hsa:{entrez}" if entrez else f"custom:{nid}"
        display_name = _s(row.get("display_name")) or nid

        entry = ET.SubElement(pathway, "entry")
        entry.set("id",   str(eid))
        entry.set("name", kgml_name)
        entry.set("type", "gene")

        g = ET.SubElement(entry, "graphics")
        g.set("name",   display_name)
        g.set("type",   "rectangle")
        g.set("x",      str(x))
        g.set("y",      str(y))
        g.set("width",  str(NODE_W))
        g.set("height", str(NODE_H))

        ensembl = _s(row.get("ensembl_id", ""))
        if ensembl:
            g.set("ensembl_id", ensembl)

    # Relations + evidence rows
    evidence_rows = []
    valid_ids = set(nodes_df["node_id"])
    n_expr = n_inh = n_bind = n_other = 0
    skipped = []

    for _, row in edges_df.iterrows():
        src = row["source"]
        tgt = row["target"]

        # 5: skip edges whose endpoints are outside the entry map (log them)
        if src not in entry_id or tgt not in entry_id:
            skipped.append(f"{src} -> {tgt}")
            continue

        e1          = entry_id[src]
        e2          = entry_id[tgt]
        # Support both old schema (interaction/level/remarks/evidence) and
        # new schema (primary_relationship / regulation_level / notes / source_urls)
        interaction = _s(row.get("interaction", "")) or _s(row.get("primary_relationship", ""))
        level       = _s(row.get("regulation_level", "")) or _s(row.get("level", ""))
        remarks     = _s(row.get("notes", "")) or _s(row.get("remarks", ""))
        evidence    = _s(row.get("source_urls", "")) or _s(row.get("evidence", ""))
        direction_raw = _s(row.get("direction", ""))

        ktype            = _kgml_type(level, row)
        ksubtype, kvalue = _kgml_subtype(interaction, row)
        if direction_raw:
            direction = direction_raw
        else:
            direction = ("positive" if "activat" in interaction.lower()
                         else "negative" if "inhibit" in interaction.lower()
                         else "unknown")

        # Count
        if ksubtype in ("expression", "activation"):
            n_expr += 1
        elif ksubtype in ("inhibition", "repression", "ubiquitination"):
            n_inh += 1
        elif ksubtype == "binding/association":
            n_bind += 1
        else:
            n_other += 1
        rel = ET.SubElement(pathway, "relation")
        rel.set("entry1", str(e1))
        rel.set("entry2", str(e2))
        rel.set("type",   ktype)

        sub = ET.SubElement(rel, "subtype")
        sub.set("name",  ksubtype)
        sub.set("value", kvalue)

        # Evidence row
        src_node = node_lookup.loc[src] if src in node_lookup.index else {}
        tgt_node = node_lookup.loc[tgt] if tgt in node_lookup.index else {}

        # For outgoing MYRF edges, source_bucket "Central node" is uninformative;
        # use the TARGET node's source_bucket to show where the relationship was sourced.
        is_myrf_src = src in ("MYRF_gene", "MYRF_protein")
        bucket_node = tgt_node if is_myrf_src else src_node

        # Synthetic edge detected by unique interaction name (not by remarks text)
        is_synthetic = interaction == "gene_to_protein_expression"

        # Confidence: prefer row-level column, then fall back to node metadata
        row_confidence = _s(row.get("confidence", ""))
        node_confidence = _s(src_node.get("evidence_strength", "")) if hasattr(src_node, "get") else ""
        confidence_val = row_confidence or node_confidence

        # evidence_type: prefer row-level column then infer
        row_ev_type = _s(row.get("evidence_type", ""))
        if is_synthetic:
            ev_type_val = "synthetic"
        elif row_ev_type:
            ev_type_val = row_ev_type
        else:
            ev_type_val = "alias_note" if remarks else ""

        evidence_rows.append({
            "source":          src,
            "target":          tgt,
            "source_entry_id": e1,
            "target_entry_id": e2,
            "interaction":     interaction,
            "kgml_type":       ktype,
            "kgml_subtype":    ksubtype,
            "direction":       direction,
            "level":           level,
            "confidence":      confidence_val,
            "evidence_type":   ev_type_val,
            "source_urls":     evidence,
            "notes":           remarks,
            "needs_manual_review": _s(row.get("needs_manual_review", "")),
            "section":         _s(bucket_node.get("source_bucket", "")) if hasattr(bucket_node, "get") else "",
        })

    if skipped:
        print(f"  Warning: {len(skipped)} edges skipped (endpoint not in node table):")
        for s in skipped:
            print(f"    {s}")

    # ── Write KGML ───────────────────────────────────────────────────────────
    header = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<!-- Custom MYRF-centred pathway (not an official KEGG pathway). '
        'Generated by 04_export_kgml.py. -->\n'
    )
    xml_body = _prettify(pathway)
    # _prettify starts with the minidom XML declaration; strip it
    if xml_body.startswith("<?xml"):
        xml_body = xml_body[xml_body.index("\n") + 1:]

    with open(OUT_KGML, "w", encoding="utf-8") as f:
        f.write(header)
        f.write(xml_body)

    os.makedirs(os.path.dirname(OUT_CUSTOM_KGML), exist_ok=True)
    shutil.copyfile(OUT_KGML, OUT_CUSTOM_KGML)
    # ── Write evidence CSV ───────────────────────────────────────────────────
    ev_df = pd.DataFrame(evidence_rows, columns=[
        "source", "target", "source_entry_id", "target_entry_id",
        "interaction", "kgml_type", "kgml_subtype",
        "direction", "level", "confidence", "evidence_type",
        "source_urls", "notes", "needs_manual_review", "section",
    ])
    # QUOTE_NONNUMERIC wraps every string field in double-quotes so that URLs
    # and adjacent text columns (notes, section) have unambiguous boundaries
    # in any plain-text viewer or editor.
    ev_df.to_csv(OUT_EV, index=False, encoding="utf-8", quoting=csv.QUOTE_NONNUMERIC)

    # ── Write entry ID map ───────────────────────────────────────────────────
    map_rows = []
    for nid in ordered:
        row = node_lookup.loc[nid]
        x, y = positions[nid]
        map_rows.append({
            "entry_id":    entry_id[nid],
            "node_id":     nid,
            "display_name": _s(row.get("display_name")) or nid,
            "kgml_name":   f"hsa:{ENTREZ_MAP[nid]}" if nid in ENTREZ_MAP else f"custom:{nid}",
            "kgml_type":   "gene",
            "x":           x,
            "y":           y,
            "width":       NODE_W,
            "height":      NODE_H,
            "role":        _s(row.get("pathway_position", "")),
            "confidence":  _s(row.get("evidence_strength", "")),
        })
    pd.DataFrame(map_rows).to_csv(OUT_MAP, index=False, encoding="utf-8")

    # ── Summary ───────────────────────────────────────────────────────────────
    n_candidate = ev_df[
        ev_df["confidence"].str.lower().isin(["moderate", "weak", "candidate"])
    ].shape[0]

    print(f"\n{'='*62}")
    print("SUMMARY")
    print(f"{'='*62}")
    print(f"  Nodes loaded             : {len(nodes_df)}")
    print(f"  Edges loaded             : {len(edges_df)}")
    print(f"  KGML entries written     : {len(ordered)}")
    print(f"  KGML relations written   : {len(evidence_rows)}")
    print(f"    expression/activation  : {n_expr}")
    print(f"    inhibition/repression  : {n_inh}")
    print(f"    binding/association    : {n_bind}")
    print(f"    other                  : {n_other}")
    print(f"    candidate/weak/moderate: {n_candidate}")
    print(f"\n  Outputs:")
    for path in (OUT_KGML, OUT_EV, OUT_MAP):
        tag = "[OK]" if os.path.exists(path) else "[MISSING]"
        print(f"    {tag} {os.path.abspath(path)}")

    print(
        "\nWARNING: This is a custom KGML-like pathway and may not be accepted\n"
        "directly by all KEGG/Pathview workflows without adaptation.\n"
        "(Gene names use  custom:SYMBOL  ->  substitute  hsa:ENTREZID  for\n"
        " production use with KEGG/Pathview.)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
