"""
08_preview_kegg_graph.py  —  Step 3 of 4
Generate a KEGG-style DOT / SVG preview of the curated pathway.

Inputs  (from script 07):
    data/final_nodes.csv
    data/final_edges.csv

Outputs:
    output/myrf_final_pathway.dot
    output/myrf_final_pathway.svg

Visual encoding
───────────────
Nodes
  Shape  : gene    → box (rounded)
           group   → box (bold border)
           compound→ ellipse
           map     → diamond
  Fill   : core_hub            → #FFD700  gold
           upstream_regulator  → #AED6F1  blue
           downstream_target   → #A9DFBF  green
           candidate_context   → #D5DBDB  grey

Edges (by kegg_edge_subtype)
  expression / activation : green  solid  normal-arrow       penwidth by confidence
  repression / inhibition : red    solid  tee-arrow
  binding/association     : blue   solid  no-arrow (dir=none)
  indirect effect         : grey   dashed normal-arrow
  state change            : purple dotted normal-arrow
  ubiquitination          : orange solid  normal-arrow  label=+u
  phosphorylation         : orange solid  normal-arrow  label=+p
  dephosphorylation       : teal   solid  normal-arrow  label=-p
  others/unknown          : grey   dashed normal-arrow  label=?

Confidence → penwidth
  Strong   → 2.0
  Moderate → 1.2
  Weak     → 0.8  (line also dashed if not already)

needs_manual_review = True → label suffix " [review]"
"""

import os
import subprocess
import sys

import pandas as pd
import pydot

# ── Paths ─────────────────────────────────────────────────────────────────────

DATA_DIR  = "data"
OUT_DIR   = "output"
IN_NODES  = os.path.join(DATA_DIR, "final_nodes.csv")
IN_EDGES  = os.path.join(DATA_DIR, "final_edges.csv")
OUT_DOT   = os.path.join(OUT_DIR,  "myrf_final_pathway.dot")
OUT_SVG   = os.path.join(OUT_DIR,  "myrf_final_pathway.svg")

# ── Visual constants ──────────────────────────────────────────────────────────

NODE_FILL = {
    "core_hub":            "#FFD700",
    "upstream_regulator":  "#AED6F1",
    "downstream_target":   "#A9DFBF",
    "candidate_context":   "#D5DBDB",
}
NODE_SHAPE = {
    "gene":     "box",
    "group":    "box",
    "compound": "ellipse",
    "map":      "diamond",
}

EDGE_STYLE: dict[str, dict] = {
    "expression":          {"color": "#27AE60", "style": "solid",  "arrowhead": "normal", "label": ""},
    "activation":          {"color": "#27AE60", "style": "solid",  "arrowhead": "normal", "label": ""},
    "repression":          {"color": "#C0392B", "style": "solid",  "arrowhead": "tee",    "label": ""},
    "inhibition":          {"color": "#C0392B", "style": "solid",  "arrowhead": "tee",    "label": ""},
    "binding/association": {"color": "#2980B9", "style": "solid",  "arrowhead": "none",   "label": ""},
    "indirect effect":     {"color": "#7F8C8D", "style": "dashed", "arrowhead": "normal", "label": ""},
    "state change":        {"color": "#8E44AD", "style": "dotted", "arrowhead": "normal", "label": ""},
    "ubiquitination":      {"color": "#E67E22", "style": "solid",  "arrowhead": "normal", "label": "+u"},
    "phosphorylation":     {"color": "#E67E22", "style": "solid",  "arrowhead": "normal", "label": "+p"},
    "dephosphorylation":   {"color": "#16A085", "style": "solid",  "arrowhead": "normal", "label": "-p"},
    "glycosylation":       {"color": "#E67E22", "style": "solid",  "arrowhead": "normal", "label": "+g"},
    "methylation":         {"color": "#8E44AD", "style": "solid",  "arrowhead": "normal", "label": "+m"},
    "others/unknown":      {"color": "#95A5A6", "style": "dashed", "arrowhead": "normal", "label": "?"},
}
DEFAULT_EDGE_STYLE = {"color": "#95A5A6", "style": "dashed", "arrowhead": "normal", "label": "?"}

CONF_PENWIDTH = {"strong": "2.0", "moderate": "1.2", "weak": "0.8"}


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


_UNICODE_SUBS = str.maketrans({
    "—": " - ", "–": " - ", "→": "->",
    "←": "<-",  "’": "'",   "‘": "'",
    "“": '"',   "”": '"',
})

def _dot_safe(text: str) -> str:
    """Replace Unicode chars that break Graphviz SVG on Windows."""
    text = text.translate(_UNICODE_SUBS)
    return text.encode("ascii", errors="ignore").decode("ascii")


def _quoted(s: str) -> str:
    return f'"{_dot_safe(s)}"'


# ── Layout helpers ────────────────────────────────────────────────────────────

def _compute_groups(edges_df: pd.DataFrame, valid_ids: set[str]) -> dict[str, str]:
    """Derive upstream/center/downstream from edge topology."""
    to_myrf   = set(edges_df.loc[edges_df["target"].isin({"MYRF_gene", "MYRF_protein"}), "source"])
    from_myrf = set(edges_df.loc[edges_df["source"].isin({"MYRF_gene", "MYRF_protein"}), "target"])
    groups: dict[str, str] = {}
    for nid in valid_ids:
        if nid in ("MYRF_gene", "MYRF_protein"):
            groups[nid] = "core_hub"
        elif nid in to_myrf:
            groups[nid] = "upstream_regulator"
        elif nid in from_myrf:
            groups[nid] = "downstream_target"
        else:
            groups[nid] = "candidate_context"
    return groups


# ── Legend builder ────────────────────────────────────────────────────────────

def _build_legend(present_subtypes: set[str]) -> str:
    """Return an HTML-like label string for the legend node."""
    rows = [
        '<<TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="4">',
        '<TR><TD COLSPAN="2" ALIGN="CENTER"><B>Legend</B></TD></TR>',
        '<TR><TD COLSPAN="2" ALIGN="CENTER"><I>Node fills</I></TD></TR>',
        '<TR><TD BGCOLOR="#FFD700">MYRF (core hub)</TD>'
        '<TD BGCOLOR="#AED6F1">Upstream regulator</TD></TR>',
        '<TR><TD BGCOLOR="#A9DFBF">Downstream target</TD>'
        '<TD BGCOLOR="#D5DBDB">Candidate / context</TD></TR>',
        '<TR><TD COLSPAN="2" ALIGN="CENTER"><I>Edge subtypes present</I></TD></TR>',
    ]

    # Use &gt; for '>' to avoid breaking graphviz HTML label parser
    desc = {
        "expression":          "--&gt; green solid   : transcriptional activation",
        "repression":          "--| red solid        : transcriptional repression",
        "activation":          "--&gt; green solid   : protein activation",
        "inhibition":          "--| red solid        : protein inhibition",
        "binding/association": "--- blue solid       : protein cooperation/binding",
        "indirect effect":     "..&gt; dashed        : green positive; red tee negative",
        "ubiquitination":      "+u  orange tee       : ubiquitination/degradation",
        "phosphorylation":     "+p  orange solid     : phosphorylation",
        "dephosphorylation":   "-p  teal solid       : dephosphorylation",
        "others/unknown":      "?   grey dashed      : unknown / needs review",
    }
    for subtype in sorted(present_subtypes):
        text = desc.get(subtype, subtype)
        rows.append(f'<TR><TD COLSPAN="2" ALIGN="LEFT">{text}</TD></TR>')

    rows += [
        '<TR><TD COLSPAN="2" ALIGN="CENTER"><I>Line width = confidence</I></TD></TR>',
        '<TR><TD>Strong : thick (2.0)</TD><TD>Moderate : medium (1.2)</TD></TR>',
        '<TR><TD COLSPAN="2">Weak : thin (0.8) + dashed</TD></TR>',
        '<TR><TD COLSPAN="2">[ review ] = needs manual review</TD></TR>',
        '</TABLE>>',
    ]
    return "\n".join(rows)


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> int:
    os.makedirs(OUT_DIR, exist_ok=True)

    for path in (IN_NODES, IN_EDGES):
        if not os.path.exists(path):
            sys.exit(f"ERROR: {path} not found. Run script 07 first.")

    nodes_df = pd.read_csv(IN_NODES)
    edges_df = pd.read_csv(IN_EDGES)
    print(f"Loaded {len(nodes_df)} nodes from {IN_NODES}")
    print(f"Loaded {len(edges_df)} edges from {IN_EDGES}")

    # Use role from CSV; fall back to topology-derived group for colour
    valid_ids = set(nodes_df["node_id"])
    topo_groups = _compute_groups(edges_df, valid_ids)

    dot_g = pydot.Dot(
        graph_type="digraph",
        rankdir="LR",
        splines="true",
        nodesep="0.5",
        ranksep="1.8",
        fontname="Helvetica",
        bgcolor="white",
    )
    dot_g.set_node_defaults(fontname="Helvetica", fontsize="11", margin="0.15,0.08")
    dot_g.set_edge_defaults(fontname="Helvetica", fontsize="9")

    # ── Add nodes ─────────────────────────────────────────────────────────────
    pnodes: dict[str, pydot.Node] = {}
    for _, row in nodes_df.iterrows():
        nid      = row["node_id"]
        role     = _s(row.get("role")) or topo_groups.get(nid, "candidate_context")
        ktype    = _s(row.get("kegg_node_type")) or "gene"
        fill     = NODE_FILL.get(role, "#D5DBDB")
        shape    = NODE_SHAPE.get(ktype, "box")
        is_hub   = role == "core_hub"
        is_group = ktype == "group"

        # Tooltip: biological_entity_type + ensembl
        tip_parts = [
            p for p in (
                _dot_safe(_s(row.get("biological_entity_type"))),
                _dot_safe(_s(row.get("ensembl_id"))),
            ) if p
        ]
        tooltip = _quoted(" | ".join(tip_parts) if tip_parts else nid)

        # Style string
        style = "filled,rounded" if shape == "box" else "filled"
        if is_group:
            style = "filled"

        pn = pydot.Node(
            nid,
            label=_quoted(nid),
            shape=shape,
            style=_quoted(style),
            fillcolor=_quoted(fill),
            color=_quoted("#444444"),
            penwidth="3.0" if is_group else ("2.5" if is_hub else "1.0"),
            fontsize="13"  if is_hub else "11",
            fontname=_quoted("Helvetica-Bold" if is_hub else "Helvetica"),
            tooltip=tooltip,
        )
        pnodes[nid] = pn
        dot_g.add_node(pn)

    # ── Rank subgraphs ────────────────────────────────────────────────────────
    upstream   = [n for n, g in topo_groups.items() if g == "upstream_regulator"]
    center     = [n for n, g in topo_groups.items() if g == "core_hub"]
    downstream = [n for n, g in topo_groups.items() if g == "downstream_target"]

    for rank_val, members in [("source", upstream), ("same", center), ("sink", downstream)]:
        if members:
            sub = pydot.Subgraph(rank=rank_val)
            for nid in members:
                sub.add_node(pydot.Node(nid))
            dot_g.add_subgraph(sub)

    # ── Add edges ─────────────────────────────────────────────────────────────
    present_subtypes: set[str] = set()
    for _, row in edges_df.iterrows():
        src         = row["source"]
        tgt         = row["target"]
        subtype     = _s(row.get("kegg_edge_subtype")) or "others/unknown"
        conf        = _s(row.get("confidence")).lower()
        needs_rev   = str(row.get("needs_manual_review", "")).lower() in ("true", "1", "yes")
        ev_value    = _s(row.get("kegg_edge_value"))

        present_subtypes.add(subtype)
        es = dict(EDGE_STYLE.get(subtype, DEFAULT_EDGE_STYLE))
        dirn = _s(row.get("direction")).lower()

        # Preserve the sign of indirect and degradative mechanisms visually.
        if subtype == "indirect effect":
            if dirn == "negative":
                es.update(color="#C0392B", arrowhead="tee")
            elif dirn == "positive":
                es.update(color="#27AE60", arrowhead="normal")
        elif subtype == "ubiquitination" and dirn == "negative":
            es.update(arrowhead="tee")

        # Confidence → penwidth; Weak → also dash
        penwidth = CONF_PENWIDTH.get(conf, "1.2")
        style    = es["style"]
        if conf == "weak" and style == "solid":
            style = "dashed"

        # Edge label: subtype symbol + [review] suffix
        edge_label = es["label"]
        if needs_rev:
            edge_label = (edge_label + " [review]").strip()

        # dir attribute for binding/association (undirected look)
        extra: dict = {}
        if es["arrowhead"] == "none":
            extra["dir"] = "none"

        pe = pydot.Edge(
            src, tgt,
            label=_quoted(edge_label),
            color=_quoted(es["color"]),
            fontcolor=_quoted(es["color"]),
            style=_quoted(style),
            arrowhead=es["arrowhead"],
            arrowsize="0.9",
            penwidth=penwidth,
            **{k: _quoted(str(v)) for k, v in extra.items()},
        )
        dot_g.add_edge(pe)

    # ── Legend ────────────────────────────────────────────────────────────────
    legend_label = _build_legend(present_subtypes)
    legend_node  = pydot.Node(
        "legend_box",
        label=legend_label,
        shape="none",
        margin="0",
        fontname=_quoted("Helvetica"),
        fontsize="9",
    )
    dot_g.add_node(legend_node)

    # ── Write DOT ─────────────────────────────────────────────────────────────
    dot_g.write_raw(OUT_DOT)
    print(f"[OK] DOT written: {os.path.abspath(OUT_DOT)}")

    # ── Render SVG ────────────────────────────────────────────────────────────
    try:
        dot_g.write_svg(OUT_SVG)
        # Check for non-ASCII encoding issue
        with open(OUT_SVG, "rb") as f:
            bad = sum(1 for b in f.read() if b > 127)
        if bad:
            print(f"WARNING: SVG has {bad} non-ASCII bytes — may not render in browser")
        else:
            print(f"[OK] SVG written: {os.path.abspath(OUT_SVG)}")
    except Exception as e:
        print(f"WARNING: SVG render failed (Graphviz dot must be in PATH).\n  {e}")

    # ── Summary ───────────────────────────────────────────────────────────────
    up_n  = sum(1 for g in topo_groups.values() if g == "upstream_regulator")
    dn_n  = sum(1 for g in topo_groups.values() if g == "downstream_target")
    rev_n = int(edges_df.get("needs_manual_review", pd.Series(dtype=bool)).sum())

    print(f"\n{'='*50}")
    print("GRAPH SUMMARY")
    print(f"{'='*50}")
    print(f"  Nodes: {len(nodes_df)}")
    print(f"    upstream_regulator: {up_n}")
    print(f"    core_hub:           {len(center)}")
    print(f"    downstream_target:  {dn_n}")
    print(f"  Edges: {len(edges_df)}")
    for subtype, cnt in edges_df["kegg_edge_subtype"].value_counts().items():
        print(f"    {subtype}: {cnt}")
    print(f"  Needs manual review: {rev_n}")
    print(f"  DOT: {os.path.abspath(OUT_DOT)}")
    print(f"  SVG: {os.path.abspath(OUT_SVG)}")
    print(f"{'='*50}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
