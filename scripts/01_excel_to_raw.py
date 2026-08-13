"""
06_excel_to_raw.py  —  Step 1 of 4
Read the curated Excel master sheet and write two verifiable intermediate
CSV files.  No KEGG classification is applied here.

Input:
    raw_data/myrf-final-consolidation-documentation.xlsx
        Sheet: Master_Gene_List
        Sheet: References

Outputs  (verify these before running script 07):
    data/raw_edges.csv          — one row per pathway edge, MYRF-split applied
    data/raw_master_gene_list.csv — normalised Master_Gene_List for node metadata

Intermediate checks to make:
  - raw_edges.csv:  correct genes? right source/target per row?
                    protein-level edges routed to MYRF_protein?
                    transcription-level edges routed to MYRF_gene?
                    synthetic MYRF_gene->MYRF_protein row present?
  - raw_master_gene_list.csv: all 57 genes, all columns?
"""

import os
import sys

import csv
import pandas as pd

# ── Paths ─────────────────────────────────────────────────────────────────────

EXCEL_PATH     = os.path.join("raw_data", "myrf-final-consolidation-documentation.xlsx")
DATA_DIR       = "data"
OUT_RAW_EDGES  = os.path.join(DATA_DIR, "raw_edges.csv")
OUT_RAW_MASTER = os.path.join(DATA_DIR, "raw_master_gene_list.csv")


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


# ── Load & normalise Excel ────────────────────────────────────────────────────

def load_excel(path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    xl = pd.ExcelFile(path, engine="openpyxl")
    master = xl.parse("Master_Gene_List")
    refs   = xl.parse("References")
    return master, refs


def normalise_master(df: pd.DataFrame) -> pd.DataFrame:
    """Rename columns by position to canonical internal names."""
    col_map = {
        df.columns[0]:  "Gene_Symbol",
        df.columns[1]:  "Aliases_or_Notes",
        df.columns[2]:  "Original_Label",
        df.columns[3]:  "Ensembl_ID",
        df.columns[4]:  "Source_Bucket",
        df.columns[5]:  "Pathway_Position",
        df.columns[6]:  "Primary_Relationship_to_MYRF",
        df.columns[7]:  "Direction",
        df.columns[8]:  "Regulation_Level",
        df.columns[9]:  "All_Relationships",
        df.columns[10]: "Evidence_Strength",
        df.columns[11]: "Pathway_Role",
        df.columns[12]: "Notes",
        df.columns[13]: "References",
        df.columns[14]: "Reference_Count",
    }
    return df.rename(columns=col_map)


# ── Build raw edges from Master_Gene_List rows ────────────────────────────────

def build_raw_edges(master: pd.DataFrame) -> pd.DataFrame:
    """
    One edge per Master_Gene_List row that has a Pathway_Position.
      Upstream   → source=gene,  target=MYRF
      Downstream → source=MYRF,  target=gene
    """
    rows = []
    for _, r in master.iterrows():
        sym = _s(r["Gene_Symbol"])
        pos = _s(r["Pathway_Position"])
        if not pos or pos == "Center":
            continue

        if pos == "Upstream":
            src, tgt = sym, "MYRF"
        elif pos == "Downstream":
            src, tgt = "MYRF", sym
        else:
            continue

        rows.append({
            "source":               src,
            "target":               tgt,
            "gene_symbol":          sym,
            "pathway_position":     pos,
            "primary_relationship": _s(r["Primary_Relationship_to_MYRF"]),
            "direction":            _s(r["Direction"]),
            "regulation_level":     _s(r["Regulation_Level"]),
            "confidence":           _s(r["Evidence_Strength"]),
            "source_bucket":        _s(r["Source_Bucket"]),
            "source_urls":          _s(r["References"]),
            "notes":                _s(r["Notes"]),
        })
    return pd.DataFrame(rows)


# ── Deduplication ─────────────────────────────────────────────────────────────

def _merge_unique(values: pd.Series, separator: str = " ; ") -> str:
    """Merge delimited evidence without discarding later workbook rows."""
    merged: list[str] = []
    for value in values:
        for item in _s(value).split(";"):
            item = item.strip()
            if item and item not in merged:
                merged.append(item)
    return separator.join(merged)


def _strongest(values: pd.Series) -> str:
    rank = {"excluded": -1, "weak": 0, "moderate": 1, "strong": 2}
    cleaned = [_s(v) for v in values if _s(v)]
    return max(cleaned, key=lambda v: rank.get(v.lower(), 0), default="")


def deduplicate_edges(edges: pd.DataFrame) -> pd.DataFrame:
    """Merge duplicate rows while retaining every citation and note."""
    keys = [
        "source", "target", "regulation_level",
        "primary_relationship", "direction",
    ]
    aggregated = (
        edges.groupby(keys, as_index=False, sort=False, dropna=False)
        .agg({
            "gene_symbol": "first",
            "pathway_position": "first",
            "confidence": _strongest,
            "source_bucket": _merge_unique,
            "source_urls": _merge_unique,
            "notes": _merge_unique,
        })
    )
    return aggregated[edges.columns].reset_index(drop=True)

# ── MYRF split ────────────────────────────────────────────────────────────────

def apply_myrf_split(edges: pd.DataFrame) -> pd.DataFrame:
    """
    Route MYRF endpoints to MYRF_gene or MYRF_protein based on regulation_level:
      Incoming (target=MYRF):
        protein      → MYRF_protein
        transcription or unknown → MYRF_gene
      Outgoing (source=MYRF):
        all → MYRF_protein  (the active cleaved TF is the source of transcription)

    Then append the biological-fact synthetic edge MYRF_gene → MYRF_protein.
    """
    edges = edges.copy()

    def _route_target(row) -> str:
        if row["target"] != "MYRF":
            return row["target"]
        return "MYRF_protein" if _lower(row["regulation_level"]) == "protein" else "MYRF_gene"

    def _route_source(row) -> str:
        if row["source"] != "MYRF":
            return row["source"]
        return "MYRF_protein"

    edges["target"] = edges.apply(_route_target, axis=1)
    edges["source"] = edges.apply(_route_source, axis=1)

    synth = {
        "source":               "MYRF_gene",
        "target":               "MYRF_protein",
        "gene_symbol":          "MYRF",
        "pathway_position":     "Center",
        "primary_relationship": "gene_to_protein_expression",
        "direction":            "positive",
        "regulation_level":     "gene",
        "confidence":           "Strong",
        "source_bucket":        "Primary literature",
        "source_urls":          "https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001625",
        "notes": (
            "MYRF is translated as an ER-membrane precursor and undergoes "
            "autoproteolytic cleavage; the N-terminal trimer enters the nucleus."
        ),
    }
    return pd.concat([edges, pd.DataFrame([synth])], ignore_index=True)


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> int:
    os.makedirs(DATA_DIR, exist_ok=True)

    if not os.path.exists(EXCEL_PATH):
        sys.exit(f"ERROR: Excel not found: {EXCEL_PATH}")

    print(f"Reading: {EXCEL_PATH}")
    master_raw, refs_df = load_excel(EXCEL_PATH)
    master = normalise_master(master_raw)
    print(f"  Master_Gene_List: {len(master)} rows x {len(master.columns)} cols")

    # Write normalised master as-is for script 07 to use for node metadata
    master.to_csv(OUT_RAW_MASTER, index=False, encoding="utf-8", quoting=csv.QUOTE_NONNUMERIC)
    print(f"  Written: {OUT_RAW_MASTER}")

    # Build and process edges
    raw   = build_raw_edges(master)
    print(f"\n  Edges from MGL rows      : {len(raw)}")

    dedup = deduplicate_edges(raw)
    dropped = len(raw) - len(dedup)
    print(f"  After pre-classification dedup : {len(dedup)}"
          + (f"  (-{dropped} exact duplicates)" if dropped else ""))

    split = apply_myrf_split(dedup)
    print(f"  After MYRF split + synthetic   : {len(split)}")

    # Print MYRF routing breakdown
    to_gene  = split[(split["target"] == "MYRF_gene")  & (split["source"] != "MYRF_gene")]
    to_prot  = split[(split["target"] == "MYRF_protein") & (split["source"] != "MYRF_gene")]
    from_prot = split[split["source"] == "MYRF_protein"]
    synth_row = split[(split["source"] == "MYRF_gene") & (split["target"] == "MYRF_protein")]

    print(f"\n  MYRF routing:")
    print(f"    -> MYRF_gene    (transcription/unknown level) : {len(to_gene)} edges")
    print(f"       sources: {', '.join(sorted(to_gene['source'].unique()))}")
    print(f"    -> MYRF_protein (protein level)              : {len(to_prot)} edges")
    print(f"       sources: {', '.join(sorted(to_prot['source'].unique()))}")
    print(f"    from MYRF_protein -> downstream              : {len(from_prot)} edges")
    print(f"    MYRF_gene -> MYRF_protein (synthetic)        : {len(synth_row)}")

    split.to_csv(OUT_RAW_EDGES, index=False, encoding="utf-8", quoting=csv.QUOTE_NONNUMERIC)
    print(f"\n[OK] {OUT_RAW_EDGES}")
    print(f"[OK] {OUT_RAW_MASTER}")
    print("\nVerify these files before running scripts/02_build_pathway_csv.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
