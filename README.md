# MYRF-centered oligodendrocyte differentiation pathway

This repository builds a custom KEGG-style pathway centered on MYRF from a curated Excel workbook. It separates the MYRF gene locus from the active MYRF protein, classifies literature-supported regulatory relationships, and exports CSV, DOT/SVG, and KGML artifacts.

## Current audited model

The 2026-08-14 literature audit contains:

- 31 nodes
- 32 edges
- 23 strong and 9 moderate edges
- 0 unresolved manual-review flags

The authoritative visual is output/myrf_final_pathway.svg. The authoritative edge register is data/final_edges.csv, and its primary-literature cross-check is data/literature_audit.csv.

The full round-two result is in output/literature_audit_round2.md.

## Source of truth

raw_data/myrf-final-consolidation-documentation.xlsx is the curated source workbook.

- Master_Gene_List contains the source pathway rows.
- References is rebuilt from the citations used in the master sheet.
- Revision_Log records the first-round corrections and their rationale.
- Literature_Audit contains the 32 retained edges and their primary evidence.

Rows rejected during audit are retained in the workbook with Pathway_Position = Excluded so the decisions remain traceable.

## Rebuild and audit

Create the environment and run the pipeline from the repository root:

    conda env create -f environment.yml
    conda activate myrf-pathway
    python scripts/01_excel_to_raw.py
    python scripts/02_build_pathway_csv.py
    python scripts/03_preview_kegg_graph.py
    python scripts/04_export_kgml.py
    python scripts/07_audit_literature.py

The final audit verifies the workbook, literature register, classified CSVs, SVG/DOT, KGML, evidence export, and synchronized custom_kegg/hsa99999.xml.

## Evidence boundaries

Moderate edges intentionally preserve limitations in the primary literature:

- MAG, ANGPTL2, and FYN increase MYRF protein through a signaling axis whose final mechanism is unresolved.
- MAPK1 and MAPK3 represent jointly tested ERK1/2 activity; individual contributions were not resolved.
- SOX8 is supported by partial genetic redundancy with SOX10, not direct MYRF enhancer binding.
- MPZ and EGR2 repression is supported in heterologous reporter assays rather than in-vivo gene repression.

## KGML and Pathview

output/myrf_pathway.kgml is a custom KGML-like pathway, not an official KEGG pathway. The exporter automatically keeps custom_kegg/hsa99999.xml synchronized. The Python audit verifies XML structure and relation parity.

Pathview rendering remains optional and requires an external R/Bioconductor environment. Legacy Pathview PDFs were removed because they represented an older 34-relation model and no longer matched the audited graph.
