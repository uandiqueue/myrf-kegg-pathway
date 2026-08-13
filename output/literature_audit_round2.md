# MYRF pathway literature audit — round 2

Generated: 2026-08-14T03:14:59+08:00

Result: **PASS**

## Scope

This second pass checks the curated edge register against the primary-literature audit table and then verifies that the workbook, CSVs, DOT/SVG, KGML, evidence export, and synchronized Pathview copy all encode the same pathway.

## Automated checks

| Check | Result | Detail |
|---|---:|---|
| edge count | PASS | found 32; expected 32 |
| node count | PASS | found 31; expected 31 |
| subtype distribution | PASS | {'expression': 14, 'indirect effect': 8, 'repression': 7, 'inhibition': 1, 'binding/association': 1, 'ubiquitination': 1} |
| duplicate classified edges | PASS | no duplicate source-target-subtype rows expected |
| manual-review flags | PASS | found 0 |
| primary evidence type | PASS | types=['Primary literature'] |
| source URLs populated | PASS | every retained edge must cite a source |
| audit register exact match | PASS | actual=32 audit=32 |
| primary citations match | PASS | each audit citation occurs on its edge |
| removed overstatements absent | PASS | [] |
| evidence export complete | PASS | rows=32 blank_interactions=0 |
| evidence review flags | PASS | evidence export must preserve resolved review status |
| KGML entry count | PASS | found 31 |
| KGML relation count | PASS | found 32 |
| KGML references resolve | PASS | all relation IDs must resolve |
| KGML subtype parity | PASS | {'inhibition': 1, 'indirect effect': 8, 'expression': 14, 'binding/association': 1, 'repression': 7, 'ubiquitination': 1} |
| Pathview KGML synchronized | PASS | output and custom_kegg copies must be byte-identical |
| SVG edge count | PASS | found 32 |
| negative indirect rendering | PASS | GPR37 and HES5 must render as red dashed inhibitory edges |
| SOX10 association rendering | PASS | MYRF-SOX10 association must be visually undirected |
| FBXW7 degradation rendering | PASS | FBXW7 ubiquitination must retain negative/degradative direction |
| workbook audit sheets | PASS | sheets=['Master_Gene_List', 'Literature_Audit', 'Revision_Log', 'References'] |
| MAPK1 identifier | PASS | MAPK1 must not reuse the MAPK3 Ensembl ID |
| excluded source rows | PASS | found 4 |
| portable reference URLs | PASS | References sheet must not contain institution-only proxy URLs |
| workbook audit row count | PASS | found 32 |

## Final model

- Nodes: 31
- Edges: 32
- Strong edges: 23
- Moderate edges: 9
- Manual-review flags: 0

Subtype counts: expression=14, indirect effect=8, repression=7, inhibition=1, binding/association=1, ubiquitination=1

## Literature-reviewed edge register

| Source | Target | Semantics | Confidence | Primary source | Evidence scope |
|---|---|---|---:|---|---|
| TMEM98 | MYRF_protein | inhibition (negative; protein) | Strong | [primary source](https://pubmed.ncbi.nlm.nih.gov/30249802/) | TMEM98 binding blocks MYRF autocleavage and nuclear translocation. |
| MAG | MYRF_protein | indirect effect (positive; protein) | Moderate | [primary source](https://pmc.ncbi.nlm.nih.gov/articles/PMC9976406/) | MAG-FYN signaling increases MYRF protein but the final mechanism is unresolved. |
| MYRF_protein | MAG | expression (positive; transcription) | Strong | [primary source](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001625) | MYRF ChIP-seq occupancy and enhancer assays support direct MAG activation. |
| GPR37 | MYRF_gene | indirect effect (negative; transcription) | Strong | [primary source](https://www.nature.com/articles/ncomms10884) | GPR37 loss elevates ERK activity and MYRF while reducing HES5. |
| MYRF_protein | TF | expression (positive; transcription) | Strong | [primary source](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001625) | MYRF ChIP-seq occupancy and enhancer assays support transferrin activation. |
| MYRF_protein | CNTN2 | expression (positive; transcription) | Strong | [primary source](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001625) | MYRF ChIP-seq occupancy and enhancer assays support CNTN2 activation. |
| MYRF_protein | PLP1 | expression (positive; transcription) | Strong | [primary source](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001625) | MYRF ChIP-seq occupancy and enhancer assays support PLP1 activation. |
| MYRF_protein | MBP | expression (positive; transcription) | Strong | [primary source](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001625) | MYRF ChIP-seq occupancy and enhancer assays support MBP activation. |
| SOX10 | MYRF_gene | expression (positive; transcription) | Strong | [primary source](https://pmc.ncbi.nlm.nih.gov/articles/PMC3814293/) | SOX10 directly binds and activates the MYRF ECR9 enhancer. |
| SOX10 | MYRF_protein | binding/association (neutral; protein) | Strong | [primary source](https://pmc.ncbi.nlm.nih.gov/articles/PMC3814293/) | MYRF and SOX10 physically associate and cooperatively activate selected regulatory elements. |
| MYRF_protein | AATK | expression (positive; transcription) | Strong | [primary source](https://pmc.ncbi.nlm.nih.gov/articles/PMC7026603/) | AATK is a bound and activated oligodendrocyte-specific MYRF-SOX10 target. |
| HES5 | MYRF_gene | indirect effect (negative; transcription) | Moderate | [primary source](https://www.nature.com/articles/ncomms10884) | HES5 suppresses the SOX10 differentiation axis and is inversely linked to MYRF. |
| MYRF_protein | MOG | expression (positive; transcription) | Strong | [primary source](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001625) | MYRF ChIP-seq and expression perturbation support MOG activation. |
| MYRF_protein | DUSP15 | expression (positive; transcription) | Strong | [primary source](https://onlinelibrary.wiley.com/doi/10.1002/glia.23044) | Functional MYRF and SOX10 sites cooperatively activate the DUSP15 promoter. |
| SOX8 | MYRF_gene | indirect effect (positive; transcription) | Moderate | [primary source](https://pmc.ncbi.nlm.nih.gov/articles/PMC3814293/) | SOX8 partially compensates for SOX10 loss without evidence of direct MYRF enhancer binding. |
| ANGPTL2 | MYRF_protein | indirect effect (positive; protein) | Moderate | [primary source](https://pmc.ncbi.nlm.nih.gov/articles/PMC9976406/) | ANGPTL2 activates the MAG-FYN axis and increases downstream MYRF protein. |
| FYN | MYRF_protein | indirect effect (positive; protein) | Moderate | [primary source](https://pmc.ncbi.nlm.nih.gov/articles/PMC9976406/) | FYN activity is associated with increased MYRF protein without direct transcriptional evidence. |
| MAPK1 | MYRF_gene | indirect effect (positive; transcription) | Moderate | [primary source](https://pmc.ncbi.nlm.nih.gov/articles/PMC4244469/) | Combined ERK1-ERK2 manipulation alters MYRF expression with individual contributions unresolved. |
| MAPK3 | MYRF_gene | indirect effect (positive; transcription) | Moderate | [primary source](https://pmc.ncbi.nlm.nih.gov/articles/PMC4244469/) | Combined ERK1-ERK2 manipulation alters MYRF expression with individual contributions unresolved. |
| MYRF_protein | RFFL | expression (positive; transcription) | Strong | [primary source](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001625) | MYRF occupancy and motif-dependent enhancer assays support direct RFFL activation. |
| MYRF_protein | NFASC | expression (positive; transcription) | Strong | [primary source](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001625) | MYRF occupancy and motif-dependent enhancer assays support direct NFASC activation. |
| MYRF_protein | GJC2 | expression (positive; transcription) | Strong | [primary source](https://pmc.ncbi.nlm.nih.gov/articles/PMC3814293/) | MYRF and SOX10 cooperatively activate the GJC2 regulatory element. |
| MYRF_protein | GJB1 | expression (positive; transcription) | Strong | [primary source](https://pmc.ncbi.nlm.nih.gov/articles/PMC3814293/) | MYRF and SOX10 cooperatively activate the GJB1 regulatory element. |
| MYRF_protein | MPZ | repression (negative; transcription) | Moderate | [primary source](https://pmc.ncbi.nlm.nih.gov/articles/PMC3814293/) | MYRF reduces SOX10-responsive MPZ promoter activity in heterologous reporter assays. |
| MYRF_protein | EGR2 | repression (negative; transcription) | Moderate | [primary source](https://pmc.ncbi.nlm.nih.gov/articles/PMC3814293/) | MYRF reduces SOX10-responsive EGR2 MSE activity in heterologous reporter assays. |
| MYRF_protein | TGFB2 | repression (negative; transcription) | Strong | [primary source](https://pmc.ncbi.nlm.nih.gov/articles/PMC7026603/) | MYRF functionally represses SOX10-dependent TGFB2 without binding the tested regulatory region. |
| MYRF_protein | ID4 | repression (negative; transcription) | Strong | [primary source](https://pmc.ncbi.nlm.nih.gov/articles/PMC7026603/) | MYRF functionally represses SOX10-dependent ID4 without binding the tested regulatory region. |
| MYRF_protein | WNT7A | repression (negative; transcription) | Strong | [primary source](https://pmc.ncbi.nlm.nih.gov/articles/PMC7026603/) | MYRF functionally represses SOX10-dependent WNT7A without binding the tested regulatory region. |
| MYRF_protein | CSPG4 | repression (negative; transcription) | Strong | [primary source](https://pmc.ncbi.nlm.nih.gov/articles/PMC7026603/) | MYRF functionally represses SOX10-dependent CSPG4 during oligodendrocyte differentiation. |
| MYRF_protein | PDGFRA | repression (negative; transcription) | Strong | [primary source](https://pmc.ncbi.nlm.nih.gov/articles/PMC7026603/) | MYRF functionally represses SOX10-dependent PDGFRA during oligodendrocyte differentiation. |
| FBXW7 | MYRF_protein | ubiquitination (negative; protein) | Strong | [primary source](https://pmc.ncbi.nlm.nih.gov/articles/PMC12370921/) | FBXW7 binds MYRF and promotes its ubiquitination and proteasomal degradation. |
| MYRF_gene | MYRF_protein | expression (positive; gene) | Strong | [primary source](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001625) | MYRF precursor expression and autocleavage produce the active nuclear N-terminal transcription factor. |

## Explicit uncertainty boundaries

- MAG, ANGPTL2, and FYN are retained as moderate indirect effects on MYRF protein abundance; the terminal mechanism is unresolved.
- MAPK1 and MAPK3 are retained separately for identifier compatibility, but the experiment supports combined ERK1/2 activity rather than resolved individual contributions.
- SOX8 is an indirect moderate edge based on partial genetic redundancy with SOX10; direct MYRF enhancer binding was not shown.
- MPZ and EGR2 repression is limited to heterologous reporter assays and is therefore moderate.

## Removed or collapsed representations

- Removed the duplicate MAG→MYRF unknown-mechanism edge.
- Removed the unsupported direct FYN→MYRF transcriptional edge.
- Collapsed directional SOX10↔MYRF protein activation into one undirected binding/association edge.
- Excluded the duplicate PLP1 WmN1 row because that enhancer did not show MYRF–SOX10 synergy.
