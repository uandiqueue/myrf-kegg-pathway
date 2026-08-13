# Literature-to-pathway curation web app

## Status

This document defines a side-project scope for a human-in-the-loop web application that converts uploaded scientific literature into an evidence-backed pathway graph and validated XML export.

The product is intended to assist expert curation. It must not present LLM-generated relationships as verified scientific conclusions without human approval.

## Product assessment

- Practicality as a human-reviewed application: high.
- Expected usefulness to researchers: high.
- Practicality as a fully autonomous pathway generator: low.
- Primary value: reducing literature triage and manual data entry while preserving evidence provenance and expert control.

The application should be described as a **literature-to-pathway curation workbench**, not an automatic pathway-discovery system.

## Problem statement

Constructing a literature-supported biological pathway currently requires researchers to:

1. locate and read relevant papers;
2. extract evidence and experimental context;
3. identify and normalize biological entities;
4. infer carefully qualified relationships;
5. record those relationships in structured tables;
6. resolve duplicates, conflicts and unsupported claims;
7. build and inspect a graph;
8. generate a pathway exchange format;
9. audit the graph against the literature.

The proposed application should automate the repetitive preparation work while keeping interpretation, correction and approval with the researcher.

## Goals

The application should:

- accept one or more uploaded scientific papers;
- preserve document identity and provenance;
- extract candidate biological entities and experimental contexts;
- help users select a target entity and biological context;
- find relevant evidence passages within the uploaded corpus;
- infer candidate pathway relationships from those passages;
- display each candidate edge beside its evidence;
- let users edit, approve, reject, merge or defer candidates;
- generate an inspectable pathway graph from approved edges;
- validate the approved graph;
- export deterministic, schema-checked XML and evidence artifacts;
- retain a complete audit history.

## Non-goals

The initial product will not:

- claim that it found every biologically relevant pathway;
- treat co-expression or correlation as causal regulation without supporting evidence;
- approve scientific conclusions without human review;
- search the entire scientific literature automatically;
- generate quantitative kinetic models from qualitative evidence;
- guarantee correct interpretation of unreadable figures, scanned pages or missing supplements;
- redistribute uploaded papers or official KEGG data;
- replace systematic review, peer review or domain-expert judgment.

## Terminology

### Target entity

Use **target entity** instead of **target compound** because a project may center on:

- a chemical compound;
- a gene;
- a protein or modified protein state;
- RNA;
- a protein complex;
- a phenotype;
- a biological process.

The MYRF project, for example, must distinguish the `MYRF_gene` locus from active `MYRF_protein` relationships.

### Biological context

Use **biological context** instead of a single free-text **system** value. Context should be structured into:

- organism;
- tissue;
- cell type;
- disease, treatment or physiological state;
- developmental stage;
- experimental model;
- in-vivo, ex-vivo or in-vitro setting;
- intervention, dose and time point when available.

### Candidate relationship

The application should return **candidate relationships supported by parseable evidence**, not promise to find “all pathways.” A pathway is assembled from reviewed candidate relationships.

## Primary user workflow

```text
Create project and define scope
  -> upload one or more papers
  -> parse documents and report coverage
  -> extract and normalize entities and contexts
  -> user resolves entity matches and selects a target
  -> retrieve target-relevant evidence
  -> infer candidate relationships
  -> assemble a candidate pathway
  -> user reviews every edge against its evidence
  -> validate approved entities and relationships
  -> generate and inspect the graph
  -> run a second consistency audit
  -> export XML and evidence package
```

## Functional scope

### 1. Project setup

The user creates a project and defines:

- project title and research question;
- target domain or pathway scope;
- inclusion and exclusion criteria;
- allowed organisms or biological contexts;
- preferred identifier namespaces;
- intended export format;
- whether the project requires one or two reviewers.

### 2. Paper upload and ingestion

The user can upload one or multiple PDF papers.

For every document, the application should retain:

- original filename;
- file hash;
- title, authors, journal and year when detected;
- DOI, PMID or another persistent identifier;
- upload time and uploader;
- extraction status;
- page count;
- warnings for unreadable pages, missing text, figures or supplements.

Each paper must be processed separately before cross-paper merging so provenance is not lost.

### 3. Document parsing

The ingestion layer should attempt to extract:

- section-aware text;
- page and paragraph locations;
- tables and captions;
- figure captions;
- supplementary references;
- reference list;
- text bounding boxes when supported.

The system should display a coverage report containing:

- successfully parsed pages and sections;
- pages that required OCR;
- unreadable content;
- tables or figures not interpreted;
- supplements that were referenced but not uploaded.

The LLM must not silently infer content from material it could not parse.

### 4. Entity and context extraction

The LLM proposes entities and contexts found in each paper.

Entity records should contain:

- exact mention from the paper;
- preferred label;
- entity type;
- candidate normalized identifiers;
- aliases;
- species;
- document locations;
- match confidence;
- user-selected canonical identity.

Recommended identifier sources include:

- ChEBI, InChIKey or PubChem CID for chemical entities;
- HGNC, NCBI Gene or Ensembl for genes;
- UniProt for proteins;
- NCBI Taxonomy for organisms;
- established biomedical ontologies for tissues, cell types and diseases.

The user must resolve ambiguous matches before pathway generation. The original paper wording must remain available after normalization.

### 5. Target selection

The user selects:

- one primary target entity;
- one or more accepted biological contexts;
- whether related forms or aliases should be included;
- the permitted evidence distance from the target;
- whether indirect relationships are allowed.

The application should show why each proposed entity/context match was suggested.

### 6. Evidence retrieval

The application retrieves evidence units related to the selected target and context.

An evidence unit should contain:

- paper identifier and file hash;
- exact excerpt;
- page, section and paragraph;
- figure, table or supplement reference when applicable;
- entities mentioned;
- experimental context;
- assay and perturbation;
- measured result;
- extraction confidence;
- any parsing warning.

Evidence retrieval should combine deterministic search with semantic retrieval. The application should show relevant negative results and unresolved passages rather than returning only confident-looking findings.

### 7. Candidate relationship extraction

The LLM converts evidence units into structured candidate edges.

Every candidate should contain:

- source entity;
- target entity;
- direction: positive, negative, neutral or unknown;
- regulation level: gene/transcription, RNA, protein, compound or process;
- mechanism: expression, repression, activation, inhibition, indirect effect, binding/association, phosphorylation, dephosphorylation, ubiquitination, degradation, reaction or unknown;
- directness: direct, indirect, associative or unresolved;
- experimental context;
- Strong, Moderate or Weak evidence assessment;
- limitations and alternative interpretations;
- one or more linked evidence units;
- LLM explanation;
- extraction model and version;
- status: Candidate by default.

The LLM must not convert correlation, co-expression, pathway membership or a statement in the introduction into a causal edge unless the cited experiment supports it.

### 8. Evidence review interface

The main review screen should show the candidate graph and an edge review panel side by side.

For each edge, the user must be able to:

- open the supporting passage in document context;
- inspect page, figure, table and supplement references;
- compare evidence from multiple papers;
- see contradictory evidence;
- edit source, target, direction, level, mechanism and confidence;
- add notes;
- approve;
- reject with a reason;
- mark as duplicate;
- merge equivalent evidence;
- mark as ambiguous;
- request additional literature.

Required review states:

- Candidate
- Approved
- Edited and approved
- Rejected
- Duplicate
- Ambiguous
- Needs additional literature

Only Approved and Edited-and-approved edges may enter a release export.

### 9. Conflict and duplicate handling

The application should detect:

- duplicate source-target-mechanism edges;
- the same observation extracted more than once;
- opposing directions under the same context;
- differences caused by organism, cell type, dose or time point;
- gene-versus-protein routing conflicts;
- direct and indirect representations of the same experiment;
- identifier collisions.

Contradictory findings must remain visible. The application must not collapse them into a single consensus edge without a recorded user decision.

### 10. Graph generation

The graph must be generated from the approved canonical data, not directly from LLM prose.

The graph view should distinguish:

- target/core nodes;
- upstream and downstream entities;
- gene and protein forms;
- positive and negative effects;
- binding/association;
- direct and indirect relationships;
- evidence strength;
- ambiguous or pending candidates when the user enables them.

Before XML export, the user should inspect the graph for:

- missing or duplicated nodes;
- wrong source-target direction;
- misleading arrowheads;
- incorrect biological level;
- hidden contradictions;
- unreadable layout;
- approved relationships without visible evidence.

### 11. Automated validation

The validator should check:

- all edge endpoints resolve to known nodes;
- normalized identifiers are valid and not conflicting;
- every approved edge has at least one approved evidence unit;
- no rejected or ambiguous edge enters the release graph;
- relationship type and direction are compatible;
- gene and protein forms are routed correctly;
- duplicate classified edges are absent;
- graph and export counts match;
- provenance fields are populated;
- the selected export schema validates;
- output is reproducible from the approved snapshot.

### 12. Second-round audit

After graph generation, the application should present a release audit containing:

- paper count and parsing coverage;
- approved, rejected, duplicate and ambiguous counts;
- node and edge counts;
- relationship subtype distribution;
- edges missing primary evidence;
- unresolved conflicts;
- graph-versus-export parity;
- schema-validation result;
- reviewer names and timestamps;
- input, model, prompt and exporter versions;
- a cryptographic hash of the approved release snapshot.

Release should be blocked when required evidence, review or schema checks fail.

### 13. Deterministic export

The LLM must not write the final XML directly.

The approved pathway should first be stored in a canonical internal representation, such as structured JSON or relational graph tables. Deterministic exporters should then generate:

- custom KGML for compatibility with the current MYRF project;
- evidence CSV;
- canonical JSON;
- node/edge CSV files;
- SVG or another inspectable graph format;
- an audit report.

Potential later exporters:

- BioPAX Level 3 for rich pathway exchange;
- SBML Qual for qualitative regulatory models;
- SBML Core for reaction-based or quantitative models;
- GraphML for general graph tooling.

Custom KGML output must be labelled as custom rather than represented as an official KEGG pathway.

## Canonical data model

### Paper

- `paper_id`
- `file_hash`
- `filename`
- `title`
- `authors`
- `year`
- `doi`
- `pmid`
- `parse_status`
- `coverage_warnings`

### Entity

- `entity_id`
- `raw_mention`
- `preferred_label`
- `entity_type`
- `namespace`
- `accession`
- `aliases`
- `species`
- `normalization_status`
- `reviewed_by`

### Context

- `context_id`
- `organism`
- `tissue`
- `cell_type`
- `disease_or_state`
- `developmental_stage`
- `experimental_model`
- `intervention`
- `dose`
- `time_point`

### Evidence unit

- `evidence_id`
- `paper_id`
- `excerpt`
- `page`
- `section`
- `paragraph_or_bbox`
- `figure_or_table`
- `context_id`
- `assay`
- `perturbation`
- `measured_result`
- `extraction_confidence`
- `parsing_warning`

### Relationship

- `relationship_id`
- `source_entity_id`
- `target_entity_id`
- `direction`
- `regulation_level`
- `mechanism`
- `directness`
- `confidence`
- `context_id`
- `limitations`
- `status`
- linked `evidence_ids`

### Review decision

- `decision_id`
- `object_type`
- `object_id`
- `previous_value`
- `new_value`
- `decision`
- `reason`
- `reviewer`
- `timestamp`

### Release

- `release_id`
- `project_id`
- approved entity and relationship snapshot;
- paper hashes;
- model and prompt versions;
- exporter and schema versions;
- validator results;
- release hash;
- creator and timestamp.

## MVP

### Included

- authenticated single-user projects;
- upload of one or more born-digital PDFs;
- section-aware text extraction;
- entity and context proposal;
- user entity matching and target selection;
- paragraph and table-caption evidence retrieval;
- candidate relationship extraction;
- edge-by-edge evidence review;
- edit, approve, reject, duplicate and ambiguous states;
- canonical evidence graph;
- graph preview;
- deterministic custom KGML, JSON, CSV and SVG export;
- automated parity and provenance checks;
- downloadable audit report.

### Deferred

- high-quality scientific figure interpretation;
- robust OCR for poor scans;
- automatic supplement downloading;
- external literature searching;
- collaborative organizations and complex permissions;
- two-reviewer adjudication;
- contradiction-resolution recommendations;
- BioPAX and SBML export;
- quantitative or kinetic model generation;
- official KEGG content integration;
- public pathway publishing.

## Suggested delivery phases

### Phase 1: Evidence extraction prototype

- PDF ingestion and coverage reporting;
- entity/context extraction;
- evidence retrieval;
- structured candidate edge output;
- no XML generation yet.

### Phase 2: Review workbench

- entity matching;
- edge review states;
- evidence viewer;
- editing, deduplication and conflict display;
- canonical graph storage;
- complete decision log.

### Phase 3: Graph and release pipeline

- graph preview;
- deterministic custom KGML exporter;
- evidence and audit exports;
- schema, parity and provenance validation;
- versioned release snapshots.

### Phase 4: Broader interoperability

- improved OCR, tables, figures and supplements;
- multi-user review;
- external ontology and literature integrations;
- BioPAX, SBML Qual and GraphML exporters;
- benchmark-driven model improvements.

## Quality and safety requirements

- No approved edge without traceable evidence.
- No final XML generated directly from unreviewed LLM output.
- No silent omission of unreadable document content.
- No silent merging of contradictory evidence.
- No conversion of association to causation without support.
- No loss of original wording after entity normalization.
- No mutation of an approved release; corrections create a new release.
- Uploaded papers must be access-controlled and deletable.
- Model, prompt and extraction versions must be recorded.
- User edits must always override model suggestions.

## Risks and mitigations

| Risk | Consequence | Mitigation |
|---|---|---|
| Poor PDF parsing | Missing or distorted evidence | Coverage report, OCR warnings and page-level inspection |
| Hallucinated mechanism | Unsupported pathway edge | Evidence-required structured extraction and mandatory approval |
| Entity ambiguity | Wrong compound, gene or protein | Ontology matching with explicit user resolution |
| Context collapse | Incorrect generalization | Structured organism, tissue, cell, state, dose and time fields |
| Multi-paper contradiction | False consensus | Preserve evidence separately and expose conflicts |
| Duplicate extraction | Inflated pathway support | Evidence-level and edge-level deduplication |
| Gene/protein conflation | Incorrect biological routing | Explicit entity forms and regulation levels |
| LLM output drift | Non-reproducible candidates | Record model/prompt version and freeze approved releases |
| Invalid XML | Tool incompatibility | Deterministic exporter and schema validation |
| Copyright or licensing issues | Improper document/data use | Private processing, no redistribution and licensing review for external databases |
| Sensitive unpublished papers | Data exposure | Encryption, retention controls and explicit deletion |

## Acceptance criteria for the MVP

The MVP is acceptable when:

1. a user can upload multiple papers and see parsing coverage for each one;
2. the system proposes entities and contexts with traceable source locations;
3. the user can resolve ambiguous target identities;
4. every candidate edge links to at least one exact evidence passage;
5. the user can edit, approve, reject, merge or defer every edge;
6. only approved edges appear in the release graph;
7. contradictions and context differences remain visible;
8. graph generation precedes XML export and can block release;
9. XML is generated deterministically from the approved canonical graph;
10. graph, evidence table and XML contain the same approved relationships;
11. the exported release includes paper hashes, decisions and validator results;
12. rerunning export from the same approved release produces an equivalent result.

## Benchmark using the MYRF project

The audited MYRF pathway can serve as the first golden evaluation set.

The application should be tested on whether it can:

- recover the 32 approved relationships as evidence-linked candidates;
- avoid or flag the four excluded source rows;
- distinguish `MYRF_gene` from `MYRF_protein`;
- model FBXW7 as protein ubiquitination/degradation rather than gene repression;
- model MYRF-SOX10 as association rather than two directional activation edges;
- retain moderate uncertainty for MAG, FYN, ANGPTL2, SOX8, MAPK1 and MAPK3;
- preserve evidence and audit decisions through graph and XML export;
- generate graph and XML artifacts with exact relationship parity.

The benchmark should measure candidate recall, unsupported-edge rate, identifier accuracy, context accuracy, reviewer edit rate and final artifact parity. It should not reward the model for producing a larger graph.

## Standards and references

- [KEGG Markup Language documentation](https://www.kegg.jp/kegg/xml/docs/)
- [BioPAX Level 3 documentation](https://www.biopax.org/release/biopax-level3-documentation.pdf)
- [SBML specifications](https://sbml.org/documents/specifications/)
- [ChEBI tools and REST API](https://www.ebi.ac.uk/chebi/tools)

## Final product principle

> The LLM prepares a traceable candidate pathway; the researcher decides what the evidence supports.
