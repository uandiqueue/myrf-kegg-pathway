# 04_test_kgml_parse.R
# Validates the structural integrity of the custom MYRF KGML file.
# Writes results to output/pathview_test_report.txt.
# Does NOT modify the KGML or any CSV files.

KGML_PATH   <- "output/myrf_pathway.kgml"
REPORT_PATH <- "output/pathview_test_report.txt"

# ── Report helpers ────────────────────────────────────────────────────────────
dir.create("output", showWarnings = FALSE)
rpt <- file(REPORT_PATH, open = "wt", encoding = "UTF-8")

log <- function(...) {
  msg <- paste0(...)
  message(msg)
  cat(msg, "\n", file = rpt)
}

log("=================================================================")
log("MYRF Custom KGML Pathway  -  Parse & Compatibility Test Report")
log("Generated : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
log("KGML file : ", normalizePath(KGML_PATH, mustWork = FALSE))
log("=================================================================")
log("")

# ── Part 1: File existence ────────────────────────────────────────────────────
log("--- 1. File Check ---")
if (!file.exists(KGML_PATH)) {
  log("ERROR: KGML file not found: ", KGML_PATH)
  close(rpt)
  stop("KGML file not found. Run scripts/04_export_kgml.py first.")
}
log("File found  : ", KGML_PATH)
log("File size   : ", file.info(KGML_PATH)$size, " bytes")
log("")

# ── Part 2: xml2 parse ────────────────────────────────────────────────────────
log("--- 2. xml2 Parse Test ---")
xml2_ok      <- FALSE
n_entries    <- 0L
n_relations  <- 0L
entry_id_set <- integer(0)

tryCatch({
  suppressPackageStartupMessages(library(xml2))

  doc  <- read_xml(KGML_PATH)
  root <- xml_name(doc)
  log("Root element         : ", root)
  log("Root check           : ", ifelse(root == "pathway", "PASSED", paste0("FAILED (got '", root, "')")))
  log("pathway/@name        : ", xml_attr(doc, "name"))
  log("pathway/@title       : ", xml_attr(doc, "title"))
  log("pathway/@org         : ", xml_attr(doc, "org"))
  log("pathway/@number      : ", xml_attr(doc, "number"))
  log("")

  # ── Entries ──
  entries     <- xml_find_all(doc, ".//entry")
  n_entries   <- length(entries)
  entry_ids   <- suppressWarnings(as.integer(xml_attr(entries, "id")))
  entry_names <- xml_attr(entries, "name")
  entry_types <- xml_attr(entries, "type")
  entry_id_set <- na.omit(entry_ids)

  log("Entry count          : ", n_entries)
  log("Entry ID range       : ", min(entry_id_set), " – ", max(entry_id_set))
  log("Unique entry types   : ", paste(sort(unique(entry_types)), collapse = ", "))
  log("")

  custom_n <- sum(grepl("^custom:", entry_names))
  hsa_n    <- sum(grepl("^hsa:",    entry_names))
  log("  custom: entries    : ", custom_n,
      ifelse(custom_n > 0, "  ← not standard KEGG; Pathview will not auto-map these", ""))
  log("  hsa:    entries    : ", hsa_n)
  log("")

  # ── Graphics ──
  gfx        <- xml_find_all(doc, ".//graphics")
  miss_x     <- sum(is.na(xml_attr(gfx, "x")))
  miss_y     <- sum(is.na(xml_attr(gfx, "y")))
  miss_name  <- sum(is.na(xml_attr(gfx, "name")))
  miss_ensembl <- sum(is.na(xml_attr(gfx, "ensembl_id")))

  log("Graphics elements    : ", length(gfx))
  log("  missing x          : ", miss_x)
  log("  missing y          : ", miss_y)
  log("  missing name       : ", miss_name)
  log("  have ensembl_id    : ", length(gfx) - miss_ensembl,
      " / ", length(gfx))
  log("  Graphics check     : ",
      ifelse(miss_x + miss_y + miss_name == 0, "PASSED", "ISSUES FOUND"))
  log("")

  # ── Relations ──
  rels        <- xml_find_all(doc, ".//relation")
  n_relations <- length(rels)
  rel_e1      <- suppressWarnings(as.integer(xml_attr(rels, "entry1")))
  rel_e2      <- suppressWarnings(as.integer(xml_attr(rels, "entry2")))
  rel_types   <- xml_attr(rels, "type")

  log("Relation count       : ", n_relations)
  log("Relation types:")
  for (t in sort(unique(rel_types))) {
    log("  ", t, " : ", sum(rel_types == t, na.rm = TRUE))
  }
  log("")

  # ── Subtypes ──
  subs      <- xml_find_all(doc, ".//subtype")
  sub_names <- xml_attr(subs, "name")
  sub_vals  <- xml_attr(subs, "value")
  log("Subtype counts:")
  for (s in sort(unique(sub_names))) {
    idx <- sub_names == s
    log("  ", formatC(s, width = 22, flag = "-"),
        ": ", sum(idx), "  (value: ",
        paste(unique(sub_vals[idx]), collapse = " | "), ")")
  }
  log("")

  # ── Entry reference validation ──
  bad_e1 <- setdiff(na.omit(rel_e1), entry_id_set)
  bad_e2 <- setdiff(na.omit(rel_e2), entry_id_set)
  log("Relation ID reference check:")
  if (length(bad_e1) == 0 && length(bad_e2) == 0) {
    log("  All entry1/entry2 IDs resolve to known entries : PASSED")
  } else {
    if (length(bad_e1) > 0)
      log("  Unresolved entry1 IDs : ", paste(bad_e1, collapse = ", "))
    if (length(bad_e2) > 0)
      log("  Unresolved entry2 IDs : ", paste(bad_e2, collapse = ", "))
  }
  log("")

  # ── MYRF-specific checks ──
  log("MYRF-specific checks:")
  has_gene <- any(grepl("MYRF_gene",    entry_names))
  has_prot <- any(grepl("MYRF_protein", entry_names))
  gene_id  <- entry_ids[grepl("MYRF_gene",    entry_names)]
  prot_id  <- entry_ids[grepl("MYRF_protein", entry_names)]
  log("  MYRF_gene entry    : ", ifelse(has_gene, paste0("FOUND  (id=", gene_id, ")"), "MISSING"))
  log("  MYRF_protein entry : ", ifelse(has_prot, paste0("FOUND  (id=", prot_id, ")"), "MISSING"))

  if (has_gene && has_prot) {
    link_exists <- any((rel_e1 == gene_id & rel_e2 == prot_id), na.rm = TRUE)
    log("  MYRF_gene→protein  : ", ifelse(link_exists, "FOUND", "MISSING"))
  }
  log("")

  xml2_ok <- TRUE
  log("xml2 parse result    : PASSED")

}, error = function(e) {
  log("xml2 parse result    : FAILED")
  log("Error: ", conditionMessage(e))
})

log("")

# ── Part 3: XML package parse ─────────────────────────────────────────────────
log("--- 3. XML Package Parse Test ---")
XML_ok <- FALSE
tryCatch({
  suppressPackageStartupMessages(library(XML))
  xdoc  <- xmlParse(KGML_PATH)
  xroot <- xmlRoot(xdoc)
  log("XML package parse    : PASSED")
  log("Root node name       : ", xmlName(xroot))
  log("Entries  (XML pkg)   : ", length(getNodeSet(xdoc, "//entry")))
  log("Relations (XML pkg)  : ", length(getNodeSet(xdoc, "//relation")))
  XML_ok <- TRUE
  free(xdoc)
}, error = function(e) {
  log("XML package parse    : FAILED")
  log("Error: ", conditionMessage(e))
})

log("")

# ── Part 4: KEGGgraph parse ───────────────────────────────────────────────────
log("--- 4. KEGGgraph Parse Test ---")
kegg_ok    <- FALSE
kegg_nodes <- 0L
kegg_edges <- 0L
tryCatch({
  suppressPackageStartupMessages(library(KEGGgraph))
  # KEGGgraph::parseKGML expects standard KEGG format (hsa:ENTREZID entries)
  kg <- withCallingHandlers(
    parseKGML(KGML_PATH),
    warning = function(w) {
      log("  KEGGgraph warning  : ", conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  kegg_nodes <- length(nodes(kg))
  kegg_edges <- length(edges(kg))
  log("KEGGgraph parse      : PASSED")
  log("  Nodes              : ", kegg_nodes)
  log("  Edges              : ", kegg_edges)
  kegg_ok <- TRUE
}, error = function(e) {
  log("KEGGgraph parse      : FAILED  (expected for custom:SYMBOL entry names)")
  log("  Error              : ", conditionMessage(e))
})

log("")

# ── Summary ───────────────────────────────────────────────────────────────────
log("--- Summary: KGML Parse ---")
log("xml2 parse           : ", ifelse(xml2_ok,   "PASSED", "FAILED"))
log("XML pkg parse        : ", ifelse(XML_ok,    "PASSED", "FAILED"))
log("KEGGgraph parse      : ", ifelse(kegg_ok,   "PASSED", "FAILED (custom entry names)"))
log("Entries validated    : ", n_entries)
log("Relations validated  : ", n_relations)
log("")
log("Note: KEGGgraph failure on custom:SYMBOL entries is EXPECTED.")
log("      xml2/XML parse success means the KGML is well-formed XML.")
log("")

close(rpt)
cat("\nPart 1 complete. Report written to:", REPORT_PATH, "\n")
cat("xml2 parse:", ifelse(xml2_ok, "PASSED", "FAILED"), "\n")
cat("Entries:", n_entries, " | Relations:", n_relations, "\n")
