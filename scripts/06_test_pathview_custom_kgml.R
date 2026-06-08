# 05_test_pathview_custom_kgml.R
# Tests whether Pathview can render the custom MYRF KGML in graphviz mode.
# Appends results to output/pathview_test_report.txt.
# Does NOT modify source CSVs or the master KGML.

KGML_SRC    <- "output/myrf_pathway.kgml"
CUSTOM_DIR  <- "custom_kegg"
CUSTOM_KGML <- file.path(CUSTOM_DIR, "hsa99999.xml")
REPORT_PATH <- "output/pathview_test_report.txt"

RSCRIPT     <- normalizePath(file.path(R.home("bin"), "Rscript"), mustWork = FALSE)

# ── Report helpers ────────────────────────────────────────────────────────────
dir.create("output",     showWarnings = FALSE)
dir.create(CUSTOM_DIR,   showWarnings = FALSE)

rpt <- file(REPORT_PATH, open = "at", encoding = "UTF-8")

log <- function(...) {
  msg <- paste0(...)
  message(msg)
  cat(msg, "\n", file = rpt)
}

log("")
log("=================================================================")
log("Part 2: Pathview Compatibility Test")
log("Generated : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
log("=================================================================")
log("")

# ── Setup: copy KGML ─────────────────────────────────────────────────────────
log("--- Setup ---")
if (!file.exists(KGML_SRC)) {
  log("ERROR: Source KGML not found: ", KGML_SRC)
  close(rpt)
  stop("Run 04_export_kgml.py first.")
}
file.copy(KGML_SRC, CUSTOM_KGML, overwrite = TRUE)
log("Copied KGML  : ", KGML_SRC, " -> ", CUSTOM_KGML)
log("KGML size    : ", file.info(CUSTOM_KGML)$size, " bytes")
log("")

# ── Gene data: Symbol → Entrez mapping ───────────────────────────────────────
log("--- Gene Data Preparation ---")

# All non-custom gene symbols from the pathway
GENE_SYMBOLS <- c(
  "MBP", "MAG", "MOG", "PLP1", "CNTN2", "AATK", "RFFL", "NFASC",
  "GJC2", "GJB1", "DUSP15", "TF",
  "SOX10", "SOX8", "FYN", "MAPK3", "MAPK1", "FBXW7", "GSK3B", "TMEM98",
  "MPZ", "EGR2", "TGFB2", "ID4", "WNT7A", "CSPG4", "PDGFRA", "MYRF"
)

entrez_vec <- NULL
entrez_ok  <- FALSE

tryCatch({
  suppressPackageStartupMessages({
    library(org.Hs.eg.db)
    library(AnnotationDbi)
  })
  entrez_raw <- mapIds(
    org.Hs.eg.db,
    keys     = GENE_SYMBOLS,
    column   = "ENTREZID",
    keytype  = "SYMBOL",
    multiVals = "first"
  )
  entrez_mapped <- entrez_raw[!is.na(entrez_raw)]
  log("Entrez mapping: ", length(entrez_mapped), " / ", length(GENE_SYMBOLS), " symbols resolved")
  for (sym in names(entrez_mapped)) {
    log("  ", formatC(sym, width = 10, flag = "-"), "-> ", entrez_mapped[sym])
  }
  entrez_vec <- setNames(rep(1, length(entrez_mapped)), entrez_mapped)
  entrez_ok  <- TRUE
}, error = function(e) {
  log("Entrez mapping failed: ", conditionMessage(e))
})

# Symbol-keyed gene.data (for gene.idtype = "SYMBOL" mode)
gene_data_sym <- setNames(rep(1, length(GENE_SYMBOLS)), GENE_SYMBOLS)
log("")
log("gene.data (SYMBOL, ", length(gene_data_sym), " genes): ",
    paste(head(names(gene_data_sym), 6), collapse = ", "), ", ...")
if (!is.null(entrez_vec))
  log("gene.data (Entrez, ", length(entrez_vec), " genes): ",
      paste(head(names(entrez_vec), 6), collapse = ", "), ", ...")
log("")

# ── Native PNG check ─────────────────────────────────────────────────────────
log("--- Native PNG Check ---")
PNG_PATH <- file.path(CUSTOM_DIR, "hsa99999.png")
if (file.exists(PNG_PATH)) {
  log("Native PNG found  : ", PNG_PATH)
  log("kegg.native=TRUE  : will be attempted")
  TRY_NATIVE <- TRUE
} else {
  log("Native PNG        : NOT FOUND at ", PNG_PATH)
  log("Native Pathview mode was skipped because custom KGML does not have",
      " a matching KEGG-style PNG background.")
  TRY_NATIVE <- FALSE
}
log("")

# ── Helper: run pathview and capture all output ───────────────────────────────
run_pathview <- function(gene.data, pathway.id, species, kegg.dir,
                         kegg.native, out.suffix, gene.idtype = "ENTREZID") {
  suppressPackageStartupMessages(library(pathview))

  captured <- character(0)

  result <- tryCatch(
    withCallingHandlers(
      pathview(
        gene.data   = gene.data,
        pathway.id  = pathway.id,
        species     = species,
        kegg.dir    = kegg.dir,
        kegg.native = kegg.native,
        out.suffix  = out.suffix,
        gene.idtype = gene.idtype,
        same.layer  = FALSE,
        plot.col.key = FALSE
      ),
      message = function(m) {
        captured <<- c(captured, paste0("[msg]  ", trimws(conditionMessage(m))))
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        captured <<- c(captured, paste0("[warn] ", trimws(conditionMessage(w))))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      captured <<- c(captured, paste0("[err]  ", conditionMessage(e)))
      NULL
    }
  )

  list(result = result, log = captured, ok = !is.null(result))
}

# ── Test 1: kegg.native=FALSE, gene.idtype="SYMBOL" ──────────────────────────
log("--- Pathview Test 1: kegg.native=FALSE  |  gene.idtype=SYMBOL ---")
t1 <- run_pathview(
  gene.data   = gene_data_sym,
  pathway.id  = "99999",
  species     = "hsa",
  kegg.dir    = CUSTOM_DIR,
  kegg.native = FALSE,
  out.suffix  = "myrf_custom_symbol",
  gene.idtype = "SYMBOL"
)
for (line in t1$log) log("  ", line)
log("Result: ", ifelse(t1$ok, "RENDERED (pathview returned non-NULL)", "NO OUTPUT (NULL returned)"))
log("")

# ── Test 2: kegg.native=FALSE, Entrez IDs ────────────────────────────────────
t2_ok <- FALSE
if (entrez_ok && !is.null(entrez_vec)) {
  log("--- Pathview Test 2: kegg.native=FALSE  |  Entrez IDs ---")
  t2 <- run_pathview(
    gene.data   = entrez_vec,
    pathway.id  = "99999",
    species     = "hsa",
    kegg.dir    = CUSTOM_DIR,
    kegg.native = FALSE,
    out.suffix  = "myrf_custom_entrez",
    gene.idtype = "ENTREZID"
  )
  for (line in t2$log) log("  ", line)
  log("Result: ", ifelse(t2$ok, "RENDERED (pathview returned non-NULL)", "NO OUTPUT (NULL returned)"))
  t2_ok <- t2$ok
  log("")
}

# ── Test 3: kegg.native=TRUE (only if PNG exists) ────────────────────────────
t3_ok <- FALSE
if (TRY_NATIVE) {
  log("--- Pathview Test 3: kegg.native=TRUE ---")
  t3 <- run_pathview(
    gene.data   = gene_data_sym,
    pathway.id  = "99999",
    species     = "hsa",
    kegg.dir    = CUSTOM_DIR,
    kegg.native = TRUE,
    out.suffix  = "myrf_custom_native",
    gene.idtype = "SYMBOL"
  )
  for (line in t3$log) log("  ", line)
  log("Result: ", ifelse(t3$ok, "RENDERED", "FAILED"))
  t3_ok <- t3$ok
  log("")
}

# ── Output files ─────────────────────────────────────────────────────────────
log("--- Output Files ---")
pv_files <- list.files(".", pattern = "hsa99999.*myrf_custom.*", full.names = FALSE)
if (length(pv_files) > 0) {
  log("Pathview output files found:")
  for (f in pv_files) log("  ", f)
} else {
  log("No pathview output image files detected in working directory.")
}
log("")

# ── Check if KGML is parseable (quick re-check for classification) ─────────
kgml_parseable <- tryCatch({
  suppressPackageStartupMessages(library(xml2))
  d <- read_xml(CUSTOM_KGML)
  length(xml_find_all(d, ".//entry")) > 0
}, error = function(e) FALSE)

pv_rendered <- t1$ok || t2_ok || t3_ok

# ── Classification ────────────────────────────────────────────────────────────
log("=================================================================")
log("RESULT CLASSIFICATION")
log("=================================================================")
log("")

if      ( kgml_parseable &&  pv_rendered) {
  cls  <- "A"
  desc <- "KGML parses successfully AND Pathview renders successfully."
} else if ( kgml_parseable && !pv_rendered) {
  cls  <- "B"
  desc <- "KGML parses successfully BUT Pathview does not render."
} else if (!kgml_parseable) {
  cls  <- "C"
  desc <- "KGML does not parse successfully."
} else {
  cls  <- "D"
  desc <- "Pathview cannot use custom KGML without additional adaptation."
}

log("Classification : ", cls)
log("Description    : ", desc)
log("")

if (!pv_rendered) {
  log("Likely reasons Pathview could not render:")
  log("")
  log("  1. Non-KEGG entry name prefix")
  log("     Entry names use 'custom:SYMBOL' (e.g. custom:MBP).")
  log("     Pathview expects 'hsa:ENTREZID' (e.g. hsa:4547).")
  log("     Pathview strips the species prefix and interprets the rest as an")
  log("     Entrez gene ID. 'custom:' is not a recognised prefix, so no genes")
  log("     can be mapped and the pathway graph may not be built.")
  log("")
  log("  2. Non-standard node IDs")
  log("     'custom:MYRF_gene' and 'custom:MYRF_protein' are synthetic split")
  log("     nodes with no counterpart in KEGG or Entrez. Pathview will ignore")
  log("     them or raise an error when attempting ID conversion.")
  log("")
  log("  3. Missing official KEGG pathway ID")
  log("     'hsa99999' does not exist in KEGG. If Pathview attempts to download")
  log("     it (e.g. as a fallback), the request will fail with a 404 / empty")
  log("     response.")
  log("")
  log("  4. Missing native PNG background")
  log("     kegg.native=TRUE requires custom_kegg/hsa99999.png. No such file")
  log("     exists, so native mode cannot render.")
  log("")
  log("  5. Unsupported entry name format in Pathview's node.map() pipeline")
  log("     Pathview's internal gene-to-node mapping (node.map()) calls")
  log("     id2eg() which handles 'hsa:', 'ko:', 'ec:' prefixes only.")
  log("     'custom:' causes the lookup to return empty results or an error.")
}

log("")
log("--- Recommended Next Steps ---")
log("")
log("  Option A  (Pathview-compatible KGML)")
log("    1. Add Entrez IDs to data/final_nodes.csv.")
log("       Use: AnnotationDbi::mapIds(org.Hs.eg.db, keys=SYMBOL,")
log("                                  column='ENTREZID', keytype='SYMBOL')")
log("    2. Regenerate KGML with entry name = 'hsa:ENTREZID' for real genes.")
log("    3. Represent MYRF_gene / MYRF_protein as a single 'hsa:745' entry")
log("       (MYRF Entrez ID per org.Hs.eg.db); use node.data for split display.")
log("    4. Re-run this script to confirm Pathview renders.")
log("")
log("  Option B  (Custom rendering pipeline, no Pathview)")
log("    1. Use KEGGgraph + igraph to build the graph object.")
log("    2. Render with ggraph / Rgraphviz directly from the KGML XML.")
log("    3. Full control over custom:MYRF_gene / MYRF_protein split nodes.")
log("    4. No dependency on official KEGG pathway IDs.")
log("")
log("  Option C  (Python pipeline, already in place)")
log("    scripts/03_preview_graph.py already produces a correct SVG via")
log("    pydot/Graphviz using the same curated data. This is the most")
log("    robust visualisation path for this non-standard pathway.")

log("")
log("=================================================================")
log("END OF REPORT")
log("=================================================================")

close(rpt)

# ── Console summary ───────────────────────────────────────────────────────────
cat("\n--- Final Summary ---\n")
cat("KGML parseable  :", ifelse(kgml_parseable, "YES", "NO"), "\n")
cat("Pathview rendered:", ifelse(pv_rendered,   "YES", "NO"), "\n")
cat("Classification  :", cls, "-", desc, "\n")
cat("\nPathview output files in working directory:\n")
if (length(pv_files) > 0) for (f in pv_files) cat(" ", f, "\n") else cat("  (none)\n")
cat("\nFull report appended to:", REPORT_PATH, "\n")
cat("\nWARNING: This is a custom KGML-like pathway and may not be accepted\n")
cat("directly by all KEGG/Pathview workflows without adaptation.\n")
