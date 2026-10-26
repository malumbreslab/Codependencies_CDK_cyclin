# Manuscript:   A reference map of CDK-cyclin dependencies across human cancers
# Project:      Codependencies
# Part:         H. Extra analysis
# Author:       Agustín Sánchez-Belmonte
# Date:         26/10/26
# Last update:  26/10/26


# ---------------- INDEX ---------------- 
# 1. Confounder-adjusted CDK-cyclin codependencies
# 2. Validation of the cluster score classifier (CSC) in DepMap
# 3. Palbociclib (PRISM): differential expression of resistant vs sensitive cell lines within each dependency cluster
# 4. Palbociclib sensitivity vs CDK4/CDK6 genetic dependency (DE models)
# 5. Palbociclib AUC vs CRISPR gene effect of CDK4, CDK6 and CCND1
# 6. Expression correlates of palbociclib sensitivity vs CDK4/CDK6/CCND1 dependency (gene-wise Spearman correlations)
# =============================================================================


# 0. Configuration -------------------------------------------------------------

# Analyses to run
RUN <- list(
  adjusted_codependency = TRUE,   # 1
  classifier_validation = TRUE,   # 2
  prism_de_by_cluster   = TRUE,   # 3
  prism_drug_vs_crispr  = TRUE,   # 4
  prism_auc_vs_crispr   = TRUE,   # 5
  prism_cor_expression  = TRUE    # 6
)

# Input and output directories (relative to the repository root)
DIRS <- list(
  depmap = "data/DepMap_25Q3",
  prism  = "data/PRISM_19Q4",
  tables = "data/supplementary_tables",
  out    = "results"
)

FILES <- list(
  # DepMap Public 25Q3
  expression  = file.path(DIRS$depmap, "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"),
  gene_effect = file.path(DIRS$depmap, "CRISPRGeneEffect.csv"),
  model       = file.path(DIRS$depmap, "Model.csv"),
  copy_number = file.path(DIRS$depmap, "OmicsCNGeneWGS.csv"),   # gene-level relative copy number
  # PRISM Repurposing 19Q4, secondary screen
  prism_auc       = file.path(DIRS$prism, "secondary-screen-dose-response-curve-parameters.csv"),
  prism_lfc       = file.path(DIRS$prism, "secondary-screen-replicate-collapsed-logfold-change.csv"),
  prism_treatment = file.path(DIRS$prism, "secondary-screen-replicate-collapsed-treatment-info.csv"),
  prism_cell_info = file.path(DIRS$prism, "secondary-screen-cell-line-info.csv"),
  # Supplementary tables of the paper
  table_s3 = file.path(DIRS$tables, "TableS3_Cluster_assigment_cell_lines.xlsx"),
  table_s9 = file.path(DIRS$tables, "TableS10_DepMap_Cluster_markers.xlsx")
)

# Dependency clusters (Table S3) and their colours
F_SUB    <- paste0("F", 1:6)
CLUSTERS <- c("A", "B", "C", "D", "E", F_SUB)
MACRO    <- c("A", "B", "C", "D", "E", "F")
CLUSTER_COLORS <- c("A"  = "#CD3333", "B"  = "#1874CD", "C"  = "#CD9B1D", "D"  = "#7A378B",
                    "E"  = "#EE9572",
                    "F1" = "#17301D", "F2" = "#264F30", "F3" = "#3E804F", "F4" = "#72BB84",
                    "F5" = "#B0D9BA", "F6" = "#CFE8D5")

# Publication figure style (Analysis 2)
FONT        <- "Helvetica"   # "Arial" if required by the journal
BASE_PT     <- 7             # base text size (pt)
SHOW_TITLES <- FALSE         # TRUE: working versions with title and subtitle
W1 <- 89; W15 <- 120; W2 <- 183   # figure widths (mm): 1, 1.5 and 2 columns
COL_MAIN <- "#534AB7"; COL_LIGHT <- "#9F98E8"; COL_GREY <- "#B4B2A9"
COL_HEAT_LO <- "#F4F8FC"; COL_HEAT_HI <- "#185FA5"


# 1. Packages ------------------------------------------------------------------
# Bioconductor packages are attached first so that the dplyr functions with the
# same names (select, filter, first, rename, desc, ...) take precedence, and
# randomForest is attached before caret/ggplot2 so that ggplot2::margin is not
# masked. AnnotationDbi and org.Hs.eg.db are used through `::`.
suppressPackageStartupMessages({
  library(genefu)          # Bioconductor: PAM50 subtyping
  library(limma)           # Bioconductor
  library(randomForest)
  library(caret)
  library(data.table)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(forcats)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(patchwork)
  library(corrplot)
})


# 2. Shared helpers ------------------------------------------------------------

## 2.1 General utilities ----
msg       <- function(...) cat(sprintf(...), "\n", sep = "")
norm_name <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
safe      <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
is_ach    <- function(v) is.character(v) && mean(grepl("^ACH-[0-9]+$", head(v[!is.na(v)], 500))) > 0.8
has_file  <- function(f) !is.null(f) && nzchar(f) && file.exists(f)
`%||%`    <- function(a, b) if (is.null(a)) b else a

# DepMap column names "GENE (EntrezID)" -> "GENE"
clean_names <- function(x) sub(" .*", "", x)

## 2.2 Figure style and saving ----
pt2mm <- function(pt) pt / ggplot2::.pt

pdf_device <- function() {
  if (capabilities("cairo")) grDevices::cairo_pdf else
    function(...) grDevices::pdf(..., useDingbats = FALSE)
}

# Publication figures: width and height in mm; vector PDF + 600-dpi PNG
save_fig_mm <- function(p, file_base, width, height) {
  ggsave(paste0(file_base, ".pdf"), p, width = width, height = height,
         units = "mm", device = pdf_device())
  ggsave(paste0(file_base, ".png"), p, width = width, height = height,
         units = "mm", dpi = 600, bg = "white")
  invisible(p)
}

# Exploratory figures (Analyses 3-6): width and height in inches; PNG + PDF
save_plot <- function(g, file_base, width, height) {
  ggsave(paste0(file_base, ".png"), g, width = width, height = height, dpi = 300)
  ggsave(paste0(file_base, ".pdf"), g, width = width, height = height)
}

theme_paper <- function(base_size = BASE_PT, base_family = FONT, grid = c("none", "x", "y")) {
  grid <- match.arg(grid)
  th <- theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text              = element_text(colour = "black"),
      axis.text         = element_text(size = base_size - 1, colour = "black"),
      axis.title        = element_text(size = base_size),
      axis.line         = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks        = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks.length = unit(1, "mm"),
      strip.background  = element_blank(),
      strip.text        = element_text(size = base_size, face = "bold",
                                       margin = ggplot2::margin(1, 0, 1, 0, "mm")),
      legend.title      = element_text(size = base_size - 1),
      legend.text       = element_text(size = base_size - 1),
      legend.key.size   = unit(3, "mm"),
      legend.background = element_blank(),
      legend.margin     = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      plot.title        = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle     = element_text(size = base_size - 1, colour = "grey30"),
      plot.tag          = element_text(size = base_size + 2, face = "bold", family = base_family),
      plot.margin       = ggplot2::margin(2, 2, 2, 2, "mm"),
      panel.spacing     = unit(2, "mm")
    )
  if (grid == "x") th <- th + theme(panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey90"))
  if (grid == "y") th <- th + theme(panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey90"))
  th
}

theme_heatmap <- function(...) {
  theme_paper(...) + theme(axis.line = element_blank(), axis.ticks = element_blank())
}

# Titles only in working versions (SHOW_TITLES = TRUE); final figures carry them in the legend
ttl <- function(title = NULL, subtitle = NULL) {
  if (SHOW_TITLES) labs(title = title, subtitle = subtitle) else labs(title = NULL, subtitle = NULL)
}

## 2.3 DepMap 25Q3 matrices (Analyses 1 and 2) ----
# Returns cell lines x genes data frames (row names = ModelID) and Model.csv
load_depmap <- function(files) {
  # Expression: log2(TPM + 1), protein-coding genes, default profile of each model
  msg("Reading %s ...", basename(files$expression))
  expr_raw <- read.csv(files$expression, check.names = FALSE)
  expr_raw <- expr_raw[expr_raw$IsDefaultEntryForModel == "Yes", ]
  colnames(expr_raw) <- clean_names(colnames(expr_raw))
  stopifnot(!anyDuplicated(expr_raw$ModelID))

  expression_data <- expr_raw |>
    dplyr::select(-1, -SequencingID, -IsDefaultEntryForMC,
                  -IsDefaultEntryForModel, -ModelConditionID)
  rownames(expression_data) <- expression_data$ModelID
  expression_data <- expression_data |> dplyr::select(-ModelID)
  expression_data <- expression_data[, !duplicated(colnames(expression_data))]
  rm(expr_raw)

  # CRISPR gene effect (Chronos)
  msg("Reading %s ...", basename(files$gene_effect))
  dep_raw <- read.csv(files$gene_effect, check.names = FALSE)
  colnames(dep_raw) <- clean_names(colnames(dep_raw))
  colnames(dep_raw)[1] <- "ModelID"
  stopifnot(!anyDuplicated(dep_raw$ModelID))

  gene_effect_data <- dep_raw |> dplyr::select(-ModelID)
  rownames(gene_effect_data) <- dep_raw$ModelID
  gene_effect_data <- gene_effect_data[, !duplicated(colnames(gene_effect_data))]
  rm(dep_raw)

  # Cell line annotation
  model_metadata <- read.csv(files$model, check.names = FALSE)
  stopifnot(all(c("ModelID", "OncotreeLineage") %in% colnames(model_metadata)))

  list(expression = expression_data, gene_effect = gene_effect_data, model = model_metadata)
}

# Gene-level copy number (cell lines x genes, row names = ModelID). Handles both
# the 25Q3 layout (ModelID + profile columns) and the older one (ModelID first).
load_copy_number <- function(file) {
  msg("Reading %s ...", basename(file))
  cn_raw <- read.csv(file, check.names = FALSE)
  if ("IsDefaultEntryForModel" %in% colnames(cn_raw))
    cn_raw <- cn_raw[cn_raw$IsDefaultEntryForModel == "Yes", ]
  colnames(cn_raw) <- clean_names(colnames(cn_raw))
  if (!"ModelID" %in% colnames(cn_raw)) colnames(cn_raw)[1] <- "ModelID"
  stopifnot(!anyDuplicated(cn_raw$ModelID))

  meta_cols <- c("ModelID", "SequencingID", "ModelConditionID",
                 "IsDefaultEntryForMC", "IsDefaultEntryForModel")
  cn <- cn_raw[, !(colnames(cn_raw) %in% meta_cols | colnames(cn_raw) == ""), drop = FALSE]
  rownames(cn) <- cn_raw$ModelID
  cn[, !duplicated(colnames(cn)), drop = FALSE]
}

## 2.4 PRISM and DepMap loaders (Analyses 3-6) ----
# These functions read the files given in a cfg list (see PRISM_CFG below).

# Cluster order: palette order first, then the rest (numeric or alphabetical)
cluster_order <- function(x, palette = NULL) {
  x <- unique(as.character(x[!is.na(x)]))
  first <- intersect(names(palette), x)
  rest  <- setdiff(x, first)
  rest  <- if (length(rest) && all(!is.na(suppressWarnings(as.numeric(rest)))))
    rest[order(as.numeric(rest))] else sort(rest)
  c(first, rest)
}

cluster_palette <- function(levels, palette) {
  cols <- unname(palette[levels])
  miss <- levels[is.na(cols)]
  if (length(miss)) {
    msg("  NOTE: clusters without a colour in cfg$cluster_colors: %s (drawn in grey)",
        paste(miss, collapse = ", "))
    cols[is.na(cols)] <- grey.colors(length(miss), start = 0.3, end = 0.7)
  }
  setNames(cols, levels)
}

dose_labels <- function(x) formatC(signif(x, 2), format = "fg", digits = 2)

load_annotation <- function(cfg) {
  if (has_file(cfg$model_file)) {
    m <- fread(cfg$model_file)
    lincol <- intersect(c("OncotreeLineage", "lineage"), names(m))[1]
    if (is.na(lincol)) stop("Model.csv has no OncotreeLineage column.")
    out <- data.table(ModelID = m$ModelID, lineage = as.character(m[[lincol]]))
    for (cc in c("CCLEName", "StrippedCellLineName", "CellLineName"))
      set(out, j = cc, value = if (cc %in% names(m)) as.character(m[[cc]]) else NA_character_)
  } else if (has_file(cfg$cell_info_file)) {
    msg("WARNING: no Model.csv; using %s (primary_tissue) as lineage.", basename(cfg$cell_info_file))
    m <- fread(cfg$cell_info_file)
    out <- data.table(ModelID = m$depmap_id, lineage = as.character(m$primary_tissue),
                      CCLEName = m$ccle_name,
                      StrippedCellLineName = sub("_.*$", "", m$ccle_name),
                      CellLineName = NA_character_)
  } else stop("Model.csv (or secondary-screen-cell-line-info.csv) is required for lineage/breast.")
  out <- out[!is.na(ModelID) & !duplicated(ModelID)]
  out[, is_breast := grepl("^breast$", lineage, ignore.case = TRUE)]
  out
}

# Header without fread(nrows = 0), which fails in some data.table versions
read_header <- function(f) {
  h <- tryCatch(names(fread(f, nrows = 1, showProgress = FALSE)), error = function(e) NULL)
  if (is.null(h)) {
    l <- readLines(f, n = 1, warn = FALSE)
    h <- gsub('^"|"$', "", strsplit(l, ",", fixed = TRUE)[[1]])
  }
  h
}

# Official DepMap AUC (same values as in the rest of the paper)
load_official_auc <- function(cfg) {
  p <- fread(cfg$auc_file, showProgress = FALSE)
  p <- p[tolower(name) == tolower(cfg$drug_name) & screen_id == cfg$screen_id &
           !is.na(auc) & !is.na(depmap_id) & depmap_id != ""]
  if (!nrow(p)) stop("No AUC for ", cfg$drug_name, " in screen ", cfg$screen_id)
  ndup <- p[, .N, by = depmap_id][N > 1, .N]
  msg("  official AUC: %s %s, %d lines, batches=%s%s", cfg$drug_name, cfg$screen_id,
      uniqueN(p$depmap_id), paste(unique(p$broad_id), collapse = ","),
      if (ndup) sprintf(" | %d lines with duplicated measurements -> median", ndup) else "")
  msg("  AUC range: %.3f - %.3f", min(p$auc), max(p$auc))
  p[, .(auc_official = median(auc)), by = depmap_id]
}

# Dose-response curve fitting (4-parameter log-logistic)
# v(x) = L + (U - L) / (1 + 10^(s * (x - x0))),  x = log10(dose uM),  L = f * U
curve4 <- function(p, x) { U <- p[2]; L <- p[1] * U; L + (U - L) / (1 + 10^(p[3] * (x - p[4]))) }

fit_curve4 <- function(x, v, lo, hi) {
  ok <- is.finite(x) & is.finite(v); x <- x[ok]; v <- v[ok]
  if (length(unique(round(x, 2))) < 4) return(NULL)
  sse <- function(p) sum((v - curve4(p, x))^2)
  best <- NULL
  for (f0 in c(0.1, 0.6)) for (s0 in c(0.7, 2)) for (x00 in quantile(x, c(.25, .5, .75))) {
    o <- try(optim(c(f0, min(1.4, max(0.6, max(v))), s0, x00), sse, method = "L-BFGS-B",
                   lower = c(0, 0.5, 0.1, lo - 2), upper = c(1, 1.5, 5, hi + 2)), silent = TRUE)
    if (!inherits(o, "try-error") && (is.null(best) || o$value < best$value)) best <- o
  }
  if (is.null(best)) return(NULL)
  tss <- sum((v - mean(v))^2)
  list(par = best$par, r2 = if (tss > 0) 1 - best$value / tss else NA_real_)
}

fit_auc <- function(x, v, lo, hi) {
  f <- fit_curve4(x, v, lo, hi)
  if (is.null(f)) return(list(auc = NA_real_, r2 = NA_real_))
  list(auc = mean(curve4(f$par, seq(lo, hi, length.out = 200))), r2 = f$r2)
}

trapz_auc <- function(x, v, lo, hi) {
  d <- data.table(x = round(x, 2), v = pmin(pmax(v, 0), 1.5))[, .(v = mean(v)), by = x][order(x)]
  d <- d[x >= round(lo, 2) & x <= round(hi, 2)]
  if (nrow(d) < 2) return(NA_real_)
  sum(diff(d$x) * (head(d$v, -1) + tail(d$v, -1)) / 2) / (max(d$x) - min(d$x))
}

# AUC recomputed from the PRISM log-fold-changes (cfg$auc_source = "lfc") and
# per-dose viability points for the dose-response plots
compute_auc <- function(cfg, annot) {
  hdr <- read_header(cfg$lfc_file)
  meta <- NULL
  if (has_file(cfg$treatment_file)) {
    ti <- fread(cfg$treatment_file)
    ccol <- intersect(c("column_name", "column_id", "col_name"), names(ti))[1]
    if (!is.na(ccol)) {
      hit <- startsWith(as.character(ti$broad_id), cfg$drug_brd)
      if ("name" %in% names(ti)) hit <- hit | tolower(ti$name) %in% tolower(cfg$drug_name)
      meta <- ti[hit, .(col = get(ccol), broad_id, dose = as.numeric(dose), screen_id)]
      meta <- meta[col %in% hdr]
    }
  }
  if (is.null(meta) || !nrow(meta)) {
    cols <- grep(paste0("^", cfg$drug_brd), hdr, value = TRUE)
    if (!length(cols)) stop("No columns for ", cfg$drug_name, " (", cfg$drug_brd, ") in the LFC file.")
    pp <- tstrsplit(cols, "::", fixed = TRUE)
    meta <- data.table(col = cols, broad_id = pp[[1]], dose = as.numeric(pp[[2]]),
                       screen_id = if (length(pp) >= 3) pp[[3]] else "single")
  }
  msg("  %d %s columns: batches=%s | screens=%s", nrow(meta), cfg$drug_name,
      paste(unique(meta$broad_id), collapse = ","), paste(unique(meta$screen_id), collapse = ","))

  raw <- fread(cfg$lfc_file, select = c(1L, match(meta$col, hdr)))
  setnames(raw, 1, "row_id")
  long <- melt(raw, id.vars = "row_id", variable.name = "col", value.name = "lfc",
               variable.factor = FALSE)[!is.na(lfc)]
  long <- merge(long, meta, by = "col")

  # Rows -> ModelID (ACH-); otherwise via CCLE name
  rid <- unique(long$row_id)
  map <- data.table(row_id = rid,
                    depmap_id = ifelse(grepl("ACH-[0-9]+", rid), sub("^.*(ACH-[0-9]+).*$", "\\1", rid), NA))
  if (anyNA(map$depmap_id)) {
    k <- match(norm_name(map$row_id), norm_name(annot$CCLEName))
    map[is.na(depmap_id), depmap_id := annot$ModelID[k[is.na(map$depmap_id)]]]
  }
  long <- merge(long, map[!is.na(depmap_id)], by = "row_id")

  # A single screen: cfg$screen_id if present, otherwise the one with most lines
  sc <- long[, .(n = uniqueN(depmap_id)), by = screen_id][order(-(screen_id == cfg$screen_id), -n)]
  scr <- sc$screen_id[1]
  long <- long[screen_id == scr]
  long[, `:=`(x = log10(dose), v = 2^lfc)]
  msg("  screen used: %s (%d lines). Doses (uM): %s", scr, uniqueN(long$depmap_id),
      paste(signif(sort(unique(long$dose)), 2), collapse = ", "))

  lo <- min(long$x); hi_r <- log10(max(long$dose[long$dose <= cfg$max_dose_uM])); hi_f <- max(long$x)
  res <- long[, {
    r  <- .SD[dose <= cfg$max_dose_uM]
    fr <- fit_auc(r$x, r$v, lo, hi_r)
    ff <- fit_auc(x, v, lo, hi_f)
    .(auc = fr$auc, r2 = fr$r2, auc_trapz = trapz_auc(r$x, r$v, lo, hi_r),
      auc_full_10uM = ff$auc, n_points = nrow(r))
  }, by = depmap_id]
  res[is.na(auc), auc := auc_trapz]            # fallback if the fit fails
  res[, screen_id := scr]

  if (has_file(cfg$cell_info_file)) {
    ci <- fread(cfg$cell_info_file)
    if (all(c("depmap_id", "passed_str_profiling") %in% names(ci))) {
      bad <- ci[toupper(as.character(passed_str_profiling)) == "FALSE", unique(depmap_id)]
      if (length(bad)) { res <- res[!depmap_id %in% bad]; msg("  %d lines removed (failed STR)", length(bad)) }
    }
  }
  msg("  LFC-based AUC for %d lines | cor(fit, trapezoid)=%.2f | cor(without 10uM, with 10uM)=%.2f",
      nrow(res), cor(res$auc, res$auc_trapz, use = "complete"),
      cor(res$auc, res$auc_full_10uM, use = "complete"))
  list(auc = res[!is.na(auc)], points = long)
}

load_expr <- function(f, profile_map_file = NULL) {
  e <- fread(f, check.names = FALSE)
  if ("IsDefaultEntryForModel" %in% names(e))
    e <- e[toupper(as.character(IsDefaultEntryForModel)) %in% c("YES", "TRUE")]
  idcol <- names(e)[vapply(e, is_ach, logical(1))][1]
  if (is.na(idcol)) {
    prcol <- names(e)[vapply(e, function(v) is.character(v) &&
                               mean(grepl("^PR-", head(v, 500))) > 0.8, logical(1))][1]
    if (is.na(prcol) || is.null(profile_map_file))
      stop("Expression matrix has no ModelID (ACH-) and no cfg$profile_map_file to map ProfileIDs.")
    pm <- fread(profile_map_file)
    if ("ProfileType" %in% names(pm)) pm <- pm[tolower(ProfileType) == "rna"]
    e[, ModelID__ := pm$ModelID[match(e[[prcol]], pm$ProfileID)]]
    e <- e[!is.na(ModelID__)]
    idcol <- "ModelID__"
  }
  ids <- e[[idcol]]
  gene_cols <- grep("\\([0-9A-Za-z.]+\\)$", names(e), value = TRUE)
  if (length(gene_cols) < 100)
    gene_cols <- setdiff(names(e)[vapply(e, is.numeric, logical(1))], c("", "V1", "index"))
  m <- t(as.matrix(e[, ..gene_cols]))
  colnames(m) <- ids
  rownames(m) <- sub(" \\(.*\\)$", "", rownames(m))
  m <- m[!duplicated(rownames(m)), !duplicated(colnames(m)), drop = FALSE]
  storage.mode(m) <- "double"
  m
}

map_ids <- function(ids, annot) {
  ids <- trimws(as.character(ids))
  mid <- ifelse(!is.na(ids) & grepl("^ACH-[0-9]+$", ids), ids, NA_character_)
  key <- norm_name(ids)
  for (cc in c("CCLEName", "StrippedCellLineName", "CellLineName")) {
    todo <- is.na(mid) & !is.na(ids) & ids != ""
    if (!any(todo)) break
    ref <- norm_name(annot[[cc]]); ref[is.na(annot[[cc]])] <- NA
    mid[todo] <- annot$ModelID[match(key[todo], ref, incomparables = NA)]
  }
  mid
}

# Reads a sheet skipping title rows: the header is the first row with
# >= 2 non-empty cells and one of them containing "cluster"
read_sheet <- function(f, sh) {
  raw <- readxl::read_excel(f, sheet = sh, col_names = FALSE, col_types = "text",
                            n_max = 25, .name_repair = "minimal")
  if (!nrow(raw)) return(NULL)
  ncell  <- rowSums(!is.na(raw))
  has_cl <- apply(raw, 1, function(r) any(grepl("cluster", r, ignore.case = TRUE)))
  hdr <- which(ncell >= 2 & has_cl)[1]
  if (is.na(hdr)) hdr <- which(ncell >= 2)[1]
  if (is.na(hdr)) return(NULL)
  t <- as.data.table(readxl::read_excel(f, sheet = sh, skip = hdr - 1, col_types = "text",
                                        .name_repair = "unique"))
  t[, which(colSums(!is.na(t)) > 0), with = FALSE]
}

guess_id_col <- function(t, annot) {
  rates <- vapply(t, function(v) {
    v <- head(v[!is.na(v)], 300); if (!length(v)) 0 else mean(!is.na(map_ids(v, annot)))
  }, numeric(1))
  if (!length(rates) || max(rates) < 0.3) NA_character_ else names(t)[which.max(rates)]
}

cluster_cols <- function(t, idc) {
  others <- setdiff(names(t), idc)
  cc <- grep("cluster", others, ignore.case = TRUE, value = TRUE)
  if (!length(cc))   # no "cluster" in the name: columns with 2-30 distinct values
    cc <- others[vapply(t[, ..others], function(v) uniqueN(na.omit(v)) %in% 2:30, logical(1))]
  cc
}

pick <- function(x, prefer, avoid = "$^") {
  if (!length(x)) return(NA_character_)
  p <- grep(prefer, x, ignore.case = TRUE, value = TRUE); if (length(p)) return(p[1])
  a <- x[!grepl(avoid, x, ignore.case = TRUE)];          if (length(a)) return(a[1])
  x[1]
}

load_clusters_any <- function(cfg, annot) {
  f <- cfg$cluster_file
  if (grepl("\\.xlsx?$", f, ignore.case = TRUE)) {
    if (!requireNamespace("readxl", quietly = TRUE))
      stop("Package readxl is missing: install.packages('readxl')")
    sheets <- readxl::excel_sheets(f)
    tabs <- setNames(lapply(sheets, function(sh) suppressMessages(read_sheet(f, sh))), sheets)
  } else {
    tabs <- list(table = fread(f, colClasses = "character"))
  }
  BR <- "breast|brca|mama"
  info <- rbindlist(lapply(names(tabs), function(sh) {
    t <- tabs[[sh]]
    if (is.null(t) || !nrow(t))
      return(data.table(sheet = sh, n_rows = 0L, id_col = NA_character_, cluster_cols = ""))
    idc <- if (!is.null(cfg$cluster_id_col) && cfg$cluster_id_col %in% names(t)) cfg$cluster_id_col
    else guess_id_col(t, annot)
    cc  <- if (is.na(idc)) character(0) else cluster_cols(t, idc)
    data.table(sheet = sh, n_rows = nrow(t), id_col = idc, cluster_cols = paste(cc, collapse = " | "))
  }))
  msg("Sheets in the cluster file:")
  print(info, row.names = FALSE)
  valid <- info[!is.na(id_col) & cluster_cols != "", sheet]
  if (!length(valid)) stop("No sheet with a cell line column + a cluster column. ",
                           "Set cfg$cluster_sheet, cfg$cluster_id_col and cfg$cluster_col.")

  sh_pan <- cfg$cluster_sheet %||% pick(valid, "pan", BR)
  sh_br  <- cfg$cluster_sheet_breast %||% {
    b <- grep(BR, setdiff(valid, sh_pan), ignore.case = TRUE, value = TRUE); if (length(b)) b[1] else sh_pan }
  cols_of <- function(sh) strsplit(info[sheet == sh, cluster_cols], " | ", fixed = TRUE)[[1]]
  col_pan <- cfg$cluster_col %||% pick(cols_of(sh_pan), "pan", BR)
  col_br  <- cfg$cluster_col_breast %||% {
    cb <- grep(BR, cols_of(sh_br), ignore.case = TRUE, value = TRUE)
    if (length(cb)) cb[1] else if (sh_br == sh_pan) col_pan else pick(cols_of(sh_br), BR) }

  build <- function(sh, col, label) {
    t <- tabs[[sh]]; idc <- info[sheet == sh, id_col]
    if (!col %in% names(t)) stop("Column '", col, "' not found in sheet '", sh, "'")
    ids <- t[[idc]]; mid <- map_ids(ids, annot)
    out <- data.table(ModelID = mid, cluster = trimws(as.character(t[[col]])))
    lost <- !is.na(ids) & is.na(mid)
    out <- out[!is.na(ModelID) & !is.na(cluster) & cluster != "" & !duplicated(ModelID)]
    msg("Clusters %s: sheet '%s', id='%s', cluster='%s' -> %d lines, %d clusters (%s)%s",
        label, sh, idc, col, nrow(out), uniqueN(out$cluster),
        paste(cluster_order(out$cluster, cfg$cluster_colors), collapse = ", "),
        if (any(lost)) sprintf("\n  WARNING: %d not mapped to a ModelID (e.g. %s)", sum(lost),
                               paste(head(ids[lost], 5), collapse = ", ")) else "")
    out
  }
  pan <- build(sh_pan, col_pan, "pan-cancer")
  br  <- if (sh_br == sh_pan && col_br == col_pan) {
    msg("Clusters breast: same as pan-cancer, restricted to breast lines"); pan
  } else build(sh_br, col_br, "breast")
  list(pancancer = pan, breast = br)
}

load_crispr <- function(f, genes) {
  hdr <- read_header(f)
  idx <- vapply(genes, function(g) {
    i <- grep(paste0("^", g, "( \\(|$)"), hdr); if (length(i)) i[1] else NA_integer_ }, integer(1))
  if (anyNA(idx)) stop("Genes not found in ", basename(f), ": ", paste(genes[is.na(idx)], collapse = ", "))
  x <- fread(f, select = c(1L, unname(idx)))
  setnames(x, c("depmap_id", genes))
  if (!is_ach(x$depmap_id)) stop("The first column of ", basename(f), " does not contain ModelIDs (ACH-)")
  x[!duplicated(depmap_id)]
}

# Locates the DepMap input files in cfg$depmap_dir (unless they are set
# explicitly in cfg) and checks that every file needed exists.
#   need: any of "auc", "lfc", "expr", "crispr"
#         (Model.csv and the cluster file are always needed)
resolve_inputs <- function(cfg, need) {
  dd <- cfg$depmap_dir
  if (!dir.exists(dd)) stop("depmap_dir does not exist: ", dd)
  files <- list.files(dd, pattern = "\\.(csv|tsv)(\\.gz)?$", ignore.case = TRUE)

  if (is.null(cfg$model_file)) {
    mf <- grep("^Model\\.csv(\\.gz)?$", files, ignore.case = TRUE, value = TRUE)
    if (!length(mf)) stop("Model.csv not found in ", dd)
    cfg$model_file <- file.path(dd, mf[1])
  }

  if ("expr" %in% need) {
    if (is.null(cfg$expr_file)) {
      ex <- files[grepl("expression", files, ignore.case = TRUE) &
                    grepl("tpm", files, ignore.case = TRUE) &
                    grepl("logp1", files, ignore.case = TRUE) &
                    !grepl("transcript|isoform", files, ignore.case = TRUE)]
      if (!length(ex)) stop("No log-TPM expression file found in ", dd, ". Expression files: ",
                            paste(grep("expression", files, TRUE, value = TRUE), collapse = ", "))
      score <- 4 * (!grepl("profile", ex, ignore.case = TRUE)) +    # model level (ACH-)
        2 * grepl("proteincoding", ex, ignore.case = TRUE) +  # protein-coding genes
        1 * (!grepl("stranded", ex, ignore.case = TRUE))      # not the stranded subset
      ex <- ex[order(-score, ex)]
      cfg$expr_file <- file.path(dd, ex[1])
      if (length(ex) > 1) msg("  (other expression files: %s)", paste(ex[-1], collapse = ", "))
    }
    if (grepl("profile", basename(cfg$expr_file), ignore.case = TRUE) && is.null(cfg$profile_map_file)) {
      pm <- grep("^OmicsDefaultModelProfiles\\.csv$", files, ignore.case = TRUE, value = TRUE)
      if (length(pm)) cfg$profile_map_file <- file.path(dd, pm[1])
    }
  }

  if ("crispr" %in% need && is.null(cfg$crispr_file)) {
    cf <- grep("^CRISPRGeneEffect\\.csv(\\.gz)?$", files, ignore.case = TRUE, value = TRUE)
    if (!length(cf)) cf <- grep("CRISPR.*GeneEffect", files, ignore.case = TRUE, value = TRUE)
    cf <- cf[!grepl("uncorrected|raw", cf, ignore.case = TRUE)]
    if (!length(cf)) stop("CRISPRGeneEffect.csv not found in ", dd, "; set cfg$crispr_file")
    cfg$crispr_file <- file.path(dd, cf[1])
  }

  key <- c(auc = "auc_file", lfc = "lfc_file", expr = "expr_file", crispr = "crispr_file")
  req <- c(unname(key[need]), "model_file", "cluster_file")
  for (f in req) if (!has_file(cfg[[f]])) stop("File not found (", f, "): ", cfg[[f]])
  show <- unique(c(req, if (!is.null(cfg$profile_map_file)) "profile_map_file",
                   if (!is.null(cfg$lfc_file)) "lfc_file"))
  msg("Input files:")
  for (f in show) msg("  %-16s %s", f, if (has_file(cfg[[f]])) cfg[[f]] else "(not available)")
  cfg
}

## 2.5 PRISM statistics ----

assign_groups <- function(auc, cfg) {
  q <- quantile(auc, c(cfg$q_sens, cfg$q_res), na.rm = TRUE, names = FALSE)
  thr_s <- min(q[1], cfg$auc_sens_max)
  list(group = fifelse(auc <= thr_s, "sensitive", fifelse(auc >= q[2], "resistant", "intermediate")),
       q_low = q[1], q_high = q[2], thr_sens = thr_s, thr_res = q[2])
}

collapse_lineage <- function(lin, min_n) {
  lin <- as.character(lin); lin[is.na(lin) | lin == ""] <- "Unknown"
  tab <- table(lin); lin[lin %in% names(tab)[tab < min_n]] <- "Other"
  if (length(unique(lin)) < 2) NULL else factor(lin)
}

# Design ~ main (+ lineage if estimable and leaving >= 3 residual df)
make_design <- function(main, lineage, adjust, min_lineage_n) {
  if (adjust) {
    lf <- collapse_lineage(lineage, min_lineage_n)
    if (!is.null(lf)) {
      X <- model.matrix(~ main + lf)
      if (qr(X)$rank == ncol(X) && nrow(X) - ncol(X) >= 3)
        return(list(X = X, lineage_adjusted = TRUE))
    }
  }
  list(X = model.matrix(~ main), lineage_adjusted = FALSE)
}

fit_de <- function(E, X, coef, expr_type, min_expr_samples) {
  if (expr_type == "counts") {
    y <- edgeR::DGEList(E)
    y <- y[edgeR::filterByExpr(y, design = X), , keep.lib.sizes = FALSE]
    y <- edgeR::calcNormFactors(y)
    fit <- eBayes(lmFit(voom(y, X), X), robust = TRUE)
  } else {
    keep <- rowSums(E >= 1, na.rm = TRUE) >= min_expr_samples   # TPM >= 1
    fit <- eBayes(lmFit(E[keep, , drop = FALSE], X), trend = TRUE, robust = TRUE)
  }
  tt <- topTable(fit, coef = coef, number = Inf, sort.by = "P")
  out <- data.table(gene = rownames(tt), tt)
  # Rank = logFC * -log10(padj)  (padj floored to avoid Inf when it underflows to 0)
  out[, Rank := logFC * -log10(pmax(adj.P.Val, .Machine$double.xmin))]
  setcolorder(out, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "Rank"))
  out
}

# One limma fit with several numeric predictors (+ lineage); returns one table per coefficient
fit_multi <- function(E, vars, lineage, coefs, cfg) {
  X  <- model.matrix(~ ., data = vars)
  lf <- collapse_lineage(lineage, cfg$min_lineage_n)
  adj <- FALSE
  if (!is.null(lf)) {
    XL <- cbind(X, model.matrix(~ lf)[, -1, drop = FALSE])
    if (qr(XL)$rank == ncol(XL) && nrow(XL) - ncol(XL) >= 3) { X <- XL; adj <- TRUE }
  }
  keep <- rowSums(E >= 1, na.rm = TRUE) >= max(cfg$min_n_group, floor(0.1 * ncol(E)))
  fit  <- eBayes(lmFit(E[keep, , drop = FALSE], X), trend = TRUE, robust = TRUE)
  out <- lapply(coefs, function(cf) {
    tt <- topTable(fit, coef = cf, number = Inf, sort.by = "P")
    d  <- data.table(gene = rownames(tt), tt)
    d[, Rank := logFC * -log10(pmax(adj.P.Val, .Machine$double.xmin))]
    setcolorder(d, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "Rank"))
    d
  })
  list(res = setNames(out, coefs), lineage_adjusted = adj)
}

compare_models <- function(a, b, fdr) {
  m <- merge(a[, .(gene, t_drug = t, logFC_drug = logFC, padj_drug = adj.P.Val)],
             b[, .(gene, t_crispr = t, logFC_crispr = logFC, padj_crispr = adj.P.Val)], by = "gene")
  m[, category := fcase(padj_drug < fdr & padj_crispr < fdr & sign(t_drug) == sign(t_crispr), "Shared (same direction)",
                        padj_drug < fdr & padj_crispr < fdr, "Opposite direction",
                        padj_drug < fdr, "Drug only",
                        padj_crispr < fdr, "CRISPR only",
                        default = "NS")]
  m
}

classify_concordance <- function(m, fdr, lfc_thr) {
  m[, `:=`(sig_bin  = adj.P.Val_bin < fdr & abs(logFC_bin) >= lfc_thr,
           sig_cont = adj.P.Val_cont < fdr)]
  m[, category := fcase(sig_bin & sig_cont & sign(logFC_bin) == sign(slope_cont), "Both (concordant)",
                        sig_bin & sig_cont, "Both (discordant)",
                        sig_bin, "Binary only",
                        sig_cont, "Continuous only",
                        default = "NS")]
  m[, c("sig_bin", "sig_cont") := NULL]
  m
}

cor_stats <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y); n <- sum(ok)
  if (n < 4) return(data.table(n = n, rho = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method))
  r  <- unname(ct$estimate); se <- 1 / sqrt(max(n - 3, 1))
  data.table(n = n, rho = r, lo = tanh(atanh(r) - 1.96 * se), hi = tanh(atanh(r) + 1.96 * se), p = ct$p.value)
}

# Spearman = Pearson on ranks; optional lineage adjustment = partial correlation
# on ranks (lineage regressed out of the ranks). p-values: t approximation.
spearman_genes <- function(E, Y, lineage = NULL, cfg) {
  keep <- rowSums(E >= 1, na.rm = TRUE) >= max(3, floor(0.1 * ncol(E))) & apply(E, 1, sd, na.rm = TRUE) > 0
  E  <- E[keep, , drop = FALSE]
  R  <- t(apply(E, 1, rank, ties.method = "average"))
  RY <- apply(Y, 2, rank, ties.method = "average")
  X  <- matrix(1, ncol(E), 1)
  if (!is.null(lineage)) {
    lf <- as.character(lineage); lf[is.na(lf) | lf == ""] <- "Unknown"
    tab <- table(lf); lf[lf %in% names(tab)[tab < cfg$min_lineage_n]] <- "Other"
    if (length(unique(lf)) > 1) X <- model.matrix(~ factor(lf))
  }
  Q   <- qr.Q(qr(X))
  Rr  <- R  - (R %*% Q) %*% t(Q)                 # residual ranks (centred)
  RYr <- RY - Q %*% (t(Q) %*% RY)
  r   <- (Rr %*% RYr) / outer(sqrt(rowSums(Rr^2)), sqrt(colSums(RYr^2)))
  df  <- ncol(E) - 2 - (ncol(X) - 1)
  tt  <- r * sqrt(df / pmax(1 - r^2, 1e-12))
  p   <- 2 * pt(-abs(tt), df)
  out <- rbindlist(lapply(colnames(Y), function(v)
    data.table(gene = rownames(E), phenotype = v, rho = r[, v], p = p[, v])))
  out[, padj := p.adjust(p, "BH"), by = phenotype]
  if (cfg$flip) out[, rho := -rho]
  out[, `:=`(n = ncol(E), lineage_adjusted = ncol(X) > 1)]
  out
}


# One row per gene and CRISPR gene: drug rho/padj next to CRISPR rho/padj + category
pair_up <- function(cr, cfg) {
  drug <- cr[phenotype == "AUC", .(gene, rho_drug = rho, p_drug = p, padj_drug = padj)]
  rbindlist(lapply(cfg$genes, function(g0) {
    m <- merge(drug, cr[phenotype == g0, .(gene, rho_crispr = rho, p_crispr = p, padj_crispr = padj)], by = "gene")
    sd_ <- m$padj_drug < cfg$fdr & abs(m$rho_drug) >= cfg$rho_min
    sc_ <- m$padj_crispr < cfg$fdr & abs(m$rho_crispr) >= cfg$rho_min
    m[, category := fcase(sd_ & sc_ & sign(rho_drug) == sign(rho_crispr), "Shared (same direction)",
                          sd_ & sc_, "Opposite direction", sd_, "Drug only", sc_, "CRISPR only", default = "NS")]
    cbind(crispr_gene = g0, m)
  }))
}


## 2.6 PRISM plots ----
GROUP_COLS <- c(sensitive = "#2471A3", intermediate = "grey55", resistant = "#C0392B")
CONC_COLS <- c("Both (concordant)" = "#6A3D9A", "Binary only" = "#E69F00",
               "Continuous only" = "#009E73", "Both (discordant)" = "#D55E00", "NS" = "grey85")
COMP_COLS <- c("Shared (same direction)" = "#6A3D9A", "Drug only" = "#E69F00",
               "CRISPR only" = "#009E73", "Opposite direction" = "#D55E00", "NS" = "grey85")

plot_volcano <- function(res, file_base, fdr, lfc_thr, title,
                         lab_up = "Up in resistant", lab_down = "Up in sensitive",
                         xlab = "logFC (resistant - sensitive)", subtitle = NULL,
                         xlim = NULL, label = FALSE, label_n = 15, label_genes = NULL) {
  res <- copy(res)
  res[, panel := factor(panel, levels = unique(panel))]   # keep cluster order (ALL last)
  res[, sig := fifelse(adj.P.Val < fdr & logFC >=  lfc_thr & logFC > 0, lab_up,
                       fifelse(adj.P.Val < fdr & logFC <= -lfc_thr & logFC < 0, lab_down, "NS"))]

  # Symmetric x axis centred on 0; genes beyond a fixed limit are drawn at the edge as triangles
  m <- if (is.null(xlim)) 1.05 * max(abs(res$logFC), na.rm = TRUE) else xlim
  res[, `:=`(beyond = abs(logFC) > m, x = pmax(pmin(logFC, m), -m), y = -log10(P.Value))]
  res <- res[order(sig != "NS")]                         # coloured points on top

  cnt  <- res[, .(n_up = sum(sig == lab_up), n_down = sum(sig == lab_down)), by = panel]
  pthr <- res[adj.P.Val < fdr, .(y = -log10(max(P.Value))), by = panel]   # padj threshold line
  nf <- uniqueN(res$panel)
  cols <- setNames(c(GROUP_COLS[["resistant"]], GROUP_COLS[["sensitive"]], "grey75"), c(lab_up, lab_down, "NS"))

  g <- ggplot(res, aes(x, y, colour = sig)) +
    geom_point(aes(shape = beyond), size = if (nf == 1) 1 else 0.5, alpha = 0.75) +
    (if (lfc_thr > 0) geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = 2, linewidth = 0.3, colour = "grey40")
     else geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60")) +
    geom_hline(data = pthr, aes(yintercept = y), linetype = 2, linewidth = 0.3, colour = "grey40") +
    geom_text(data = cnt, aes(x = -Inf, y = Inf, label = n_down), inherit.aes = FALSE,
              hjust = -0.4, vjust = 1.4, colour = GROUP_COLS[["sensitive"]], fontface = "bold", size = 3.5) +
    geom_text(data = cnt, aes(x = Inf, y = Inf, label = n_up), inherit.aes = FALSE,
              hjust = 1.4, vjust = 1.4, colour = GROUP_COLS[["resistant"]], fontface = "bold", size = 3.5) +
    facet_wrap(~ panel, scales = "free_y") +
    scale_x_continuous(limits = c(-m, m), expand = expansion(mult = 0.03)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +      # headroom for the counts
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    scale_colour_manual(values = cols, breaks = names(cols), drop = FALSE) +
    labs(title = title,
         subtitle = subtitle %||% (if (lfc_thr > 0) sprintf("Coloured: padj < %.2g and |logFC| >= %.2g", fdr, lfc_thr)
                                   else sprintf("Coloured: padj < %.2g", fdr)),
         caption = if (any(res$beyond)) sprintf("Triangles: genes beyond the axis limit (|x| > %.2g)", m) else NULL,
         x = xlab, y = "-log10(P)", colour = NULL) +
    theme_bw(10) + theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size = 2.5)))

  if (label) {
    top <- res[sig != "NS"][order(adj.P.Val), head(.SD, label_n), by = .(panel, sig)]
    goi <- res[gene %in% label_genes]
    lab <- unique(rbind(top, goi, fill = TRUE), by = c("panel", "gene"))
    if (nrow(goi))   # genes of interest: black ring around the point
      g <- g + geom_point(data = goi, shape = 21, colour = "black", fill = NA, size = 2.2, stroke = 0.5)
    if (nrow(lab)) {
      g <- if (requireNamespace("ggrepel", quietly = TRUE))
        g + ggrepel::geom_text_repel(data = lab, aes(label = gene), colour = "black", size = 2.7,
                                     max.overlaps = Inf, box.padding = 0.4, min.segment.length = 0,
                                     segment.size = 0.2, segment.colour = "grey40", show.legend = FALSE)
      else
        g + geom_text(data = lab, aes(label = gene), colour = "black", size = 2.6, vjust = -0.7,
                      check_overlap = TRUE, show.legend = FALSE)
    }
  }

  if (nf == 1) save_plot(g, file_base, width = 6, height = 5.5)
  else         save_plot(g, file_base, width = 3.6 * min(nf, 3) + 1, height = 3.4 * ceiling(nf / 3) + 1.3)
}

plot_concordance <- function(m, file_base, title) {
  m <- copy(m)
  m[, panel := factor(panel, levels = unique(panel))]
  m <- m[order(category != "NS")]
  m[, category := factor(category, levels = names(CONC_COLS))]
  r <- m[, .(lab = sprintf("rho = %.2f", cor(t_bin, t_cont, method = "spearman"))), by = panel]
  nf <- uniqueN(m$panel)
  g <- ggplot(m, aes(t_bin, t_cont, colour = category)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_point(size = 0.5, alpha = 0.7) +
    geom_text(data = r, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
              hjust = -0.1, vjust = 1.5, size = 3, colour = "grey20") +
    facet_wrap(~ panel) +
    scale_colour_manual(values = CONC_COLS, drop = FALSE) +
    labs(title = title, subtitle = "Each point is a gene; t statistics of both models",
         x = "t, binary model (resistant vs sensitive)", y = "t, continuous model (AUC)", colour = NULL) +
    theme_bw(10) + theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size = 2.5)))
  if (nf == 1) save_plot(g, file_base, width = 5.5, height = 5.5)
  else         save_plot(g, file_base, width = 3.4 * min(nf, 3) + 1, height = 3.2 * ceiling(nf / 3) + 1.4)
}

plot_auc <- function(dd, thr, file_base, title, ylab, cl_levels) {
  lv <- c("pancancer", "breast")
  dd[, `:=`(context = factor(context, levels = lv), cluster = factor(cluster, levels = cl_levels),
            group = factor(group, levels = names(GROUP_COLS)))]
  thr[, context := factor(context, levels = lv)]
  g <- ggplot(dd, aes(cluster, auc, colour = group)) +
    geom_jitter(width = 0.2, height = 0, size = 1.3, alpha = 0.8) +
    geom_hline(data = thr, aes(yintercept = thr_sens), linetype = 2, colour = GROUP_COLS[["sensitive"]]) +
    geom_hline(data = thr, aes(yintercept = thr_res),  linetype = 2, colour = GROUP_COLS[["resistant"]]) +
    facet_wrap(~ context, scales = "free_x") +
    scale_colour_manual(values = GROUP_COLS) +
    labs(title = title, x = "Cluster", y = ylab, colour = NULL) +
    theme_bw(10) + theme(panel.grid.minor = element_blank())
  save_plot(g, file_base, width = 10, height = 4.5)
}

# Relative viability per dose by sensitivity group (median and IQR)
plot_dose_response_by_group <- function(points, groups, file_base, title) {
  d <- merge(points[, .(depmap_id, dose, v)], groups[, .(depmap_id, context, group)],
             by = "depmap_id", allow.cartesian = TRUE)
  s <- d[, .(med = median(v), lo = quantile(v, .25), hi = quantile(v, .75)), by = .(context, group, dose)]
  s[, `:=`(context = factor(context, levels = c("pancancer", "breast")),
           group = factor(group, levels = names(GROUP_COLS)))]
  g <- ggplot(s, aes(dose, med, colour = group, fill = group)) +
    geom_hline(yintercept = 1, linetype = 3, colour = "grey60") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line() + geom_point(size = 1) +
    scale_x_log10(breaks = sort(unique(s$dose)), labels = dose_labels) +
    facet_wrap(~ context) +
    scale_colour_manual(values = GROUP_COLS, aesthetics = c("colour", "fill")) +
    labs(title = title, x = "Dose (\u00b5M)", y = "Relative viability (2^LFC), median and IQR",
         colour = NULL, fill = NULL) +
    theme_bw(10) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(g, file_base, width = 9, height = 4.2)
}

# All clusters in the same plot (one figure per context) + cluster x group panels
plot_dose_response_by_cluster <- function(points, groups, cfg, out_dir) {
  d <- merge(points[, .(depmap_id, dose, v)],
             groups[!is.na(cluster), .(depmap_id, context, cluster, group)],
             by = "depmap_id", allow.cartesian = TRUE)
  if (!nrow(d)) return(invisible(NULL))
  d <- d[, .(v = mean(v)), by = .(context, depmap_id, cluster, group, dose)]   # one value per line and dose
  doses <- sort(unique(d$dose)); lo <- log10(min(doses)); hi <- log10(max(doses))
  dodge <- position_dodge(width = 0.12)
  stats_all <- list()

  for (ctx in c("pancancer", "breast")) {
    dc <- d[context == ctx]
    if (!nrow(dc)) next
    n <- dc[, .(n = uniqueN(depmap_id)), by = cluster]
    small <- n[n < cfg$min_n_curve, cluster]
    if (length(small))
      msg("  [%s] clusters not drawn (< %d lines): %s", ctx, cfg$min_n_curve, paste(small, collapse = ", "))
    dc <- dc[!cluster %in% small]
    if (!nrow(dc)) next
    lv   <- cluster_order(dc$cluster, cfg$cluster_colors)
    cols <- cluster_palette(lv, cfg$cluster_colors)
    leg  <- setNames(sprintf("%s (n=%d)", lv, n[match(lv, cluster), n]), lv)

    s <- dc[, .(n = .N, mean = mean(v), sem = sd(v) / sqrt(.N), median = median(v),
                q25 = quantile(v, .25), q75 = quantile(v, .75)), by = .(cluster, dose)]
    stats_all[[ctx]] <- cbind(context = ctx, s)
    s[, cluster := factor(cluster, levels = lv)]

    g <- ggplot(s, aes(dose, colour = cluster, group = cluster)) +
      geom_hline(yintercept = 1, linetype = 3, colour = "grey60")
    if (cfg$curve_style == "fit") {
      fits <- rbindlist(lapply(lv, function(k) {
        pk <- dc[cluster == k]
        f <- fit_curve4(log10(pk$dose), pk$v, lo, hi)
        if (is.null(f)) return(NULL)
        xg <- seq(lo, hi, length.out = 200)
        data.table(cluster = k, dose = 10^xg, v = curve4(f$par, xg))
      }))
      if (nrow(fits)) {
        fits[, cluster := factor(cluster, levels = lv)]
        g <- g + geom_line(data = fits, aes(y = v), linewidth = 0.8)
      }
      g <- g + geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.08,
                             linewidth = 0.35, position = dodge) +
        geom_point(aes(y = mean), size = 1.7, position = dodge)
      ylab <- "Relative viability (2^LFC), mean \u00b1 SEM"
    } else {
      g <- g + geom_line(aes(y = median), linewidth = 0.8) + geom_point(aes(y = median), size = 1.7)
      ylab <- "Relative viability (2^LFC), median"
    }
    g <- g + scale_x_log10(breaks = doses, labels = dose_labels) +
      scale_colour_manual(values = cols, labels = leg, name = "Cluster") +
      labs(title = sprintf("%s (%s) - %s", cfg$drug_name, cfg$screen_id,
                           if (ctx == "pancancer") "pan-cancer" else "breast"),
           x = "Dose (\u00b5M)", y = ylab) +
      theme_bw(11) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
    save_plot(g, file.path(out_dir, paste0("dose_response_by_cluster_", ctx)), width = 7, height = 4.8)
  }
  if (length(stats_all))
    fwrite(rbindlist(stats_all), file.path(out_dir, "dose_response_by_cluster.csv"))

  # Within each cluster: sensitive vs intermediate vs resistant
  s2 <- d[, .(med = median(v), lo = quantile(v, .25), hi = quantile(v, .75)),
          by = .(context, cluster, group, dose)]
  nn <- d[, .(S = uniqueN(depmap_id[group == "sensitive"]), R = uniqueN(depmap_id[group == "resistant"]),
              N = uniqueN(depmap_id)), by = .(context, cluster)]
  nn[, `:=`(context = factor(context, levels = c("pancancer", "breast")),
            cluster = factor(cluster, levels = cluster_order(nn$cluster, cfg$cluster_colors)))]
  setorder(nn, context, cluster)
  nn[, panel := sprintf("%s | %s (S=%d, R=%d, n=%d)", context, cluster, S, R, N)]
  s2 <- merge(s2, nn[, .(context = as.character(context), cluster = as.character(cluster), panel)],
              by = c("context", "cluster"))
  s2[, `:=`(panel = factor(panel, levels = nn$panel), group = factor(group, levels = names(GROUP_COLS)))]
  ncol_f <- min(uniqueN(nn$cluster), 6)
  g2 <- ggplot(s2, aes(dose, med, colour = group, fill = group)) +
    geom_hline(yintercept = 1, linetype = 3, colour = "grey60") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line() + geom_point(size = 0.9) +
    scale_x_log10(breaks = doses, labels = dose_labels) +
    facet_wrap(~ panel, ncol = ncol_f) +
    scale_colour_manual(values = GROUP_COLS, aesthetics = c("colour", "fill")) +
    labs(title = sprintf("%s: dose-response by cluster and group", cfg$drug_name),
         x = "Dose (\u00b5M)", y = "Relative viability (2^LFC), median and IQR", colour = NULL, fill = NULL) +
    theme_bw(9) + theme(legend.position = "bottom", panel.grid = element_blank(),
                        axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
  save_plot(g2, file.path(out_dir, "dose_response_by_cluster_and_group"),
           width = 2.5 * ncol_f + 1, height = 2.4 * ceiling(nrow(nn) / ncol_f) + 1.2)
}

# Lines to label: every line in small panels; otherwise the most dependent,
# the most sensitive and the largest outliers from the fit (+ cfg$label_lines)
pick_labels <- function(x, n, all_max, always) {
  if (nrow(x) <= all_max) return(x)
  x <- copy(x)
  x[, resid := auc - predict(lm(auc ~ gene_effect, data = x), x)]
  k <- max(1, ceiling(n / 3))
  idx <- unique(c(head(order(x$gene_effect), k), head(order(x$auc), k), head(order(-abs(x$resid)), k),
                  which(norm_name(x$name) %in% norm_name(always %||% character(0)))))
  x[idx]
}

coef_label <- function(st, method) {
  sym <- if (method == "spearman") "rho" else "r"
  st[, lab := sprintf("%s\n%s = %.2f, p = %s, n = %d", panel, sym, rho, formatC(p, format = "g", digits = 2), n)]
}

plot_corr <- function(dd, file_base, title, subtitle, cfg, fill_by = "cluster", ncol = 3,
                      label_n = cfg$label_n, label_all_max = cfg$label_all_max, width = 11, height = 4.2) {
  st <- coef_label(dd[, cor_stats(gene_effect, auc, cfg$cor_method), by = panel], cfg$cor_method)
  lv <- cluster_order(dd$cluster, cfg$cluster_colors)
  dd <- copy(dd)[, cluster := factor(cluster, levels = lv)]
  p_lev <- if (is.factor(dd$panel)) levels(droplevels(dd$panel)) else unique(dd$panel)
  st <- st[match(p_lev, as.character(panel))]
  dd[, strip := factor(st$lab[match(as.character(panel), as.character(st$panel))], levels = st$lab)]   # statistics in the strip
  base <- ggplot(dd, aes(gene_effect, auc)) +
    geom_vline(xintercept = cfg$dep_threshold, linetype = 2, colour = "grey60", linewidth = 0.3) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "black", fill = "grey75", linewidth = 0.6) +
    geom_point(aes(fill = cluster), shape = 21, colour = "grey25", stroke = 0.2, size = 1.9, alpha = 0.9) +
    facet_wrap(~ strip, scales = "free_x", ncol = ncol) +
    scale_fill_manual(values = cluster_palette(lv, cfg$cluster_colors), na.value = "grey80",
                      name = "Cluster", drop = TRUE) +
    labs(title = title, subtitle = subtitle,
         x = "CRISPR gene effect (Chronos; lower = more dependent)",
         y = sprintf("%s AUC (PRISM %s)", cfg$drug_name, cfg$screen_id)) +
    theme_bw(10) + theme(panel.grid = element_blank(), strip.text = element_text(size = 9))
  if (fill_by == "none") base <- base + guides(fill = "none")
  save_plot(base, file_base, width, height)

  lab <- dd[, pick_labels(.SD, label_n, label_all_max, cfg$label_lines), by = .(panel, strip)]
  g <- base +
    ggrepel::geom_text_repel(data = lab, aes(label = name), size = 2.3, colour = "black",
                             max.overlaps = Inf, min.segment.length = 0, segment.size = 0.2,
                             segment.colour = "grey40", box.padding = 0.3, seed = 1)
  save_plot(g, paste0(file_base, "_labeled"), width, height)
  st
}

axis_txt <- function(cfg) {
  if (!cfg$flip) list(x = "rho (expression vs AUC)   [+ = higher in resistant]",
                      y = "rho (expression vs CRISPR gene effect)   [+ = higher in non-dependent]")
  else           list(x = "rho (expression vs AUC) x -1   [+ = higher in sensitive]",
                      y = "rho (expression vs gene effect) x -1   [+ = higher in dependent]")
}

plot_rr <- function(pp, file_base, title, subtitle, cfg, ncol = 3, label_per_cat = cfg$label_per_category,
                    width = 12, height = 5.6) {
  pp <- copy(pp)
  sim <- pp[, .(r = cor(rho_drug, rho_crispr, method = "spearman"),
                n_shared = sum(category == "Shared (same direction)")), by = panel]
  sim[, strip := sprintf("%s\nsimilarity r = %.2f | shared = %d", panel, r, n_shared)]
  lev <- if (is.factor(pp$panel)) levels(droplevels(pp$panel)) else unique(pp$panel)
  sim <- sim[match(lev, as.character(panel))]
  pp[, strip := factor(sim$strip[match(as.character(panel), as.character(sim$panel))], levels = sim$strip)]
  pp[, category := factor(category, levels = names(COMP_COLS))]
  pp <- pp[order(category != "NS")]
  lim <- max(abs(c(pp$rho_drug, pp$rho_crispr)), na.rm = TRUE) * 1.05
  ax  <- axis_txt(cfg)
  base <- ggplot(pp, aes(rho_drug, rho_crispr, colour = category)) +
    geom_abline(slope = 1, intercept = 0, linetype = 3, colour = "grey55") +
    geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey65", linewidth = 0.3) +
    geom_point(size = 0.55, alpha = 0.75) +
    facet_wrap(~ strip, ncol = ncol) +
    coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    scale_colour_manual(values = COMP_COLS, drop = FALSE) +
    labs(title = title, subtitle = subtitle, x = ax$x, y = ax$y, colour = NULL) +
    theme_bw(10) + theme(legend.position = "bottom", panel.grid.minor = element_blank(),
                         strip.text = element_text(size = 8.5)) +
    guides(colour = guide_legend(override.aes = list(size = 2.5)))
  save_plot(base, file_base, width, height)

  top <- pp[category != "NS"][order(pmin(padj_drug, padj_crispr)), head(.SD, label_per_cat), by = .(strip, category)]
  goi <- pp[gene %in% cfg$label_genes]
  lab <- unique(rbind(top, goi, fill = TRUE), by = c("strip", "gene"))
  g <- base +
    geom_point(data = goi, shape = 21, colour = "black", fill = NA, size = 1.8, stroke = 0.45) +
    ggrepel::geom_text_repel(data = lab, aes(label = gene), colour = "black", size = 2.2,
                             max.overlaps = Inf, min.segment.length = 0, segment.size = 0.2,
                             segment.colour = "grey40", box.padding = 0.3, seed = 1, show.legend = FALSE)
  save_plot(g, paste0(file_base, "_labeled"), width, height)
  sim
}


# Analysis 1. Confounder-adjusted CDK-cyclin codependencies ---------------------
# Codependency = Pearson correlation of CRISPR gene effects across cell lines,
# for every pair of CDK and cyclin genes. The adjusted codependency is the
# correlation of the residuals of each gene effect after regressing out lineage,
# a proliferation score, the global dependency burden and, for each gene, its own
# expression and copy number (i.e. a partial correlation given these confounders).

run_adjusted_codependency <- function(depmap,
                                      out_dir = file.path(DIRS$out, "01_adjusted_codependency")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  ## 1.1 Inputs: cell lines present in every dataset, rows aligned ----
  copy_number_data <- load_copy_number(FILES$copy_number)
  common_lines <- Reduce(intersect, list(
    rownames(depmap$gene_effect),
    depmap$model$ModelID,
    rownames(depmap$expression),
    rownames(copy_number_data)
  ))

  dependency  <- depmap$gene_effect[common_lines, , drop = FALSE]
  expression  <- depmap$expression[common_lines, , drop = FALSE]
  copy_number <- copy_number_data[common_lines, , drop = FALSE]
  metadata    <- depmap$model |>
    filter(ModelID %in% common_lines) |>
    column_to_rownames("ModelID")
  metadata    <- metadata[common_lines, , drop = FALSE]
  # every confounder and gene effect below is combined by row position
  stopifnot(identical(rownames(dependency), rownames(expression)),
            identical(rownames(dependency), rownames(copy_number)),
            identical(rownames(dependency), rownames(metadata)))
  msg("Cell lines available: %d", nrow(dependency))

  ## 1.2 Confounders shared by all genes ----
  # (1) Lineage (one indicator column per lineage)
  lineage_dummy <- model.matrix(~ OncotreeLineage - 1, data = metadata)

  # (2) Proliferation score: mean gene effect of essential cell-cycle genes
  prolif_genes <- c("MKI67", "PCNA", "TOP2A", "MCM2", "CCNB1", "CDK1")
  prolif_genes <- intersect(prolif_genes, colnames(dependency))
  prolif_score <- rowMeans(dependency[, prolif_genes], na.rm = TRUE)

  # (3) Global dependency burden: mean gene effect across all genes
  global_burden <- rowMeans(dependency, na.rm = TRUE)

  # Expression and copy number are added for each gene in residualize_gene()
  confounders_base <- cbind(
    lineage_dummy,
    prolif_score  = prolif_score,
    global_burden = global_burden
  )
  confounders_base <- scale(confounders_base)
  msg("Base confounder matrix: %d x %d", nrow(confounders_base), ncol(confounders_base))

  ## 1.3 Residualise the gene effect of each gene ----
  # gene effect ~ lineage + proliferation + global burden + expression + copy number
  residualize_gene <- function(gene_name, gene_effects,
                               expression, copy_number,
                               confounders_base) {
    y <- gene_effects[, gene_name]

    # Gene-specific confounders: expression and copy number of the same gene
    if (gene_name %in% colnames(expression)) {
      expr_col <- scale(expression[, gene_name])
    } else {
      expr_col <- rep(0, nrow(confounders_base))
    }
    if (gene_name %in% colnames(copy_number)) {
      cn_col <- scale(copy_number[, gene_name])
    } else {
      cn_col <- rep(0, nrow(confounders_base))
    }
    X <- cbind(confounders_base,
               expr      = expr_col,
               cn        = cn_col)

    # Complete cases only
    valid <- complete.cases(X, y)
    X_ok  <- X[valid, ]
    y_ok  <- y[valid]

    fit <- lm(y_ok ~ X_ok)
    res <- numeric(length(y))
    res[valid]  <- residuals(fit)
    res[!valid] <- NA
    res
  }

  # CDK and cyclin genes (CCNQ, SPDYE3 and SPDYE1 are not included;
  # SPDYE1 has missing gene-effect values)
  cdk_cyclins <- c("CDK1", "CDK2", "CDK3", "CDK4", "CDK5", "CDK6", "CDK7", "CDK8", "CDK9", "CDK10",
                   "CDK12", "CDK13", "CDK14", "CDK15", "CDK16", "CDK17", "CDK18", "CDK19", "CDK20",
                   "CCNA1", "CCNA2", "CCNB1", "CCNB2", "CCNB3", "CCNC", "CCND1", "CCND2", "CCND3",
                   "CCNE1", "CCNE2", "CCNF", "CCNG1", "CCNG2", "CCNH", "CCNI", "CCNI2", "CCNJ",
                   "CCNJL", "CCNK", "CCNL1", "CCNL2", "CCNO", "CCNP", "CCNT1", "CCNT2", "CCNY",
                   "CCNYL1", "CNPPD1", "CNTD1", "CABLES1", "CABLES2", "CDK5R1", "CDK5R2",
                   "SPDYA", "SPDYE4")

  residuals_matrix <- sapply(cdk_cyclins, function(g) {
    residualize_gene(g, dependency, expression,
                     copy_number, confounders_base)
  })
  rownames(residuals_matrix) <- rownames(dependency)
  msg("Residual matrix: %d x %d", nrow(residuals_matrix), ncol(residuals_matrix))

  ## 1.4 Raw and adjusted codependencies ----
  # Raw correlation (as in the original manuscript)
  corr_raw <- cor(
    dependency[, cdk_cyclins],
    use    = "pairwise.complete.obs",
    method = "pearson"
  )
  # Adjusted correlation: correlation of the residuals (= partial correlation)
  corr_adjusted <- cor(
    residuals_matrix,
    use    = "pairwise.complete.obs",
    method = "pearson"
  )

  # Check: well-known CDK-cyclin pairs should remain high after the adjustment,
  # whereas spurious pairs should move towards 0
  known_pairs <- data.frame(
    geneA = c("CDK4",  "CDK6",  "CDK2",  "CDK7",  "CDK8"),
    geneB = c("CCND1", "CCND3", "CCNE1", "CCNH",  "CCNC")
  )
  cat("\n--- Known pairs: raw vs adjusted ---\n")
  for (i in 1:nrow(known_pairs)) {
    a <- known_pairs$geneA[i]
    b <- known_pairs$geneB[i]
    cat(sprintf("%s - %s: raw=%.3f | adjusted=%.3f\n",
                a, b,
                corr_raw[a, b],
                corr_adjusted[a, b]))
  }

  # P-value of a Pearson correlation
  cor_pvalue <- function(r, n) {
    t_stat <- r * sqrt((n - 2) / (1 - r^2))
    2 * pt(abs(t_stat), df = n - 2, lower.tail = FALSE)
  }

  # One row per gene pair (upper triangle of both matrices)
  n_genes <- length(cdk_cyclins)
  results <- data.frame()
  for (i in 1:(n_genes - 1)) {
    for (j in (i + 1):n_genes) {
      gA <- cdk_cyclins[i]
      gB <- cdk_cyclins[j]

      r_raw <- corr_raw[gA, gB]
      r_adj <- corr_adjusted[gA, gB]
      # Effective N: complete pairs of residuals
      n_eff <- sum(complete.cases(
        residuals_matrix[, c(gA, gB)]
      ))

      results <- rbind(results, data.frame(
        geneA      = gA,
        geneB      = gB,
        r_raw      = r_raw,
        r_adjusted = r_adj,
        delta_r    = r_adj - r_raw,
        p_raw      = cor_pvalue(r_raw, n_eff),
        p_adjusted = cor_pvalue(r_adj, n_eff)
      ))
    }
  }

  # Benjamini-Hochberg FDR
  results$fdr_raw      <- p.adjust(results$p_raw,      method = "BH")
  results$fdr_adjusted <- p.adjust(results$p_adjusted, method = "BH")

  # Classification of each pair
  results <- results |> mutate(
    status = case_when(
      fdr_raw < 0.05 & fdr_adjusted < 0.05 ~ "Robust",
      fdr_raw < 0.05 & fdr_adjusted >= 0.05 ~ "Artifact",
      fdr_raw >= 0.05 & fdr_adjusted < 0.05 ~ "New pair",
      TRUE ~ "No sign"
    )
  )
  results$diff <- results$r_raw - results$r_adjusted
  print(table(results$status))

  ## 1.5 Figures and table ----
  # Heatmaps, raw vs adjusted (colour scale clipped to [-0.5, 0.5])
  clip <- function(m, lim = 0.5) {
    m[m > lim] <- lim
    m[m < -lim] <- -lim
    m
  }
  heat_cols <- colorRampPalette(c("#185FA5", "white", "#993C1D"))(200)

  pdf(file.path(out_dir, "heatmap_codependency_raw_vs_adjusted.pdf"), width = 10, height = 5)
  par(mfrow = c(1, 2))
  corrplot(clip(corr_raw),
           method  = "color",
           type    = "upper",
           col     = heat_cols,
           tl.cex  = 0.7,
           tl.col  = "black",
           title   = "Codependency raw",
           col.lim = c(-0.5, 0.5),
           is.corr = FALSE,
           mar     = c(0, 0, 2, 0)
  )
  corrplot(clip(corr_adjusted),
           method  = "color",
           type    = "upper",
           col     = heat_cols,
           tl.cex  = 0.7,
           tl.col  = "black",
           title   = "Codependency adjusted",
           col.lim = c(-0.5, 0.5),
           is.corr = FALSE,
           mar     = c(0, 0, 2, 0)
  )
  dev.off()

  # Scatter plot, raw vs adjusted correlation (well-known pairs labelled)
  highlight_pairs <- c("CDK4 CCND1", "CDK6 CCND3", "CDK2 CCNE1",
                       "CDK7 CCNH",  "CDK8 CCNC",  "CDK9 CCNK")

  results <- results |>
    mutate(status = factor(status,
                           levels = c("Robust", "New pair", "Artifact", "No sign"),
                           labels = c("Robust", "New pair", "Artifact", "Not significant"))) |>
    arrange(desc(status))   # grey points are drawn first, below the rest

  highlight <- results |>
    filter(paste(geneA, geneB) %in% highlight_pairs | paste(geneB, geneA) %in% highlight_pairs)

  lims <- range(c(results$r_raw, results$r_adjusted), na.rm = TRUE)

  p <- ggplot(results, aes(x = r_raw, y = r_adjusted)) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_point(aes(colour = status), alpha = 0.7, size = 1.2) +
    geom_point(data = highlight, shape = 21, size = 3,
               stroke = 0.7, colour = "black", fill = NA) +
    geom_label_repel(
      data = highlight,
      aes(label = paste(geneA, geneB, sep = "–")),
      size = 2.3, fill = alpha("white", 0.85),
      min.segment.length = 0, segment.size = 0.3,
      box.padding = 0.4, max.overlaps = Inf, seed = 42
    ) +
    scale_colour_manual(
      values = c("Robust"          = "#0F6E56",
                 "New pair"        = "#534AB7",
                 "Artifact"        = "#993C1D",
                 "Not significant" = "grey75"),
      name = NULL
    ) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    coord_equal(xlim = lims, ylim = lims) +
    labs(x = "Raw correlation (r)",
         y = "Adjusted correlation (r)") +
    theme_classic(base_size = 8) +
    theme(
      legend.position  = "bottom",
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      axis.text        = element_text(colour = "black"),
      axis.line        = element_line(linewidth = 0.3),
      axis.ticks       = element_line(linewidth = 0.3)
    )

  ggsave(file.path(out_dir, "scatter_codependency_raw_vs_adjusted.pdf"), p,
         width = 89, height = 95, units = "mm", device = pdf_device())

  write.table(results, file.path(out_dir, "codependency_raw_vs_adjusted.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  msg("Results in: %s", normalizePath(out_dir))

  invisible(list(results = results, corr_raw = corr_raw, corr_adjusted = corr_adjusted,
                 residuals = residuals_matrix))
}


# Analysis 2. Validation of the cluster score classifier (CSC) in DepMap ---------
# Cell lines are assigned to the dependency clusters of Table S3 from their
# expression profile with the scheme used to project TCGA samples (cluster
# scores = mean expression of the cluster markers minus control genes, z-scored
# within the test set, argmax). Markers are re-derived in each 70% training split
# with the recipe of Table S9 and the classifier is evaluated on the held-out 30%
# (10 repeated stratified splits), against a Random Forest benchmark, single-gene
# baselines and, in breast cancer cell lines, PAM50. The section also tests the
# prediction of CDK/cyclin gene effects from the cluster scores and the overlap
# between the Table S9 markers and those derived in each training split.

run_classifier_validation <- function(depmap,
                                      out_dir = file.path(DIRS$out, "02_classifier_validation")) {
  fig_dir <- file.path(out_dir, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  save_fig <- function(p, name, width, height) save_fig_mm(p, file.path(fig_dir, name), width, height)
  set.seed(42)

  ## 2.0 Parameters ----
  MARKER_SOURCE <- "derive"   # "derive" (main: markers derived in each training split) | "tableS9" (reference; leakage)
  N_MARKERS     <- 100
  N_REPEATS     <- 10
  TRAIN_FRAC    <- 0.70
  N_PERM        <- 1000
  CTRL_SIZE     <- 50
  N_BINS        <- 25
  SCORE_SEED    <- 0

  # Marker selection: IDENTICAL to the procedure used to build Table S9
  #   - all genes of the expression matrix, no expression filter
  #   - limma, ~0 + cluster, contrast cluster - mean of the other 10 clusters
  #   - eBayes over the 11 contrasts jointly
  #   - padj = BH over the complete p-value matrix (genes x 11 contrasts)
  #   - score = logFC * -log10(padj); top 100 genes by score
  DERIVE_MIN_MEAN <- NA             # NA = no expression filter (as in S9)
  DERIVE_CONTRAST <- "mean_others"  # as in S9
  DERIVE_RANK     <- "score"        # "score" (as in S9) | "t" | "logFC"
  AUTO_MATCH_S9   <- FALSE          # TRUE: search the recipe that best reproduces S9 (diagnostic only)
  S9_MIN_REPRO    <- 95             # minimum % of S9 reproduced with all cell lines

  # Single-gene baselines and gene effects to predict
  SIMPLE  <- c("CCND1", "CCNE1", "CDK4", "CDK6", "RB1", "ESR1")
  TARGETS <- c("CDK4", "CDK6", "CDK2", "CCND1", "CCND3", "CCNE1")

  ## 2.1 Inputs ----
  expr_all    <- as.matrix(depmap$expression)
  gene_effect <- as.matrix(depmap$gene_effect)
  metadata    <- depmap$model
  cluster_labels <- read_excel(FILES$table_s3, sheet = 1)

  common_lines <- intersect(rownames(expr_all), cluster_labels$ModelID)

  labels_df <- cluster_labels |>
    filter(ModelID %in% common_lines) |>
    distinct(ModelID, .keep_all = TRUE) |>
    mutate(cluster = as.character(cluster))

  unexpected <- setdiff(unique(labels_df$cluster), CLUSTERS)
  if (length(unexpected) > 0) stop("Unexpected clusters in Table S3: ", paste(unexpected, collapse = ", "))

  common_lines <- labels_df$ModelID
  y <- factor(labels_df$cluster, levels = CLUSTERS)
  names(y) <- common_lines

  # Scoring matrix: all genes without NA
  expr_score <- expr_all[common_lines, colSums(is.na(expr_all[common_lines, ])) == 0]

  # Gene universe for marker derivation (DERIVE_MIN_MEAN). NA = every gene of the
  # matrix, including genes with some NA (accepted by lmFit), exactly as for S9.
  # Scoring always uses expr_score (no NA).
  make_universe <- function(min_mean) {
    X <- expr_all[common_lines, , drop = FALSE]
    if (is.na(min_mean)) X else X[, which(colMeans(X, na.rm = TRUE) > min_mean), drop = FALSE]
  }
  expr_sub <- make_universe(DERIVE_MIN_MEAN)

  cat("Cell lines:", nrow(expr_score), "| Genes (scoring):", ncol(expr_score),
      "| Genes (marker derivation):", ncol(expr_sub), "\n")
  print(table(y))
  if (min(table(y)) < 10) warning("Some clusters have < 10 cell lines; their test sets will be very small.")

  # Table S9 (always loaded: reproducibility control, section 2.3)
  s9_raw <- read_excel(FILES$table_s9, sheet = 1)
  stopifnot(all(CLUSTERS %in% colnames(s9_raw)))
  s9 <- setNames(lapply(CLUSTERS, function(cl)
    head(unique(na.omit(as.character(s9_raw[[cl]]))), N_MARKERS)), CLUSTERS)
  markers_s9 <- s9
  if (MARKER_SOURCE == "tableS9")
    warning("MARKER_SOURCE = 'tableS9': markers were derived with all cell lines -> optimistic metrics (leakage).")

  ## 2.2 Functions ----
  # --- Markers: each cluster vs the rest (limma), Table S9 recipe ---
  # contrast = "mean_others": cluster - mean of the other clusters (as in S9)
  # contrast = "pooled_rest": cluster vs all other lines pooled (diagnostic only)
  # rank_by  = "score" (logFC * -log10(padj), as in S9) | "t" | "logFC"
  # Only genes with logFC > 0. With "mean_others", padj = BH over the complete
  # genes x contrasts matrix, as for S9.
  derive_markers <- function(X, labels, n,
                             contrast = DERIVE_CONTRAST, rank_by = DERIVE_RANK) {
    labels <- droplevels(labels)
    cats   <- levels(labels)
    tX     <- t(X)

    if (contrast == "mean_others") {
      design <- model.matrix(~ 0 + labels)
      colnames(design) <- cats
      fit <- lmFit(tX, design)
      contr_str <- sapply(cats, function(ct) {
        others <- setdiff(cats, ct)
        paste0(ct, " - (", paste(others, collapse = "+"), ")/", length(others))
      })
      cm <- makeContrasts(contrasts = contr_str, levels = design)
      colnames(cm) <- cats
      fit2 <- eBayes(contrasts.fit(fit, cm))
      lfc  <- fit2$coefficients
      tst  <- fit2$t
      padj <- matrix(p.adjust(fit2$p.value, method = "BH"),
                     nrow = nrow(lfc), dimnames = dimnames(lfc))
    } else if (contrast == "pooled_rest") {
      res <- lapply(cats, function(cl) {
        grp <- factor(ifelse(labels == cl, "in", "out"), levels = c("out", "in"))
        f   <- eBayes(lmFit(tX, model.matrix(~ grp)))
        cbind(lfc = f$coefficients[, 2], t = f$t[, 2], p = f$p.value[, 2])
      })
      lfc  <- sapply(res, function(r) r[, "lfc"]); colnames(lfc) <- cats
      tst  <- sapply(res, function(r) r[, "t"]);   colnames(tst) <- cats
      padj <- matrix(p.adjust(sapply(res, function(r) r[, "p"]), method = "BH"),
                     nrow = nrow(lfc), dimnames = dimnames(lfc))
    } else stop("Unknown DERIVE_CONTRAST: ", contrast)

    stat <- switch(rank_by,
                   score = lfc * -log10(padj),
                   t     = tst,
                   logFC = lfc,
                   stop("Unknown DERIVE_RANK: ", rank_by))

    setNames(lapply(cats, function(cl) {
      ok <- which(lfc[, cl] > 0 & !is.na(stat[, cl]))
      g  <- rownames(lfc)[ok]
      head(g[order(-stat[ok, cl])], n)
    }), cats)
  }

  # % of the S9 markers of each cluster reproduced by a marker list
  repro_vs_s9 <- function(mk) {
    sapply(CLUSTERS, function(cl) 100 * length(intersect(s9[[cl]], mk[[cl]])) / length(s9[[cl]]))
  }

  # --- Control genes: re-implementation of scanpy.tl.score_genes ---
  make_ctrl <- function(X_ref, markers, ctrl_size = CTRL_SIZE, n_bins = N_BINS, seed = SCORE_SEED) {
    avg     <- colMeans(X_ref)
    n_items <- round(length(avg) / (n_bins - 1))
    cuts    <- setNames(rank(avg, ties.method = "min") %/% n_items, names(avg))
    lapply(markers, function(gs) {
      gs <- intersect(gs, names(cuts))
      set.seed(seed)
      unique(unlist(lapply(sort(unique(cuts[gs])), function(b) {
        pool <- setdiff(names(cuts)[cuts == b], gs)
        sample(pool, min(ctrl_size, length(pool)))
      })))
    })
  }

  score_with_ctrl <- function(X, markers, ctrl) {
    sapply(names(markers), function(cl) {
      g <- intersect(markers[[cl]], colnames(X))
      rowMeans(X[, g, drop = FALSE]) - rowMeans(X[, ctrl[[cl]], drop = FALSE])
    })
  }

  score_self <- function(X, markers) score_with_ctrl(X, markers, make_ctrl(X, markers))

  argmax_cluster <- function(Z) factor(colnames(Z)[max.col(Z, ties.method = "first")], levels = CLUSTERS)

  to_macro <- function(f) factor(ifelse(as.character(f) %in% F_SUB, "F", as.character(f)), levels = MACRO)

  # --- Metrics ---
  summ <- function(pred, truth, lev = CLUSTERS) {
    cm <- confusionMatrix(factor(pred, levels = lev), factor(truth, levels = lev))
    bacc <- cm$byClass[, "Balanced Accuracy"]
    names(bacc) <- sub("Class: ", "", names(bacc))
    list(cm = cm,
         acc = unname(cm$overall["Accuracy"]),
         nir = unname(cm$overall["AccuracyNull"]),
         acc_p = unname(cm$overall["AccuracyPValue"]),
         bacc = bacc[lev],
         mbacc = mean(bacc, na.rm = TRUE))
  }

  mean_bacc_fast <- function(pred, truth) {
    mean(sapply(levels(truth), function(cl) {
      p <- pred == cl; t <- truth == cl
      if (!any(t) || all(t)) return(NA_real_)
      (mean(p[t]) + mean(!p[!t])) / 2
    }), na.rm = TRUE)
  }

  perm_p <- function(pred, truth, B = N_PERM) {
    obs  <- mean_bacc_fast(pred, truth)
    null <- replicate(B, mean_bacc_fast(pred, sample(truth)))
    (sum(null >= obs) + 1) / (B + 1)
  }

  nearest_centroid <- function(x_tr, y_tr, x_te) {
    cen <- tapply(x_tr, y_tr, mean, na.rm = TRUE)
    cen <- cen[!is.na(cen)]
    factor(names(cen)[vapply(x_te, function(v) which.min(abs(cen - v)), integer(1))],
           levels = CLUSTERS)
  }

  ## 2.3 Reproducibility control: does derive_markers() reproduce Table S9 with all lines? ----
  # Otherwise the validation would assess a procedure different from the one used
  # to project the TCGA samples.
  s9_genes <- unique(unlist(s9))
  cat(sprintf("\nS9 genes present in the 25Q3 expression matrix: %.1f%%\n",
              100 * mean(s9_genes %in% colnames(expr_all))))
  cat(sprintf("S9 genes present in the matrix without NA (expr_score): %.1f%%\n",
              100 * mean(s9_genes %in% colnames(expr_score))))

  if (AUTO_MATCH_S9) {
    grid <- expand.grid(min_mean = c(NA, 0, 0.5, 1),
                        contrast = c("mean_others", "pooled_rest"),
                        rank_by  = c("score", "t", "logFC"),
                        stringsAsFactors = FALSE)
    message("Searching the configuration that reproduces Table S9 (", nrow(grid), " combinations)...")
    grid$repro_mean <- NA_real_; grid$repro_min <- NA_real_
    for (k in seq_len(nrow(grid))) {
      mk <- derive_markers(make_universe(grid$min_mean[k]), y, N_MARKERS,
                           contrast = grid$contrast[k], rank_by = grid$rank_by[k])
      r  <- repro_vs_s9(mk)
      grid$repro_mean[k] <- mean(r); grid$repro_min[k] <- min(r)
      message(sprintf("  min_mean=%s | %s | %s -> %.1f%% (min %.1f%%)",
                      grid$min_mean[k], grid$contrast[k], grid$rank_by[k],
                      grid$repro_mean[k], grid$repro_min[k]))
    }
    grid <- grid[order(-grid$repro_mean), ]
    print(grid, row.names = FALSE)
    write.csv(grid, file.path(out_dir, "TableS_S9_reproduction_grid.csv"), row.names = FALSE)

    DERIVE_MIN_MEAN <- grid$min_mean[1]
    DERIVE_CONTRAST <- grid$contrast[1]
    DERIVE_RANK     <- grid$rank_by[1]
    expr_sub <- make_universe(DERIVE_MIN_MEAN)
    message(sprintf("Selected configuration: min_mean = %s | contrast = %s | rank = %s",
                    DERIVE_MIN_MEAN, DERIVE_CONTRAST, DERIVE_RANK))
  }

  all_markers   <- derive_markers(expr_sub, y, N_MARKERS)[CLUSTERS]
  s9_repro      <- repro_vs_s9(all_markers)
  s9_repro_tbl  <- tibble(cluster = CLUSTERS, pct_s9_reproduced = s9_repro)
  print(s9_repro_tbl)
  cat(sprintf("\nReproduction of Table S9 with all cell lines: %.1f%% (range %.1f–%.1f%%)\n",
              mean(s9_repro), min(s9_repro), max(s9_repro)))
  write.csv(s9_repro_tbl, file.path(out_dir, "TableS_S9_reproduction_all_lines.csv"), row.names = FALSE)

  if (mean(s9_repro) < S9_MIN_REPRO) {
    warning(sprintf(paste0(
      "derive_markers reproduces only %.1f%% of Table S9 with all cell lines (< %d%%).\n",
      "Table S9 was NOT generated with this procedure. Check the DepMap release used for S9, ",
      "the method (e.g. scanpy rank_genes_groups / Wilcoxon), the cell lines included and the gene filter."),
      mean(s9_repro), S9_MIN_REPRO), immediate. = TRUE)
  }

  ## 2.4 Repeated stratified 70/30 splits ----
  run_split <- function(seed) {
    set.seed(seed)
    tr <- createDataPartition(y, p = TRAIN_FRAC, list = FALSE)[, 1]
    tr_lines <- common_lines[tr]
    te_lines <- common_lines[-tr]
    ytr <- y[tr]; yte <- y[-tr]

    mk <- if (MARKER_SOURCE == "derive") derive_markers(expr_sub[tr_lines, ], ytr, N_MARKERS) else markers_s9
    mk <- mk[CLUSTERS]

    Xtr <- expr_score[tr_lines, , drop = FALSE]
    Xte <- expr_score[te_lines, , drop = FALSE]

    # (1) Main: TCGA scheme (scores z-scored within the test set)
    S_te_self <- score_self(Xte, mk)
    Z_te      <- scale(S_te_self)
    pred_main <- argmax_cluster(Z_te)

    # (2) Secondary: scores anchored in the training set
    ctrl_tr   <- make_ctrl(Xtr, mk)
    S_tr_anch <- score_with_ctrl(Xtr, mk, ctrl_tr)
    S_te_anch <- score_with_ctrl(Xte, mk, ctrl_tr)
    sc_tr     <- scale(S_tr_anch)
    Z_te_anch <- scale(S_te_anch, center = attr(sc_tr, "scaled:center"),
                       scale = attr(sc_tr, "scaled:scale"))
    pred_anch <- argmax_cluster(Z_te_anch)

    # (3) Benchmark: Random Forest
    feats <- intersect(unique(unlist(mk)), colnames(Xtr))
    rf    <- randomForest(x = Xtr[, feats], y = ytr, ntree = 500)

    S_tr_self <- score_self(Xtr, mk)

    list(seed = seed, train = tr_lines, test = te_lines, y_test = yte, markers = mk,
         S_tr_self = S_tr_self, S_te_self = S_te_self, Z_te = Z_te,
         pred_main = pred_main, pred_anch = pred_anch,
         pred_rf = predict(rf, Xte[, feats]),
         prob_rf = predict(rf, Xte[, feats], type = "prob"))
  }

  message("Running ", N_REPEATS, " 70/30 splits (MARKER_SOURCE = ", MARKER_SOURCE, ")...")
  splits <- map(seq_len(N_REPEATS), function(i) {
    message("  split ", i)
    run_split(seed = 41 + i)
  })

  MODELS <- c(`CSC (TCGA scheme)` = "pred_main",
              `CSC (train-anchored)` = "pred_anch",
              `Random Forest` = "pred_rf")

  ## 2.5 Metrics: 11 clusters and macro-clusters (F1-F6 merged) ----
  metrics_tbl <- map_dfr(seq_along(splits), function(i) {
    sp <- splits[[i]]
    map_dfr(names(MODELS), function(m) {
      p <- sp[[MODELS[[m]]]]
      s <- summ(p, sp$y_test)
      tibble(repeat_id = i, level = "11 clusters", model = m,
             accuracy = s$acc, nir = s$nir, acc_p_vs_nir = s$acc_p,
             mean_bacc = s$mbacc, perm_p_bacc = perm_p(p, sp$y_test))
    })
  })

  metrics_macro_tbl <- map_dfr(seq_along(splits), function(i) {
    sp <- splits[[i]]
    yt <- to_macro(sp$y_test)
    map_dfr(names(MODELS), function(m) {
      p <- to_macro(sp[[MODELS[[m]]]])
      s <- summ(p, yt, lev = MACRO)
      tibble(repeat_id = i, level = "Macro-clusters (F merged)", model = m,
             accuracy = s$acc, nir = s$nir, acc_p_vs_nir = s$acc_p,
             mean_bacc = s$mbacc, perm_p_bacc = perm_p(p, yt))
    })
  })

  metrics_all <- bind_rows(metrics_tbl, metrics_macro_tbl)
  metrics_summary <- metrics_all |>
    group_by(level, model) |>
    summarise(accuracy_mean = mean(accuracy), accuracy_sd = sd(accuracy),
              nir = mean(nir),
              n_splits_acc_gt_nir = sum(acc_p_vs_nir < 0.05),
              mean_bacc_mean = mean(mean_bacc), mean_bacc_sd = sd(mean_bacc),
              perm_p_bacc_max = max(perm_p_bacc), .groups = "drop")
  print(metrics_summary)
  write.csv(metrics_all,     file.path(out_dir, "TableS_classifier_metrics_by_split.csv"), row.names = FALSE)
  write.csv(metrics_summary, file.path(out_dir, "TableS_classifier_metrics_summary.csv"),  row.names = FALSE)

  perclass_tbl <- map_dfr(seq_along(splits), function(i) {
    sp <- splits[[i]]
    map_dfr(names(MODELS), function(m) {
      cm <- summ(sp[[MODELS[[m]]]], sp$y_test)$cm
      bc <- as.data.frame(cm$byClass)
      tibble(repeat_id = i, model = m, cluster = sub("Class: ", "", rownames(bc)),
             bacc = bc$`Balanced Accuracy`, sensitivity = bc$Sensitivity,
             precision = bc$`Pos Pred Value`)
    })
  })
  perclass_summary <- perclass_tbl |>
    group_by(model, cluster) |>
    summarise(across(c(bacc, sensitivity, precision),
                     list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE))),
              .groups = "drop") |>
    mutate(cluster = factor(cluster, levels = CLUSTERS)) |>
    arrange(model, cluster)
  print(perclass_summary, n = Inf)
  write.csv(perclass_summary, file.path(out_dir, "TableS_classifier_perclass.csv"), row.names = FALSE)

  ## 2.6 Assignment margin (TCGA scheme) ----
  margin_tbl <- map_dfr(seq_along(splits), function(i) {
    sp  <- splits[[i]]
    srt <- t(apply(sp$Z_te, 1, sort, decreasing = TRUE))
    tibble(repeat_id = i, margin = srt[, 1] - srt[, 2],
           correct = sp$pred_main == sp$y_test,
           correct_macro = to_macro(sp$pred_main) == to_macro(sp$y_test))
  }) |>
    group_by(repeat_id) |>
    mutate(margin_tertile = factor(ntile(margin, 3), labels = c("Low", "Medium", "High"))) |>
    ungroup()

  margin_summary <- margin_tbl |>
    group_by(repeat_id, margin_tertile) |>
    summarise(acc = mean(correct), acc_macro = mean(correct_macro), .groups = "drop") |>
    group_by(margin_tertile) |>
    summarise(across(c(acc, acc_macro), list(mean = mean, sd = sd)), .groups = "drop")
  print(margin_summary)
  write.csv(margin_summary, file.path(out_dir, "TableS_accuracy_by_margin.csv"), row.names = FALSE)

  ## 2.7 Baselines: single genes (nearest centroid) ----
  baseline_tbl <- map_dfr(seq_along(splits), function(i) {
    sp <- splits[[i]]
    map_dfr(intersect(SIMPLE, colnames(expr_all)), function(g) {
      p <- nearest_centroid(expr_all[sp$train, g], y[sp$train], expr_all[sp$test, g])
      s <- summ(p, sp$y_test)
      tibble(repeat_id = i, model = g, accuracy = s$acc, mean_bacc = s$mbacc)
    })
  })

  comparison_summary <- bind_rows(
    metrics_tbl |> dplyr::select(repeat_id, model, accuracy, mean_bacc),
    baseline_tbl
  ) |>
    group_by(model) |>
    summarise(accuracy_mean = mean(accuracy), accuracy_sd = sd(accuracy),
              mean_bacc_mean = mean(mean_bacc), mean_bacc_sd = sd(mean_bacc),
              .groups = "drop") |>
    arrange(desc(mean_bacc_mean))
  print(comparison_summary)
  write.csv(comparison_summary, file.path(out_dir, "TableS_model_comparison.csv"), row.names = FALSE)

  ## 2.8 PAM50: breast cancer cell lines only ----
  data(pam50.robust, package = "genefu")   # PAM50 centroids used by molecular.subtyping()

  breast_lines <- intersect(common_lines,
                            metadata$ModelID[metadata$OncotreeLineage == "Breast"])
  cat("\nBreast cancer cell lines with an assigned cluster:", length(breast_lines), "\n")

  X_br <- expr_all[breast_lines, , drop = FALSE]
  X_br <- X_br[, colSums(is.na(X_br)) == 0, drop = FALSE]

  entrez <- suppressWarnings(AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db, keys = colnames(X_br), column = "ENTREZID",
    keytype = "SYMBOL", multiVals = "first"))
  ok <- !is.na(entrez)

  annot_pam50 <- data.frame(probe = colnames(X_br)[ok],
                            Gene.Symbol = colnames(X_br)[ok],
                            EntrezGene.ID = entrez[ok],
                            row.names = colnames(X_br)[ok],
                            stringsAsFactors = FALSE)

  X_br_c <- sweep(X_br[, ok], 2, apply(X_br[, ok], 2, median), "-")

  pam50_call <- molecular.subtyping(sbt.model = "pam50", data = X_br_c,
                                    annot = annot_pam50, do.mapping = TRUE)
  pam50 <- setNames(as.character(pam50_call$subtype), rownames(X_br_c))
  print(table(PAM50 = pam50, Cluster = y[names(pam50)]))

  pam50_tbl <- map_dfr(seq_along(splits), function(i) {
    sp <- splits[[i]]
    tr <- intersect(sp$train, breast_lines)
    te <- intersect(sp$test,  breast_lines)

    map_tbl <- tibble(subtype = pam50[tr], cluster = y[tr]) |>
      count(subtype, cluster) |>
      group_by(subtype) |>
      slice_max(n, n = 1, with_ties = FALSE) |>
      ungroup()

    p_pam <- factor(map_tbl$cluster[match(pam50[te], map_tbl$subtype)], levels = CLUSTERS)
    p_csc <- sp$pred_main[match(te, sp$test)]
    yt    <- y[te]
    keep  <- !is.na(p_pam)

    tibble(repeat_id = i, n_breast_test = sum(keep),
           acc_pam50  = mean(p_pam[keep] == yt[keep]),
           acc_csc    = mean(p_csc[keep] == yt[keep]),
           bacc_pam50 = mean_bacc_fast(p_pam[keep], droplevels(yt[keep])),
           bacc_csc   = mean_bacc_fast(p_csc[keep], droplevels(yt[keep])))
  })
  print(pam50_tbl |> summarise(across(-repeat_id, list(mean = mean, sd = sd), na.rm = TRUE)))
  write.csv(pam50_tbl, file.path(out_dir, "TableS_PAM50_vs_CSC_breast.csv"), row.names = FALSE)

  ## 2.9 Prediction of CDK/cyclin gene effects ----
  FEATURE_SETS <- c(
    list(`Cluster scores (CSC)` = CLUSTERS),
    setNames(as.list(SIMPLE), paste(SIMPLE, "expression")),
    list(`All simple markers` = SIMPLE,
         `CSC + simple markers` = c(CLUSTERS, SIMPLE))
  )

  ge_tbl <- map_dfr(seq_along(splits), function(i) {
    sp <- splits[[i]]
    tr <- intersect(sp$train, rownames(gene_effect))
    te <- intersect(sp$test,  rownames(gene_effect))
    S  <- rbind(sp$S_tr_self, sp$S_te_self)
    D  <- data.frame(S[c(tr, te), , drop = FALSE],
                     expr_all[c(tr, te), SIMPLE, drop = FALSE],
                     check.names = FALSE)

    map_dfr(intersect(TARGETS, colnames(gene_effect)), function(tg) {
      ytr <- gene_effect[tr, tg]; yte <- gene_effect[te, tg]
      ok_tr <- !is.na(ytr); ok_te <- !is.na(yte)

      map_dfr(names(FEATURE_SETS), function(fs) {
        cols <- FEATURE_SETS[[fs]]
        dtr  <- data.frame(resp = ytr[ok_tr], D[tr[ok_tr], cols, drop = FALSE], check.names = FALSE)
        fit  <- lm(resp ~ ., data = dtr)
        pr   <- predict(fit, newdata = D[te[ok_te], cols, drop = FALSE])
        yo   <- yte[ok_te]
        tibble(repeat_id = i, target = tg, features = fs,
               r2_oos   = 1 - sum((yo - pr)^2) / sum((yo - mean(ytr[ok_tr]))^2),
               spearman = cor(yo, pr, method = "spearman"))
      })
    })
  })

  ge_summary <- ge_tbl |>
    group_by(target, features) |>
    summarise(r2_mean = mean(r2_oos), r2_sd = sd(r2_oos),
              rho_mean = mean(spearman), rho_sd = sd(spearman), .groups = "drop") |>
    arrange(target, desc(r2_mean))
  print(ge_summary, n = Inf)
  write.csv(ge_summary, file.path(out_dir, "TableS_gene_effect_prediction.csv"), row.names = FALSE)

  ## 2.10 Figures (publication style) ----
  # Widths in mm at column width (W1, W15, W2); 5-7 pt sans-serif text; no titles
  # (they go in the legend); gene names in italics; vector PDF + 600-dpi PNG.
  # The metrics quoted in the legends are written to FigS_captions.txt.

  # Axis labels: gene names in italics
  GENES_ITAL <- unique(c(SIMPLE, TARGETS))
  gene_lab <- function(x) {
    x  <- as.character(x)
    g  <- sub(" expression$", "", x)
    esc <- function(s) gsub("'", "\\\\'", s)
    out <- ifelse(x %in% GENES_ITAL, sprintf("italic('%s')", x),
                  ifelse(g %in% GENES_ITAL & grepl(" expression$", x),
                         sprintf("italic('%s')~'expression'", g),
                         sprintf("'%s'", esc(x))))
    parse(text = out)
  }

  MODEL_COLS <- c(`CSC (TCGA scheme)` = COL_MAIN,
                  `CSC (train-anchored)` = COL_LIGHT,
                  `Random Forest` = COL_GREY)

  mean_sd <- function(v) data.frame(y = mean(v), ymin = mean(v) - sd(v), ymax = mean(v) + sd(v))

  main_lab <- "CSC (TCGA scheme)"
  ms  <- metrics_summary |> filter(level == "11 clusters", model == main_lab)
  msm <- metrics_summary |> filter(level == "Macro-clusters (F merged)", model == main_lab)
  captions <- list()

  # Confusion matrix plot
  plot_cm <- function(tab, lev) {
    df <- as.data.frame(tab)
    colnames(df) <- c("Predicted", "True", "Freq")
    df <- df |>
      group_by(True) |> mutate(prop = Freq / sum(Freq)) |> ungroup() |>
      mutate(Predicted = factor(Predicted, levels = lev),
             True      = factor(True, levels = rev(lev)),
             diag      = as.character(Predicted) == as.character(True))
    ggplot(df, aes(Predicted, True)) +
      geom_tile(aes(fill = prop), colour = "white", linewidth = 0.3) +
      geom_tile(data = dplyr::filter(df, diag), fill = NA, colour = "black", linewidth = 0.35) +
      geom_text(data = dplyr::filter(df, Freq > 0),
                aes(label = Freq, colour = prop > 0.5),
                size = pt2mm(BASE_PT - 1.5), family = FONT) +
      scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
      scale_fill_gradient(low = COL_HEAT_LO, high = COL_HEAT_HI, limits = c(0, 1),
                          labels = percent, name = "Row %", na.value = "white",
                          guide = guide_colourbar(barwidth = unit(2, "mm"),
                                                  barheight = unit(18, "mm"))) +
      scale_x_discrete(position = "top", expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      coord_equal() +
      labs(x = "Predicted cluster", y = "True cluster") +
      theme_heatmap()
  }

  # A. 11 x 11 confusion matrix (TCGA scheme), summed over splits
  cm_sum <- Reduce(`+`, lapply(splits, function(sp) summ(sp$pred_main, sp$y_test)$cm$table))
  p_A <- plot_cm(cm_sum, CLUSTERS) +
    ttl("Cluster score classifier (TCGA scheme) — held-out test sets")
  save_fig(p_A, "FigS_confusion_matrix_csc", W15, 105)
  captions$A <- sprintf(
    "Confusion matrix of the cluster score classifier (TCGA scheme) on held-out DepMap cell lines, summed over %d stratified 70/30 splits. Numbers are counts; fill is the row-normalised proportion; diagonal outlined. Accuracy %.2f ± %.2f (no-information rate %.2f); mean balanced accuracy %.2f ± %.2f (mean ± SD across splits).",
    N_REPEATS, ms$accuracy_mean, ms$accuracy_sd, ms$nir, ms$mean_bacc_mean, ms$mean_bacc_sd)

  # B. Random Forest calibration (all clusters in one panel)
  MIN_BIN_N <- 5   # bins with fewer lines are not shown (very noisy)

  calib <- map_dfr(splits, function(sp) {
    map_dfr(CLUSTERS, function(cl) tibble(cluster = cl, prob = sp$prob_rf[, cl],
                                          obs = as.integer(sp$y_test == cl)))
  }) |>
    mutate(cluster = factor(cluster, levels = CLUSTERS),
           bin = cut(prob, seq(0, 1, 0.1), include.lowest = TRUE)) |>
    group_by(cluster, bin) |>
    summarise(mean_pred = mean(prob), mean_obs = mean(obs), n = n(), .groups = "drop") |>
    filter(n >= MIN_BIN_N)

  p_B <- ggplot(calib, aes(mean_pred, mean_obs, colour = cluster, group = cluster)) +
    geom_abline(linetype = "dashed", colour = "grey60", linewidth = 0.25) +
    geom_line(linewidth = 0.4) +
    geom_point(aes(size = n), alpha = 0.9, stroke = 0) +
    scale_colour_manual(values = CLUSTER_COLORS, name = "Cluster") +
    scale_size(range = c(0.5, 2.5), trans = "log10", name = "N lines",
               breaks = c(10, 100, 1000)) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       expand = expansion(mult = 0.01)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       expand = expansion(mult = 0.01)) +
    coord_equal() +
    guides(colour = guide_legend(ncol = 2, override.aes = list(size = 1.8, linewidth = 0.6)),
           size   = guide_legend(override.aes = list(colour = "grey40"))) +
    labs(x = "Mean predicted probability", y = "Observed fraction") +
    ttl("Random Forest calibration — held-out test sets") +
    theme_paper(grid = "y") +
    theme(legend.key.size = unit(3, "mm"))
  save_fig(p_B, "FigS_calibration_rf", 65, 60)
  captions$B <- sprintf("Calibration of the Random Forest benchmark on held-out lines (all splits pooled). Predicted probabilities binned in 0.1 intervals (bins with < %d lines omitted); point size indicates the number of lines per bin (log scale); dashed line, perfect calibration.", MIN_BIN_N)

  # C. Balanced accuracy per cluster (CSC vs Random Forest)
  p_C <- perclass_summary |>
    filter(model %in% c(main_lab, "Random Forest")) |>
    mutate(model = factor(model, levels = c(main_lab, "Random Forest"))) |>
    ggplot(aes(cluster, bacc_mean, fill = model)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50", linewidth = 0.25) +
    geom_col(position = position_dodge(0.8), width = 0.75) +
    geom_errorbar(aes(ymin = bacc_mean - bacc_sd, ymax = bacc_mean + bacc_sd),
                  position = position_dodge(0.8), width = 0.3, linewidth = 0.25) +
    scale_fill_manual(values = MODEL_COLS) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(ylim = c(0.4, 1)) +
    labs(x = "Cluster", y = "Balanced accuracy", fill = NULL) +
    ttl("Per-cluster classification performance") +
    theme_paper(grid = "y") + theme(legend.position = "top")
  save_fig(p_C, "FigS_perclass_bacc", W15, 60)
  captions$C <- "Per-cluster balanced accuracy on held-out lines (mean ± SD across splits). Dashed line, chance (0.5)."

  # D. Comparison with single-gene markers
  p_D <- comparison_summary |>
    mutate(model = fct_reorder(model, mean_bacc_mean),
           highlight = as.character(model) %in% names(MODELS)) |>
    ggplot(aes(mean_bacc_mean, model, fill = highlight)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50", linewidth = 0.25) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(xmin = mean_bacc_mean - mean_bacc_sd,
                      xmax = mean_bacc_mean + mean_bacc_sd),
                  width = 0.3, linewidth = 0.25, orientation = "y") +
    scale_fill_manual(values = c(`TRUE` = COL_MAIN, `FALSE` = COL_GREY), guide = "none") +
    scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    scale_y_discrete(labels = gene_lab) +
    labs(x = "Mean balanced accuracy", y = NULL) +
    ttl("Cluster classification: CSC vs single markers") +
    theme_paper(grid = "x")
  save_fig(p_D, "FigS_model_comparison", W1, 70)
  captions$D <- "Mean balanced accuracy on held-out lines for the cluster score classifier (CSC), Random Forest and nearest-centroid classifiers based on single genes (mean ± SD across splits). Dashed line, chance."

  # E. Accuracy by assignment-margin tertile
  p_E <- margin_summary |>
    pivot_longer(-margin_tertile, names_to = c("level", ".value"),
                 names_pattern = "(acc|acc_macro)_(mean|sd)") |>
    mutate(level = factor(recode(level, acc = "11 clusters",
                                 acc_macro = "Macro-clusters (F merged)"),
                          levels = c("11 clusters", "Macro-clusters (F merged)"))) |>
    ggplot(aes(margin_tertile, mean, fill = level)) +
    geom_col(position = position_dodge(0.8), width = 0.75) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                  position = position_dodge(0.8), width = 0.3, linewidth = 0.25) +
    scale_fill_manual(values = c(COL_MAIN, COL_GREY)) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Assignment margin (tertile)", y = "Accuracy", fill = NULL) +
    ttl("Accuracy increases with assignment confidence") +
    theme_paper(grid = "y") + theme(legend.position = "top")
  save_fig(p_E, "FigS_accuracy_by_margin", 40, 60)
  captions$E <- "Held-out accuracy of the CSC (TCGA scheme) stratified by assignment margin (z-score of the top cluster minus the second), in within-split tertiles (mean ± SD across splits)."

  # F. Prediction of gene effects (out-of-sample R2), all genes in one panel
  #    Rows = predictors; columns = gene effect predicted; colour and number = mean
  #    out-of-sample R2. Black outline = best predictor for each gene.

  ge_plot <- ge_summary |>
    mutate(features = factor(features, levels = rev(names(FEATURE_SETS))),
           target   = factor(target, levels = intersect(TARGETS, unique(target)))) |>
    group_by(target) |>
    mutate(best = r2_mean == max(r2_mean)) |>
    ungroup()

  n_feat <- length(FEATURE_SETS)
  R2_MAX <- max(0.5, ceiling(max(ge_plot$r2_mean) * 10) / 10)

  p_F <- ggplot(ge_plot, aes(target, features)) +
    geom_tile(aes(fill = r2_mean), colour = "white", linewidth = 0.3) +
    geom_tile(data = dplyr::filter(ge_plot, best), fill = NA, colour = "black", linewidth = 0.45) +
    geom_text(aes(label = sprintf("%.2f", r2_mean), colour = r2_mean > 0.6 * R2_MAX),
              size = pt2mm(BASE_PT - 1.5), family = FONT) +
    # separators: CSC | single markers | combined models
    geom_hline(yintercept = c(2.5, n_feat - 0.5), colour = "grey30", linewidth = 0.3) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
    scale_fill_gradient(low = "#F1EFE8", high = COL_MAIN, limits = c(0, R2_MAX),
                        oob = scales::squish, name = expression(italic(R)^2),
                        guide = guide_colourbar(barwidth = unit(2, "mm"),
                                                barheight = unit(18, "mm"))) +
    scale_x_discrete(position = "top", labels = gene_lab, expand = c(0, 0)) +
    scale_y_discrete(labels = gene_lab, expand = c(0, 0)) +
    labs(x = "Gene effect predicted", y = NULL) +
    ttl("Prediction of CDK/cyclin gene effects (out-of-sample R²)") +
    theme_heatmap()
  save_fig(p_F, "FigS_gene_effect_prediction", 65, 55)

  captions$F <- "Out-of-sample R² (mean across 10 splits) of linear models predicting CRISPR gene effects of CDK and cyclin genes (columns) in held-out cell lines from cluster scores, the expression of individual markers, all simple markers combined, or cluster scores plus all simple markers (rows). The best predictor for each gene is outlined; negative values are shown as 0. SDs are given in Supplementary Table S_gene_effect_prediction."

  # F (landscape). Same as F with rows = gene effect predicted and columns = predictors
  # Short predictor labels (same order as FEATURE_SETS)
  short_feat <- function(x) {
    x <- as.character(x)
    x <- sub(" expression$", "", x)
    x <- sub("^Cluster scores \\(CSC\\)$", "CSC", x)
    x <- sub("^All simple markers$", "All markers", x)
    x <- sub("^CSC \\+ simple markers$", "CSC + markers", x)
    x
  }
  FEAT_LEVELS <- short_feat(names(FEATURE_SETS))   # CSC, CCND1 ... ESR1, All markers, CSC + markers

  ge_land <- ge_summary |>
    mutate(features = factor(short_feat(features), levels = FEAT_LEVELS),
           target   = factor(target, levels = intersect(TARGETS, unique(target)))) |>
    group_by(target) |>
    mutate(best = r2_mean == max(r2_mean)) |>
    ungroup()

  R2_MAX <- max(0.5, ceiling(max(ge_land$r2_mean) * 10) / 10)

  p_F_land <- ggplot(ge_land, aes(features, fct_rev(target))) +
    geom_tile(aes(fill = pmax(r2_mean, 0)), colour = "white", linewidth = 0.3) +
    geom_tile(data = dplyr::filter(ge_land, best), fill = NA, colour = "black", linewidth = 0.45) +
    geom_text(aes(label = sprintf("%.2f", r2_mean), colour = r2_mean > 0.6 * R2_MAX),
              size = pt2mm(BASE_PT - 1.5), family = FONT) +
    # separators: CSC | single markers | combined models
    geom_vline(xintercept = c(1.5, length(FEAT_LEVELS) - 1.5),
               colour = "grey30", linewidth = 0.3) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
    scale_fill_gradient(low = "#F1EFE8", high = COL_MAIN, limits = c(0, R2_MAX),
                        oob = scales::squish, name = expression(italic(R)^2),
                        guide = guide_colourbar(barwidth = unit(18, "mm"),
                                                barheight = unit(2, "mm"),
                                                title.position = "left",
                                                title.vjust = 1)) +
    scale_x_discrete(position = "top", labels = gene_lab, expand = c(0, 0)) +
    scale_y_discrete(labels = gene_lab, expand = c(0, 0)) +
    labs(x = "Predictor", y = "Gene effect predicted") +
    ttl("Prediction of CDK/cyclin gene effects (out-of-sample R²)") +
    theme_heatmap() +
    theme(legend.position = "bottom",
          legend.margin = ggplot2::margin(-2, 0, 0, 0, "mm"),
          axis.text.x.top = element_text(angle = 30, hjust = 0, vjust = 0))

  save_fig(p_F_land, "FigS_gene_effect_prediction_landscape", 70, 50)

  captions$F_landscape <- paste0(
    "Out-of-sample R² (mean across 10 training/test splits) of linear models predicting the CRISPR ",
    "gene effect of CDK and cyclin genes (rows) in held-out cell lines from the cluster scores (CSC), ",
    "the expression of individual markers, all markers combined, or cluster scores plus all markers ",
    "(columns). The best predictor of each gene effect is outlined in black; negative values are ",
    "shown as zero in the colour scale. Standard deviations are given in Supplementary Table SX.")

  # G. Performance in each split vs chance
  chance_df <- metrics_all |> group_by(level) |> summarise(nir = mean(nir), .groups = "drop")
  chance_lines <- bind_rows(
    chance_df |> transmute(level, metric = "Accuracy", chance = nir),
    chance_df |> transmute(level, metric = "Mean balanced accuracy", chance = 0.5)
  )

  p_G <- metrics_all |>
    dplyr::select(repeat_id, level, model, accuracy, mean_bacc) |>
    pivot_longer(c(accuracy, mean_bacc), names_to = "metric", values_to = "value") |>
    mutate(metric = recode(metric, accuracy = "Accuracy", mean_bacc = "Mean balanced accuracy"),
           model = factor(model, levels = names(MODEL_COLS))) |>
    ggplot(aes(model, value, colour = model)) +
    geom_hline(data = chance_lines, aes(yintercept = chance),
               linetype = "dashed", colour = "grey50", linewidth = 0.25) +
    geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.3, fill = NA) +
    geom_point(position = position_jitter(width = 0.12, seed = 1), size = 0.8, alpha = 0.8) +
    facet_grid(metric ~ level, scales = "free_y", switch = "y") +
    scale_colour_manual(values = MODEL_COLS, guide = "none") +
    labs(x = NULL, y = NULL) +
    ttl("Held-out performance across splits") +
    theme_paper(grid = "y") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          strip.placement = "outside")
  save_fig(p_G, "FigS_metrics_across_splits", 60, 100)
  captions$G <- sprintf("Held-out accuracy and mean balanced accuracy in each of the %d splits, at the 11-cluster and macro-cluster (F1–F6 merged) level. Dashed lines: no-information rate (accuracy) or 0.5 (balanced accuracy).", N_REPEATS)

  # H. Null distribution of the permutation test (split 1)
  sp1 <- splits[[1]]
  set.seed(1)
  null_vals <- replicate(N_PERM, mean_bacc_fast(sp1$pred_main, sample(sp1$y_test)))
  obs1  <- mean_bacc_fast(sp1$pred_main, sp1$y_test)
  p_obs <- (sum(null_vals >= obs1) + 1) / (N_PERM + 1)

  p_H <- ggplot(tibble(null = null_vals), aes(null)) +
    geom_histogram(bins = 40, fill = COL_GREY, colour = NA) +
    geom_vline(xintercept = obs1, colour = COL_MAIN, linewidth = 0.5) +
    annotate("text", x = obs1, y = Inf, hjust = 1.08, vjust = 1.5,
             label = sprintf("Observed = %.3f\np = %.3g", obs1, p_obs),
             size = pt2mm(BASE_PT - 1), family = FONT, colour = COL_MAIN, lineheight = 0.9) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Mean balanced accuracy (permuted labels)", y = "Count") +
    ttl("Permutation test (CSC, TCGA scheme; split 1)") +
    theme_paper()
  save_fig(p_H, "FigS_permutation_null", 60, 55)
  captions$H <- sprintf("Null distribution of mean balanced accuracy obtained by permuting the true cluster labels of the held-out lines (%d permutations, split 1). Vertical line, observed value (p = %.3g).", N_PERM, p_obs)

  # I. Sensitivity, precision and balanced accuracy per cluster
  prev_df <- tibble(cluster = factor(names(table(y)), levels = CLUSTERS),
                    name = "Prevalence",
                    value = as.numeric(prop.table(table(y))))

  p_I <- perclass_summary |>
    filter(model == main_lab) |>
    transmute(cluster, Sensitivity = sensitivity_mean, Precision = precision_mean,
              `Balanced accuracy` = bacc_mean) |>
    pivot_longer(-cluster) |>
    bind_rows(prev_df) |>
    mutate(name = factor(name, levels = c("Prevalence", "Precision", "Sensitivity",
                                          "Balanced accuracy"))) |>
    ggplot(aes(cluster, name, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", value), colour = value > 0.6),
              size = pt2mm(BASE_PT - 1.5), family = FONT) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
    scale_fill_gradient(low = "#F1EFE8", high = COL_MAIN, limits = c(0, 1), name = NULL,
                        guide = guide_colourbar(barwidth = unit(2, "mm"),
                                                barheight = unit(14, "mm"))) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    ttl("Per-cluster performance (CSC, TCGA scheme)") +
    theme_heatmap()
  save_fig(p_I, "FigS_perclass_heatmap", 80, 30)
  captions$I <- "Per-cluster sensitivity, precision and balanced accuracy of the CSC (TCGA scheme) on held-out lines (mean across splits). Cluster prevalence is shown as the precision expected by chance."

  # J. Macro-cluster confusion matrix (F1–F6 -> F)
  cm_macro <- Reduce(`+`, lapply(splits, function(sp)
    summ(to_macro(sp$pred_main), to_macro(sp$y_test), lev = MACRO)$cm$table))
  p_J <- plot_cm(cm_macro, MACRO) +
    ttl("CSC (TCGA scheme) — macro-clusters (F1–F6 merged)")
  save_fig(p_J, "FigS_confusion_matrix_macro", W1, 80)
  captions$J <- sprintf(
    "As in the 11-cluster confusion matrix, with F1–F6 merged into a single macro-cluster F (same predictions). Accuracy %.2f ± %.2f; mean balanced accuracy %.2f ± %.2f.",
    msm$accuracy_mean, msm$accuracy_sd, msm$mean_bacc_mean, msm$mean_bacc_sd)

  # K. Heatmap of z-scores (split 1)
  z_df <- as.data.frame(sp1$Z_te) |>
    rownames_to_column("ModelID") |>
    mutate(true = sp1$y_test) |>
    pivot_longer(-c(ModelID, true), names_to = "score_cluster", values_to = "z") |>
    mutate(z = pmax(pmin(z, 3), -3),
           score_cluster = factor(score_cluster, levels = rev(CLUSTERS)))
  line_order <- z_df |>
    filter(as.character(score_cluster) == as.character(true)) |>
    arrange(true, desc(z)) |>
    pull(ModelID)

  p_K <- ggplot(z_df, aes(factor(ModelID, levels = line_order), score_cluster, fill = z)) +
    geom_raster() +
    facet_grid(~ true, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(low = "#185FA5", mid = "white", high = "#B2182B",
                         limits = c(-3, 3), breaks = c(-3, 0, 3), name = "z-score",
                         guide = guide_colourbar(barwidth = unit(2, "mm"),
                                                 barheight = unit(15, "mm"))) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "Held-out cell lines (grouped by true cluster)", y = "Cluster score") +
    ttl("Cluster z-scores in held-out cell lines (split 1)") +
    theme_heatmap() +
    theme(axis.text.x = element_blank(), panel.spacing = unit(0.4, "mm"),
          strip.text = element_text(size = BASE_PT - 1, face = "bold"))
  save_fig(p_K, "FigS_zscore_heatmap", W2, 60)
  captions$K <- "Standardised cluster scores (z, truncated at ±3) of held-out lines in split 1, grouped by true cluster and ordered by the score of their own cluster."

  # L. Added value of the cluster scores (delta R2)
  best_single <- ge_tbl |>
    filter(features %in% paste(SIMPLE, "expression")) |>
    group_by(repeat_id, target) |>
    summarise(best_single = max(r2_oos), .groups = "drop")

  dr2 <- ge_tbl |>
    filter(features %in% c("Cluster scores (CSC)", "All simple markers", "CSC + simple markers")) |>
    dplyr::select(repeat_id, target, features, r2_oos) |>
    pivot_wider(names_from = features, values_from = r2_oos) |>
    left_join(best_single, by = c("repeat_id", "target")) |>
    transmute(repeat_id, target,
              `CSC alone vs best single marker` = `Cluster scores (CSC)` - best_single,
              `CSC added to all simple markers` = `CSC + simple markers` - `All simple markers`) |>
    pivot_longer(-c(repeat_id, target), names_to = "comparison", values_to = "delta_r2")
  write.csv(dr2, file.path(out_dir, "TableS_gene_effect_delta_r2.csv"), row.names = FALSE)

  p_L <- ggplot(dr2, aes(target, delta_r2, colour = comparison)) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.25) +
    geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.6, seed = 1),
               alpha = 0.45, size = 0.7, stroke = 0) +
    stat_summary(fun.data = mean_sd, geom = "pointrange",
                 position = position_dodge(0.6), size = 0.25, linewidth = 0.35) +
    scale_colour_manual(values = c(COL_GREY, COL_MAIN)) +
    scale_x_discrete(labels = gene_lab) +
    labs(x = NULL, y = expression(Delta * " out-of-sample " * italic(R)^2), colour = NULL) +
    ttl("Contribution of cluster scores to predicting CDK/cyclin gene effects") +
    theme_paper(grid = "y") +
    theme(legend.position = "top", legend.direction = "vertical")
  save_fig(p_L, "FigS_gene_effect_delta_r2", 30, 65)
  captions$L <- "Change in out-of-sample R² attributable to cluster scores: CSC alone minus the best single-gene model, and CSC plus all single markers minus single markers alone. Small points, individual splits; large points, mean ± SD."

  # M. CSC vs PAM50 in breast cancer cell lines (paired by split)
  p_M <- pam50_tbl |>
    dplyr::select(repeat_id, `PAM50` = bacc_pam50, `CSC` = bacc_csc) |>
    pivot_longer(-repeat_id, names_to = "model", values_to = "bacc") |>
    mutate(model = factor(model, levels = c("PAM50", "CSC"))) |>
    ggplot(aes(model, bacc, group = repeat_id)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50", linewidth = 0.25) +
    geom_line(colour = "grey75", linewidth = 0.3) +
    geom_point(aes(colour = model), size = 1.4) +
    scale_colour_manual(values = c(PAM50 = COL_GREY, CSC = COL_MAIN), guide = "none") +
    scale_y_continuous(limits = c(0.4, 1)) +
    labs(x = NULL, y = "Mean balanced accuracy\n(held-out breast lines)") +
    ttl("Breast cancer cell lines: CSC vs PAM50") +
    theme_paper(grid = "y")
  save_fig(p_M, "FigS_PAM50_vs_CSC", 55, 65)
  captions$M <- sprintf("Mean balanced accuracy for cluster assignment of held-out breast cancer lines using PAM50 subtype (majority mapping learned in training) or the CSC. Each line connects one split (n = %d–%d test lines). Dashed line, chance.",
                        min(pam50_tbl$n_breast_test), max(pam50_tbl$n_breast_test))

  # PAM50. PAM50 subtype composition per cluster (breast cancer cell lines)
  #   FigS_PAM50_by_cluster_colnorm: % of each subtype within each cluster (columns sum to 100)
  #   FigS_PAM50_by_cluster_rownorm: % of each subtype falling in each cluster (rows sum to 100)
  DROP_EMPTY <- FALSE   # TRUE = drop clusters/subtypes without breast cancer lines (e.g. B and E)
  PAM50_LEV  <- c("LumA", "LumB", "Her2", "Basal", "Normal")

  pam_counts <- tibble(ModelID = names(pam50), subtype = pam50, cluster = y[names(pam50)]) |>
    dplyr::filter(!is.na(cluster), !is.na(subtype)) |>
    mutate(subtype = factor(subtype, levels = intersect(PAM50_LEV, unique(subtype))),
           cluster = factor(cluster, levels = CLUSTERS)) |>
    dplyr::count(subtype, cluster, .drop = FALSE)

  if (DROP_EMPTY) {
    keep_cl  <- pam_counts |> group_by(cluster) |> summarise(N = sum(n)) |> dplyr::filter(N > 0) |> pull(cluster)
    keep_sub <- pam_counts |> group_by(subtype) |> summarise(N = sum(n)) |> dplyr::filter(N > 0) |> pull(subtype)
    pam_counts <- pam_counts |>
      dplyr::filter(cluster %in% keep_cl, subtype %in% keep_sub) |>
      mutate(cluster = droplevels(cluster), subtype = droplevels(subtype))
  }

  # axis labels with the number of lines
  n_cl  <- pam_counts |> group_by(cluster) |> summarise(N = sum(n), .groups = "drop")
  n_sub <- pam_counts |> group_by(subtype) |> summarise(N = sum(n), .groups = "drop")
  lab_cl  <- setNames(sprintf("%s\n(n = %d)", n_cl$cluster, n_cl$N), n_cl$cluster)
  lab_sub <- setNames(sprintf("%s (n = %d)", n_sub$subtype, n_sub$N), n_sub$subtype)

  # global association test (Fisher, simulated p-value: many small cells)
  tab_ps <- xtabs(n ~ subtype + cluster, data = pam_counts)
  tab_ps <- tab_ps[rowSums(tab_ps) > 0, colSums(tab_ps) > 0, drop = FALSE]
  set.seed(1)
  fisher_p <- fisher.test(tab_ps, simulate.p.value = TRUE, B = 1e5)$p.value

  # normalize = "cluster" (by columns) | "subtype" (by rows)
  plot_pam50 <- function(normalize) {
    df <- pam_counts |>
      group_by(across(all_of(normalize))) |>
      mutate(pct = if (sum(n) > 0) 100 * n / sum(n) else NA_real_) |>
      ungroup() |>
      mutate(lab = ifelse(n == 0, "0", sprintf("%.0f%%\n(%d)", pct, n)))
    list(
      data = df,
      plot = ggplot(df, aes(cluster, subtype, fill = pct)) +
        geom_tile(colour = "white", linewidth = 0.3) +
        geom_text(aes(label = lab, colour = !is.na(pct) & pct > 50),
                  size = pt2mm(BASE_PT - 2), family = FONT, lineheight = 0.85) +
        scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
        scale_fill_gradient(low = "#F1EFE8", high = COL_MAIN, limits = c(0, 100), na.value = "white",
                            name = if (normalize == "cluster") "% of\ncluster" else "% of\nsubtype",
                            guide = guide_colourbar(barwidth = unit(2, "mm"), barheight = unit(15, "mm"))) +
        scale_x_discrete(labels = lab_cl, position = "top", expand = c(0, 0)) +
        scale_y_discrete(labels = lab_sub, limits = rev, expand = c(0, 0)) +
        labs(x = "Dependency cluster", y = "PAM50 subtype") +
        ttl(if (normalize == "cluster") "PAM50 composition of each cluster (columns sum to 100%)"
            else "Distribution of each PAM50 subtype across clusters (rows sum to 100%)",
            sprintf("Breast cancer cell lines; Fisher's exact test (simulated) p = %.2g", fisher_p)) +
        theme_heatmap()
    )
  }

  pam_col <- plot_pam50("cluster")
  pam_row <- plot_pam50("subtype")
  save_fig(pam_col$plot, "FigS_PAM50_by_cluster_colnorm", W15, 55)
  save_fig(pam_row$plot, "FigS_PAM50_by_cluster_rownorm", W15, 55)

  write.csv(bind_rows(mutate(pam_col$data, normalized_by = "cluster (columns)"),
                      mutate(pam_row$data, normalized_by = "PAM50 subtype (rows)")),
            file.path(out_dir, "TableS_PAM50_by_cluster_breast.csv"), row.names = FALSE)

  # Combined panel (columns / rows)
  p_pam_both <- (pam_col$plot / pam_row$plot) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = BASE_PT + 2, face = "bold", family = FONT))
  save_fig(p_pam_both, "FigS_PAM50_by_cluster_both", W15, 110)

  captions$PAM <- sprintf(paste0(
    "PAM50 subtype of DepMap breast cancer cell lines (n = %d) across dependency clusters. ",
    "(a) Percentage of the cell lines of each cluster assigned to each PAM50 subtype (columns sum to 100%%). ",
    "(b) Percentage of the cell lines of each PAM50 subtype falling in each cluster (rows sum to 100%%). ",
    "Numbers in parentheses are cell-line counts; white cells, clusters without breast cancer cell lines. ",
    "Association between PAM50 subtype and cluster: Fisher's exact test with Monte Carlo simulation, P = %.2g."),
    sum(pam_counts$n), fisher_p)

  # N. Assignment margin of correct vs incorrect assignments
  p_N <- margin_tbl |>
    mutate(outcome = factor(ifelse(correct, "Correct", "Incorrect"),
                            levels = c("Correct", "Incorrect"))) |>
    ggplot(aes(margin, fill = outcome, colour = outcome)) +
    geom_density(alpha = 0.4, linewidth = 0.3) +
    scale_fill_manual(values = c(Correct = COL_MAIN, Incorrect = COL_GREY)) +
    scale_colour_manual(values = c(Correct = COL_MAIN, Incorrect = "grey45")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Assignment margin (z top − z second)", y = "Density", fill = NULL, colour = NULL) +
    ttl("Correct assignments have larger margins") +
    theme_paper() + theme(legend.position = c(0.8, 0.85))
  save_fig(p_N, "FigS_margin_density", 60, 55)
  captions$N <- "Distribution of the assignment margin for correct and incorrect held-out assignments (CSC, TCGA scheme; all splits pooled)."

  # Composite supplementary figure (double column)
  fig_supp <- (p_A + p_J + plot_layout(widths = c(1.6, 1))) /
    p_I /
    (p_G | p_N) /
    (p_L | p_M + plot_layout(widths = c(1.8, 1))) +
    plot_layout(heights = c(1.5, 0.45, 1.2, 0.9)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = BASE_PT + 2, face = "bold", family = FONT))
  save_fig(fig_supp, "FigS_classifier_validation_panel", W2, 240)

  ## 2.11 Final markers for TCGA (all lines, 11 clusters) ----
  # = all_markers (section 2.3); with ~100% reproduction they match Table S9.
  final_markers <- if (MARKER_SOURCE == "derive") all_markers else markers_s9
  final_markers <- final_markers[CLUSTERS]
  max_len  <- max(lengths(final_markers))
  final_df <- as.data.frame(lapply(final_markers, function(g) c(g, rep(NA, max_len - length(g)))),
                            check.names = FALSE)
  write.csv(final_df, file.path(out_dir, "TableS_final_cluster_markers_for_TCGA.csv"), row.names = FALSE)
  cat("\nFinal markers exported:", paste(names(final_markers), lengths(final_markers),
                                                sep = "=", collapse = ", "), "\n")

  ## 2.12 Table S9 markers vs markers derived in each training split ----
  #  - pct_of_s9: % of the S9 markers recovered in each split
  #  - jaccard, pct_top20, CAT (top-k), train-vs-train stability
  #  - "All lines" reference: same recipe with all lines (should be ~100%)
  # Training markers of each split (re-derived if MARKER_SOURCE = "tableS9")
  train_markers <- lapply(splits, function(sp) {
    if (MARKER_SOURCE == "derive") sp$markers[CLUSTERS]
    else derive_markers(expr_sub[sp$train, ], y[sp$train], N_MARKERS)[CLUSTERS]
  })

  overlap_stats <- function(a, b) {
    sh <- intersect(a, b)
    tibble(n_s9 = length(a), n_other = length(b), n_shared = length(sh),
           pct_of_s9 = 100 * length(sh) / length(a),
           jaccard   = length(sh) / length(union(a, b)),
           pct_top20 = 100 * length(intersect(head(a, 20), head(b, 20))) /
             min(20, length(a), length(b)))
  }

  # --- S9 vs each split ---
  overlap_tbl <- map_dfr(seq_along(train_markers), function(i)
    map_dfr(CLUSTERS, function(cl)
      bind_cols(tibble(repeat_id = i, cluster = cl),
                overlap_stats(s9[[cl]], train_markers[[i]][[cl]]))))

  # --- S9 vs all lines (method reference) ---
  ref_all_tbl <- map_dfr(CLUSTERS, function(cl)
    bind_cols(tibble(cluster = cl), overlap_stats(s9[[cl]], all_markers[[cl]])))

  # --- Train vs train stability ---
  pairs <- combn(seq_along(train_markers), 2, simplify = FALSE)
  stab_tbl <- map_dfr(pairs, function(pr)
    map_dfr(CLUSTERS, function(cl) {
      a <- train_markers[[pr[1]]][[cl]]; b <- train_markers[[pr[2]]][[cl]]
      tibble(cluster = cl, pct = 100 * length(intersect(a, b)) / min(length(a), length(b)))
    }))

  # --- Frequency of each S9 gene across splits ---
  s9_gene_freq <- map_dfr(CLUSTERS, function(cl) {
    tibble(cluster = cl, rank_s9 = seq_along(s9[[cl]]), gene = s9[[cl]],
           in_expr_universe = s9[[cl]] %in% colnames(expr_sub),
           n_splits = sapply(s9[[cl]], function(g)
             sum(sapply(train_markers, function(m) g %in% m[[cl]]))),
           pct_splits = 100 * n_splits / length(train_markers))
  })

  # --- Summary per cluster ---
  overlap_summary <- overlap_tbl |>
    group_by(cluster) |>
    summarise(n_s9 = dplyr::first(n_s9),
              pct_of_s9_mean = mean(pct_of_s9), pct_of_s9_sd = sd(pct_of_s9),
              pct_of_s9_min = min(pct_of_s9), pct_of_s9_max = max(pct_of_s9),
              jaccard_mean = mean(jaccard), pct_top20_mean = mean(pct_top20),
              .groups = "drop") |>
    left_join(ref_all_tbl |> dplyr::select(cluster, pct_all_lines = pct_of_s9), by = "cluster") |>
    left_join(stab_tbl |> group_by(cluster) |>
                summarise(pct_train_vs_train = mean(pct), .groups = "drop"), by = "cluster") |>
    left_join(s9_gene_freq |> group_by(cluster) |>
                summarise(pct_s9_in_universe = 100 * mean(in_expr_universe),
                          pct_s9_in_all_splits = 100 * mean(n_splits == length(train_markers)),
                          pct_s9_in_80pct_splits = 100 * mean(pct_splits >= 80),
                          .groups = "drop"), by = "cluster") |>
    mutate(cluster = factor(cluster, levels = CLUSTERS)) |>
    arrange(cluster)

  overall <- overlap_tbl |> summarise(pct_of_s9_mean = mean(pct_of_s9), pct_of_s9_sd = sd(pct_of_s9))
  print(overlap_summary, width = Inf)
  cat(sprintf("\nOverall: %.1f%% ± %.1f%% of the S9 markers are recovered in the training splits (mean ± SD, %d splits x %d clusters)\n",
              overall$pct_of_s9_mean, overall$pct_of_s9_sd, length(train_markers), length(CLUSTERS)))

  write.csv(overlap_tbl,     file.path(out_dir, "TableS_markers_S9_vs_train_by_split.csv"), row.names = FALSE)
  write.csv(overlap_summary, file.path(out_dir, "TableS_markers_S9_vs_train_summary.csv"),  row.names = FALSE)
  write.csv(s9_gene_freq,    file.path(out_dir, "TableS_markers_S9_gene_frequency.csv"),    row.names = FALSE)

  # Figures
  ref_pts <- bind_rows(
    ref_all_tbl |> transmute(cluster, value = pct_of_s9, type = "All lines vs S9"),
    stab_tbl |> group_by(cluster) |> summarise(value = mean(pct), .groups = "drop") |>
      mutate(type = "Train vs train (stability)")
  ) |> mutate(cluster = factor(cluster, levels = CLUSTERS))

  # O. % of the S9 markers recovered in each split
  p_O <- overlap_tbl |>
    mutate(cluster = factor(cluster, levels = CLUSTERS)) |>
    ggplot(aes(cluster, pct_of_s9)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.3, colour = COL_MAIN, fill = NA) +
    geom_point(position = position_jitter(width = 0.12, seed = 1),
               size = 0.8, alpha = 0.7, colour = COL_MAIN, stroke = 0) +
    geom_point(data = ref_pts, aes(y = value, shape = type),
               size = 1.6, fill = "white", colour = "black", stroke = 0.4) +
    scale_shape_manual(values = c(`All lines vs S9` = 23, `Train vs train (stability)` = 4),
                       name = NULL) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Cluster", y = "S9 markers recovered\nin training split (%)") +
    ttl("Overlap between Table S9 markers and training-split markers") +
    theme_paper(grid = "y") + theme(legend.position = "top")
  save_fig(p_O, "FigS_markers_S9_vs_train", 60, 65)

  # P. Recovery according to the rank of the gene in S9
  p_P <- s9_gene_freq |>
    mutate(cluster = factor(cluster, levels = CLUSTERS),
           rank_bin = cut(rank_s9, breaks = c(0, 20, 40, 60, 80, 100),
                          labels = c("1–20", "21–40", "41–60", "61–80", "81–100"))) |>
    group_by(cluster, rank_bin) |>
    summarise(pct = mean(pct_splits), .groups = "drop") |>
    ggplot(aes(rank_bin, cluster, fill = pct)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.0f", pct), colour = pct > 60),
              size = pt2mm(BASE_PT - 1.5), family = FONT) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
    scale_fill_gradient(low = "#F1EFE8", high = COL_MAIN, limits = c(0, 100), name = "%",
                        guide = guide_colourbar(barwidth = unit(2, "mm"), barheight = unit(15, "mm"))) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(limits = rev(CLUSTERS), expand = c(0, 0)) +
    labs(x = "Rank in Table S9", y = "Cluster") +
    ttl("Recovery of S9 markers across splits by rank") +
    theme_heatmap()
  save_fig(p_P, "FigS_markers_S9_recovery_by_rank", W1, 75)

  captions$O <- sprintf(
    "Percentage of the top %d Table S9 markers of each cluster recovered among the top %d markers derived (limma, cluster vs rest) in each 70%% training split. Boxes and small points: %d splits; diamonds: markers derived with all lines; crosses: mean pairwise overlap between training splits. Overall: %.1f%% ± %.1f%%.",
    N_MARKERS, N_MARKERS, length(train_markers), overall$pct_of_s9_mean, overall$pct_of_s9_sd)
  captions$P <- "Mean percentage of splits in which Table S9 markers are recovered, grouped by their rank in Table S9."

  # Q. Heatmap split x cluster (% of S9 recovered)
  p_Q <- overlap_tbl |>
    mutate(cluster = factor(cluster, levels = CLUSTERS),
           split   = factor(repeat_id, levels = rev(sort(unique(repeat_id))))) |>
    ggplot(aes(cluster, split, fill = pct_of_s9)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.0f", pct_of_s9), colour = pct_of_s9 > 60),
              size = pt2mm(BASE_PT - 1.5), family = FONT) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
    scale_fill_gradient(low = "#F1EFE8", high = COL_MAIN, limits = c(0, 100), name = "%",
                        guide = guide_colourbar(barwidth = unit(2, "mm"), barheight = unit(15, "mm"))) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "Cluster", y = "Split") +
    ttl("S9 markers recovered in each training split (%)") +
    theme_heatmap()
  save_fig(p_Q, "FigS_markers_S9_vs_train_heatmap", W15, 60)

  # R. Concordance at the top (CAT): overlap in the top-k markers, all clusters in one panel
  K_GRID <- c(5, seq(10, N_MARKERS, by = 5))
  cat_fun <- function(a, b, k) {
    kk <- min(k, length(a), length(b))
    100 * length(intersect(head(a, kk), head(b, kk))) / kk
  }
  cat_tbl <- map_dfr(seq_along(train_markers), function(i)
    map_dfr(CLUSTERS, function(cl)
      tibble(repeat_id = i, cluster = cl, k = K_GRID,
             cat = sapply(K_GRID, function(k) cat_fun(s9[[cl]], train_markers[[i]][[cl]], k)))))
  cat_all <- map_dfr(CLUSTERS, function(cl)
    tibble(cluster = cl, k = K_GRID,
           cat = sapply(K_GRID, function(k) cat_fun(s9[[cl]], all_markers[[cl]], k))))
  write.csv(cat_tbl, file.path(out_dir, "TableS_markers_S9_CAT_by_split.csv"), row.names = FALSE)
  SHOW_CI  <- TRUE   # band with the 95% confidence interval of the mean across splits
  CI_ALPHA <- 0.15

  cat_sum <- cat_tbl |>
    group_by(cluster, k) |>
    summarise(m = mean(cat), s = sd(cat), n = n(), .groups = "drop") |>
    mutate(cluster = factor(cluster, levels = CLUSTERS),
           half_ci = qt(0.975, df = n - 1) * s / sqrt(n),   # 95% CI (Student t)
           lo = pmax(m - half_ci, 0),
           hi = pmin(m + half_ci, 100))
  write.csv(cat_sum, file.path(out_dir, "TableS_markers_S9_CAT_summary_CI.csv"), row.names = FALSE)

  # "All lines" reference: with the S9 recipe it is 100% for every cluster -> a single line
  all_ref <- cat_all |> group_by(k) |> summarise(cat = mean(cat), .groups = "drop")

  p_R <- ggplot(cat_sum, aes(k, m, colour = cluster, group = cluster)) +
  { if (SHOW_CI) geom_ribbon(aes(ymin = lo, ymax = hi, fill = cluster),
                             colour = NA, alpha = CI_ALPHA) } +
    geom_line(data = all_ref, aes(k, cat, linetype = "All cell lines"),
              inherit.aes = FALSE, colour = "black", linewidth = 0.3) +
    geom_line(linewidth = 0.5) +
    scale_colour_manual(values = CLUSTER_COLORS, name = "Cluster") +
    scale_fill_manual(values = CLUSTER_COLORS, guide = "none") +
    scale_linetype_manual(values = c(`All cell lines` = "22"), name = NULL) +
    scale_x_continuous(breaks = c(5, 25, 50, 75, 100), expand = expansion(mult = 0.01)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                       expand = expansion(mult = c(0, 0.02))) +
    guides(colour   = guide_legend(ncol = 2, order = 1, override.aes = list(linewidth = 0.8)),
           linetype = guide_legend(order = 2)) +
    labs(x = "Top-k markers", y = "Overlap with Table S9 (%)") +
    ttl("Concordance at the top between S9 and training-derived markers") +
    theme_paper(grid = "y") +
    theme(legend.key.width = unit(5, "mm"))
  save_fig(p_R, "FigS_markers_S9_CAT", 75, 50)

  captions$R <- paste0(
    "Concordance at the top (CAT): percentage of shared genes between the top-k Table S9 markers ",
    "and the top-k markers derived in each training split, as a function of k. Lines, mean across 10 splits ",
    "(one per cluster); shaded bands, 95% confidence interval of the mean (t distribution). ",
    "Dashed line: markers derived with all cell lines (100% for every cluster).")

  # S. Overlap vs cluster size
  size_df <- tibble(cluster = names(table(y)), n_lines = as.integer(table(y)))
  s_df <- overlap_summary |>
    mutate(cluster = as.character(cluster)) |>
    left_join(size_df, by = "cluster")
  rho <- suppressWarnings(cor.test(s_df$n_lines, s_df$pct_of_s9_mean, method = "spearman", exact = FALSE))

  p_S <- ggplot(s_df, aes(n_lines, pct_of_s9_mean)) +
    geom_errorbar(aes(ymin = pct_of_s9_mean - pct_of_s9_sd, ymax = pct_of_s9_mean + pct_of_s9_sd),
                  width = 0, linewidth = 0.25, colour = "grey40") +
    geom_point(size = 1.4, colour = COL_MAIN) +
    geom_text(aes(label = cluster), vjust = -0.9, size = pt2mm(BASE_PT - 1), family = FONT) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.6,
             label = sprintf("Spearman rho = %.2f\np = %.2g", rho$estimate, rho$p.value),
             size = pt2mm(BASE_PT - 1), family = FONT, lineheight = 0.9) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0, 100)) +
    labs(x = "Cell lines in cluster (log scale)", y = "S9 markers recovered (%)") +
    ttl("Marker recovery vs cluster size") +
    theme_paper(grid = "y")
  save_fig(p_S, "FigS_markers_S9_vs_cluster_size", W1, 65)

  # T. Recovery class of each S9 marker (stacked bars)
  #    Mutually exclusive classes (sum to 100% per cluster): all splits (10/10) |
  #    >50% of the splits (6-9) | at least one (1-5) | none (0).
  #    BY_RANK = TRUE additionally splits by S9 rank bin (one panel per bin).

  BY_RANK <- FALSE

  nS   <- length(train_markers)
  half <- floor(nS / 2)
  REC_LEV <- c(sprintf("All splits (%d/%d)", nS, nS),
               sprintf(">50%% of splits (%d–%d/%d)", half + 1, nS - 1, nS),
               sprintf("At least one split (1–%d/%d)", half, nS),
               sprintf("No split (0/%d)", nS),
               "Not in expression matrix")
  REC_COLS <- setNames(c(COL_MAIN, COL_LIGHT, "#DAD7F3", COL_GREY, "grey88"), REC_LEV)

  rec_df <- s9_gene_freq |>
    mutate(cluster  = factor(cluster, levels = CLUSTERS),
           rank_bin = cut(rank_s9, breaks = c(0, 20, 40, 60, 80, 100),
                          labels = paste("Rank", c("1–20", "21–40", "41–60", "61–80", "81–100"))),
           class    = case_when(!in_expr_universe  ~ REC_LEV[5],
                                n_splits == nS     ~ REC_LEV[1],
                                n_splits >  nS / 2 ~ REC_LEV[2],
                                n_splits >= 1      ~ REC_LEV[3],
                                TRUE               ~ REC_LEV[4]),
           class    = factor(class, levels = REC_LEV))

  # drop empty classes from the legend (e.g. "Not in expression matrix" = 0)
  used_lev <- REC_LEV[REC_LEV %in% unique(as.character(rec_df$class))]

  group_vars <- if (BY_RANK) c("cluster", "rank_bin") else "cluster"
  rec_sum <- rec_df |>
    dplyr::count(across(all_of(group_vars)), class) |>
    group_by(across(all_of(group_vars))) |>
    mutate(pct = 100 * n / sum(n)) |>
    ungroup() |>
    mutate(class = factor(class, levels = used_lev))
  write.csv(rec_sum, file.path(out_dir, "TableS_markers_S9_recovery_class.csv"), row.names = FALSE)

  p_T <- ggplot(rec_sum, aes(cluster, pct, fill = class)) +
    geom_col(width = 0.75, position = position_stack(reverse = TRUE)) +
    { if (BY_RANK) facet_wrap(~ rank_bin, nrow = 1) } +
    scale_fill_manual(values = REC_COLS[used_lev], name = NULL) +
    scale_y_continuous(limits = c(0, 100.01), breaks = seq(0, 100, 25),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Cluster", y = "Table S9 markers (%)") +
    ttl("Recovery of Table S9 markers across training splits") +
    theme_paper() +
    theme(legend.key.size = unit(2.5, "mm"),
          axis.text.x = if (BY_RANK) element_text(angle = 90, vjust = 0.5, hjust = 1) else element_text())
  save_fig(p_T, "FigS_markers_S9_recovery_class", 75, 40)

  captions$T <- sprintf(paste0(
    "Classification of each Table S9 marker by the number of training splits (of %d) in which it is ",
    "recovered among the top %d training-derived markers: all splits, more than half, at least one ",
    "but no more than half, or none.%s Categories are mutually exclusive."),
    nS, N_MARKERS, if (BY_RANK) " Panels show bins of Table S9 rank." else "")

  # Composite panel
  fig_mk <- (p_O + p_S + plot_layout(widths = c(1.6, 1))) /
    p_R /
    (p_T + p_P + plot_layout(widths = c(1.6, 1))) /
    p_Q +
    plot_layout(heights = c(1, 1.1, 1, 0.9)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = BASE_PT + 2, face = "bold", family = FONT))
  save_fig(fig_mk, "FigS_markers_S9_vs_train_panel", W2, 250)

  captions$Q <- "Percentage of Table S9 markers recovered in each training split and cluster."
  captions$S <- sprintf("Mean recovery of Table S9 markers (± SD across splits) versus cluster size (Spearman rho = %.2f, p = %.2g).",
                        rho$estimate, rho$p.value)

  ## 2.13 Numbers quoted in the text ----
  txt <- list(
    X_repro_all_lines = sprintf("%.1f%% (range %.1f–%.1f%%)", mean(s9_repro), min(s9_repro), max(s9_repro)),
    Z_recovered_train = sprintf("%.1f%% ± %.1f%%", mean(overlap_tbl$pct_of_s9), sd(overlap_tbl$pct_of_s9)),
    Y_train_vs_train  = sprintf("%.1f%% (cluster range %.1f–%.1f%%)", mean(stab_tbl$pct),
                                min(tapply(stab_tbl$pct, stab_tbl$cluster, mean)),
                                max(tapply(stab_tbl$pct, stab_tbl$cluster, mean))),
    W_top20_recovery  = sprintf("%.1f%%", mean(s9_gene_freq$pct_splits[s9_gene_freq$rank_s9 <= 20])),
    config            = sprintf("min_mean = %s | contrast = %s | rank = %s",
                                DERIVE_MIN_MEAN, DERIVE_CONTRAST, DERIVE_RANK)
  )
  str(txt)
  writeLines(paste(names(txt), unlist(txt), sep = ": "), file.path(out_dir, "numbers_for_text.txt"))

  ## Figure legends ----
  writeLines(paste0("(", names(captions), ") ", unlist(captions), "\n"),
             file.path(fig_dir, "FigS_captions.txt"))
  msg("Results in: %s", normalizePath(out_dir))

  invisible(list(splits = splits, metrics = metrics_summary, perclass = perclass_summary,
                 gene_effect_prediction = ge_summary, final_markers = final_markers,
                 captions = captions))
}


# Analyses 3-6. Palbociclib (PRISM) ------------------------------------------------
# Shared settings of the four PRISM analyses. Each analysis adds its own
# parameters with modifyList(); NULL entries are detected automatically.
PRISM_CFG <- list(
  # PRISM Repurposing secondary screen (official AUC, as in the rest of the paper)
  auc_file   = FILES$prism_auc,
  drug_name  = "palbociclib",
  screen_id  = "MTS010",
  auc_cap    = 1,                     # AUC > 1 -> 1 (NULL = no capping)

  # DepMap 25Q3 (files detected automatically in depmap_dir if NULL)
  depmap_dir       = DIRS$depmap,
  expr_file        = NULL,            # log2(TPM+1), protein-coding genes
  model_file       = NULL,            # Model.csv
  crispr_file      = NULL,            # CRISPRGeneEffect.csv
  profile_map_file = NULL,            # only if expression is indexed by ProfileID (PR-)
  cell_info_file   = NULL,            # optional fallback for the lineage annotation

  # Clusters (Table S3); sheet and columns are detected automatically
  cluster_file         = FILES$table_s3,
  cluster_sheet        = NULL,        # sheet for pan-cancer
  cluster_sheet_breast = NULL,        # sheet with a breast-specific clustering, if any
  cluster_id_col       = NULL,        # cell line column (ACH- IDs or names)
  cluster_col          = NULL,        # cluster column for pan-cancer
  cluster_col_breast   = NULL,        # breast cluster column, if different
  cluster_colors       = CLUSTER_COLORS
)


# Analysis 3. Differential expression, resistant vs sensitive, within clusters ----
#  Contexts: pan-cancer and breast.
#  1) AUC (cfg$auc_source):
#     "official": AUC from secondary-screen-dose-response-curve-parameters.csv,
#                 screen MTS010 (same values used in the rest of the paper).
#     "lfc":      recomputed from the log-fold-changes (4-parameter curve,
#                 0.6 nM - 2.5 uM; the 10 uM dose is excluded).
#     AUC > 1 (fitted viability above DMSO) are NOT removed: they are the most
#     resistant lines. They are capped at 1 (cfg$auc_cap).
#  2) Cutoff (per context, computed over ALL cell lines of that context):
#       sensitive    : AUC <= P25  and  AUC <= 0.80
#       resistant    : AUC >= P75
#       intermediate : excluded from the binary contrast
#     plus limma with continuous AUC (no cutoff) as a control.
#     logFC > 0 = higher expression in RESISTANT lines.
#  3) Dose-response curves (relative viability = 2^LFC, screen MTS010):
#     by sensitivity group, by cluster (all clusters in one plot) and by
#     cluster x group.

CFG_PRISM_DE <- modifyList(PRISM_CFG, list(
  # --- PRISM ---
  auc_source     = "official",          # "official" (as in the paper) | "lfc"
  lfc_file       = FILES$prism_lfc,
  # optional (same release); the analysis works without them:
  treatment_file = FILES$prism_treatment,
  cell_info_file = FILES$prism_cell_info,
  drug_brd       = "BRD-K51313569",     # Broad ID of palbociclib (without batch suffix)
  max_dose_uM    = 3,                   # only for auc_source = "lfc" (Inf = include 10 uM)
  expr_type      = "logtpm",            # "logtpm" -> limma-trend | "counts" -> voom (needs edgeR)

  # --- cutoff ---
  q_sens       = 0.25,
  q_res        = 0.75,
  auc_sens_max = 0.80,
  min_gap      = 0.10,

  # --- analysis ---
  min_n_group    = 3,
  min_lineage_n  = 3,
  include_all    = TRUE,
  run_continuous = TRUE,
  fdr            = 0.05,
  lfc_thr        = 0.1,                 # minimum |logFC| to colour/count a gene as DE

  # --- volcano plots ---
  volcano_xlim            = NULL,       # binary: NULL = symmetric from the data, or e.g. 1.5 -> [-1.5, 1.5]
  volcano_xlim_continuous = NULL,       # continuous: idem, e.g. 0.6 -> [-0.6, 0.6]
  volcano_label_n         = 15,         # significant genes labelled per side (labelled plots)
  volcano_label_genes     = c("RB1", "CDKN2A", "CDKN2B", "CDKN2D", "CCNE1", "CCNE2",
                              "CCNG1", "CCNP", "ESD", "TACSTD2"),   # always labelled (ringed)
  control_genes  = c("RB1", "CDKN2A", "CCNE1", "CCND1", "CDK6", "ESR1"),

  # --- dose-response plots ---
  curve_style  = "fit",                 # "fit" = 4-parameter curve + mean +/- SEM | "median"
  min_n_curve  = 3,                     # clusters with fewer lines are not drawn

  out_dir      = file.path(DIRS$out, "03_palbociclib_DE_by_cluster")
))

run_prism_de_by_cluster <- function(cfg) {
  cfg <- resolve_inputs(cfg, need = c(if (cfg$auc_source == "official") "auc" else "lfc", "expr"))
  dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)

  annot <- load_annotation(cfg)
  A <- NULL
  if (has_file(cfg$lfc_file)) {
    msg("Reading %s log-fold-changes ...", cfg$drug_name)
    A <- compute_auc(cfg, annot)
  }
  if (cfg$auc_source == "official") {
    auc_tab <- load_official_auc(cfg)
    if (!is.null(A)) {
      auc_tab <- merge(auc_tab, A$auc[, .(depmap_id, auc_lfc = auc)], by = "depmap_id", all.x = TRUE)
      msg("  cor(official AUC, LFC-based AUC without 10 uM) = %.2f (Spearman)",
          cor(auc_tab$auc_official, auc_tab$auc_lfc, method = "spearman", use = "complete"))
    } else auc_tab[, auc_lfc := NA_real_]
    auc_tab[, auc_raw := auc_official]
  } else {
    auc_tab <- A$auc[, .(depmap_id, auc_raw = auc, auc_official = NA_real_, auc_lfc = auc)]
  }
  auc_tab[, auc := auc_raw]
  if (!is.null(cfg$auc_cap)) {
    n_cap <- sum(auc_tab$auc > cfg$auc_cap)
    auc_tab[, auc := pmin(auc_raw, cfg$auc_cap)]
    msg("  %d lines with AUC > %.2g (%.1f%%): capped at %.2g, not removed",
        n_cap, cfg$auc_cap, 100 * n_cap / nrow(auc_tab), cfg$auc_cap)
  }

  msg("Loading expression ...")
  E <- load_expr(cfg$expr_file, cfg$profile_map_file)
  msg("  %d genes x %d cell lines", nrow(E), ncol(E))
  if (nrow(E) < 1000 || ncol(E) < 50)
    stop("The expression matrix looks incomplete (", nrow(E), " genes x ", ncol(E),
         " lines). Check cfg$expr_file.")
  cls <- load_clusters_any(cfg, annot)

  d <- auc_tab[depmap_id %in% colnames(E)]
  d <- merge(d, annot[, .(depmap_id = ModelID, ccle_name = CCLEName, lineage, is_breast)],
             by = "depmap_id", all.x = TRUE)
  n_cl <- sum(d$depmap_id %in% cls$pancancer$ModelID)
  msg("%d lines with AUC + expression (%d with a cluster, %d breast)",
      nrow(d), n_cl, sum(d$is_breast %in% TRUE))
  if (nrow(d) < 20) stop("Very few lines shared between PRISM and expression: check the IDs.")
  if (n_cl < 2 * cfg$min_n_group) stop("Almost no PRISM line has a cluster: check the cluster file.")

  summary_rows <- list(); thr_rows <- list(); grp_rows <- list(); ctrl_rows <- list()
  auc_plot_data <- list(); auc_plot_thr <- list()

  for (ctx in c("pancancer", "breast")) {
    u <- if (ctx == "breast") d[is_breast %in% TRUE] else copy(d)
    u <- merge(u, cls[[ctx]][, .(depmap_id = ModelID, cluster)], by = "depmap_id", all.x = TRUE)
    if (nrow(u) < 2 * cfg$min_n_group) { msg("[%s] only %d lines: skipped", ctx, nrow(u)); next }
    gr <- assign_groups(u$auc, cfg)
    u[, group := gr$group]
    flags <- c(if (gr$q_high - gr$q_low < cfg$min_gap) "low_dynamic_range",
               if (gr$thr_sens < gr$q_low) "sensitive_limited_by_absolute_AUC")
    thr_rows[[ctx]] <- data.table(
      drug = cfg$drug_name, context = ctx, n_lines = nrow(u), P25 = gr$q_low, P75 = gr$q_high,
      thr_sensitive = gr$thr_sens, thr_resistant = gr$thr_res,
      n_sensitive = sum(u$group == "sensitive"), n_resistant = sum(u$group == "resistant"),
      flags = paste(flags, collapse = ";"))
    msg("\n[%s] n=%d  sensitive: AUC<=%.3f (n=%d)  resistant: AUC>=%.3f (n=%d) %s",
        ctx, nrow(u), gr$thr_sens, sum(u$group == "sensitive"), gr$thr_res,
        sum(u$group == "resistant"), if (length(flags)) paste0("[", paste(flags, collapse = ","), "]") else "")
    grp_rows[[ctx]] <- u[, .(context = ctx, depmap_id, ccle_name, lineage, cluster,
                             auc, auc_raw, auc_official, auc_lfc, group)]

    uc <- u[!is.na(cluster)]
    auc_plot_data[[ctx]] <- uc[, .(context = ctx, cluster, auc, group)]
    auc_plot_thr[[ctx]]  <- data.table(context = ctx, thr_sens = gr$thr_sens, thr_res = gr$thr_res)

    clusters <- cluster_order(uc$cluster, cfg$cluster_colors)
    if (cfg$include_all) clusters <- c(clusters, "ALL")
    cdir <- file.path(cfg$out_dir, ctx)
    dir.create(cdir, showWarnings = FALSE, recursive = TRUE)
    volc <- list(); volc_c <- list(); conc <- list()

    for (k in clusters) {
      s <- if (k == "ALL") uc else uc[cluster == k]
      b <- s[group != "intermediate"]
      ns <- sum(b$group == "sensitive"); nr <- sum(b$group == "resistant")
      row <- data.table(context = ctx, cluster = k, n_cluster = nrow(s), n_sensitive = ns, n_resistant = nr,
                        n_intermediate = sum(s$group == "intermediate"), lineage_adjusted = NA,
                        n_genes = NA_integer_, n_up_resistant = NA_integer_, n_up_sensitive = NA_integer_,
                        n_continuous_sig = NA_integer_, n_robust = NA_integer_, status = "ok")
      res <- NULL; resc <- NULL

      # --- binary: resistant vs sensitive ---
      if (ns >= cfg$min_n_group && nr >= cfg$min_n_group) {
        main <- factor(b$group, levels = c("sensitive", "resistant"))
        des  <- make_design(main, b$lineage, ctx == "pancancer", cfg$min_lineage_n)
        res  <- fit_de(E[, b$depmap_id, drop = FALSE], des$X, "mainresistant",
                       cfg$expr_type, min(ns, nr))
        fwrite(res, file.path(cdir, sprintf("DE_cluster_%s.csv", safe(k))))
        row[, `:=`(lineage_adjusted = des$lineage_adjusted, n_genes = nrow(res),
                   n_up_resistant = sum(res$adj.P.Val < cfg$fdr & res$logFC >=  cfg$lfc_thr),
                   n_up_sensitive = sum(res$adj.P.Val < cfg$fdr & res$logFC <= -cfg$lfc_thr))]
        if (ns < 5 || nr < 5) row[, status := "ok_low_n"]
        volc[[k]] <- res[, .(gene, logFC, P.Value, adj.P.Val, panel = sprintf("%s (S=%d, R=%d)", k, ns, nr))]
        ctrl_rows[[paste(ctx, k)]] <- res[gene %in% cfg$control_genes,
                                          .(context = ctx, cluster = k, gene, logFC, adj.P.Val)]
      } else {
        row[, status := sprintf("binary_skipped (S=%d, R=%d; minimum %d per group)", ns, nr, cfg$min_n_group)]
      }

      # --- continuous: AUC without cutoff (logFC per +0.1 AUC) ---
      if (cfg$run_continuous && nrow(s) >= 2 * cfg$min_n_group + 1) {
        desc <- make_design(s$auc * 10, s$lineage, ctx == "pancancer", cfg$min_lineage_n)
        resc <- fit_de(E[, s$depmap_id, drop = FALSE], desc$X, "main", cfg$expr_type,
                       max(cfg$min_n_group, floor(0.1 * nrow(s))))
        fwrite(resc, file.path(cdir, sprintf("DE_continuous_cluster_%s.csv", safe(k))))
        row[, n_continuous_sig := sum(resc$adj.P.Val < cfg$fdr)]
        volc_c[[k]] <- resc[, .(gene, logFC, P.Value, adj.P.Val, panel = sprintf("%s (n=%d)", k, nrow(s)))]
      }

      # --- concordance binary vs continuous ---
      if (!is.null(res) && !is.null(resc)) {
        m <- merge(res[,  .(gene, logFC_bin = logFC, t_bin = t, adj.P.Val_bin = adj.P.Val)],
                   resc[, .(gene, slope_cont = logFC, t_cont = t, adj.P.Val_cont = adj.P.Val)], by = "gene")
        m <- classify_concordance(m, cfg$fdr, cfg$lfc_thr)
        row[, n_robust := sum(m$category == "Both (concordant)")]
        conc[[k]] <- cbind(context = ctx, cluster = k, panel = sprintf("%s (S=%d, R=%d)", k, ns, nr), m)
      }
      summary_rows[[paste(ctx, k)]] <- row
      msg("  cluster %-6s S=%-3d R=%-3d -> %s", k, ns, nr,
          if (is.na(row$n_genes)) row$status else
            sprintf("%d up-res / %d up-sens (padj<%.2g, |logFC|>=%.2g) | continuous: %s | robust: %s%s",
                    row$n_up_resistant, row$n_up_sensitive, cfg$fdr, cfg$lfc_thr,
                    row$n_continuous_sig, row$n_robust,
                    if (isTRUE(row$lineage_adjusted)) " [lineage-adjusted]" else ""))
    }
    # Per-cluster panels and the ALL comparison are saved as separate figures
    split_plot <- function(lst, base, FUN, ...) {
      per <- lst[names(lst) != "ALL"]
      if (length(per)) FUN(rbindlist(per), file.path(cdir, paste0(base, "_clusters")), ...)
      if ("ALL" %in% names(lst)) FUN(lst[["ALL"]], file.path(cdir, paste0(base, "_ALL")), ...)
    }
    # Volcanos: clusters panel + ALL without labels + ALL with labels
    volcano_set <- function(lst, base, ...) {
      split_plot(lst, base, plot_volcano, ...)
      if ("ALL" %in% names(lst))
        plot_volcano(lst[["ALL"]], file.path(cdir, paste0(base, "_ALL_labeled")), ..., label = TRUE,
                     label_n = cfg$volcano_label_n, label_genes = cfg$volcano_label_genes)
    }
    if (length(volc))
      volcano_set(volc, "volcano", fdr = cfg$fdr, lfc_thr = cfg$lfc_thr, xlim = cfg$volcano_xlim,
                  title = sprintf("%s - %s: resistant vs sensitive (binary model)", cfg$drug_name, ctx))
    if (length(volc_c))
      volcano_set(volc_c, "volcano_continuous", fdr = cfg$fdr, lfc_thr = 0, xlim = cfg$volcano_xlim_continuous,
                  title = sprintf("%s - %s: association with AUC (continuous model)", cfg$drug_name, ctx),
                  lab_up = "Higher with resistance", lab_down = "Higher with sensitivity",
                  xlab = "Change in log2(TPM+1) per +0.1 AUC",
                  subtitle = sprintf("expression ~ AUC%s; coloured: padj < %.2g",
                                     if (ctx == "pancancer") " + lineage" else "", cfg$fdr))
    if (length(conc)) {
      cc <- rbindlist(conc)
      fwrite(cc[, !"panel"][order(cluster, category != "Both (concordant)", adj.P.Val_bin)],
             file.path(cdir, "DE_concordance.csv"))
      split_plot(conc, "concordance_binary_vs_continuous", plot_concordance,
                 title = sprintf("%s - %s: binary vs continuous model", cfg$drug_name, ctx))
    }
  }

  grp <- rbindlist(grp_rows)
  fwrite(grp,                     file.path(cfg$out_dir, "cell_line_groups.csv"))
  fwrite(rbindlist(thr_rows),     file.path(cfg$out_dir, "thresholds.csv"))
  fwrite(rbindlist(summary_rows), file.path(cfg$out_dir, "comparison_summary.csv"))

  if (length(auc_plot_data)) {
    apd <- rbindlist(auc_plot_data)
    plot_auc(apd, rbindlist(auc_plot_thr), file.path(cfg$out_dir, "AUC_by_cluster"),
             sprintf("%s: AUC and cutoffs", cfg$drug_name),
             ylab = if (cfg$auc_source == "official") sprintf("AUC PRISM (%s)", cfg$screen_id) else "AUC (0.6 nM - 2.5 uM)",
             cl_levels = cluster_order(apd$cluster, cfg$cluster_colors))
  }
  if (!is.null(A)) {
    msg("\nDose-response plots ...")
    plot_dose_response_by_group(A$points, grp, file.path(cfg$out_dir, "dose_response_by_group"),
                                sprintf("%s (%s): dose-response by group", cfg$drug_name, cfg$screen_id))
    plot_dose_response_by_cluster(A$points, grp, cfg, cfg$out_dir)
  }

  # Biological control: in the ALL contrast expect high CDKN2A/CCNE1 and low RB1 in resistant lines
  if (length(ctrl_rows)) {
    ctrl <- rbindlist(ctrl_rows)
    fwrite(ctrl, file.path(cfg$out_dir, "biomarker_controls.csv"))
    msg("\nControls (cluster ALL; logFC > 0 = higher in resistant):")
    print(ctrl[cluster == "ALL", .(context, gene, logFC = round(logFC, 2), FDR = signif(adj.P.Val, 2))],
          row.names = FALSE)
  }
  msg("\nDone. Results in: %s", normalizePath(cfg$out_dir))
  invisible(list(groups = grp, summary = rbindlist(summary_rows)))
}


# Analysis 4. Palbociclib sensitivity vs genetic dependency on CDK4/CDK6 --------
#  Same cell lines in every model: PRISM AUC + CRISPR (Chronos) CDK4 and CDK6 +
#  expression (+ a cluster assignment, as in the ALL comparison of Analysis 3).
#  Orientation (both phenotypes point the same way):
#     AUC high          = resistant to palbociclib
#     Chronos high (~0) = NOT dependent on the gene
#  -> positive coefficient = higher expression in resistant / non-dependent lines
#  Predictors are z-scored, so every coefficient is "change in log2(TPM+1) per
#  1 SD of the phenotype" and drug and CRISPR effects are comparable.
#  Models (limma-trend, robust eBayes, adjusted for lineage):
#     drug_AUC         expression ~ AUC
#     CRISPR_CDK4      expression ~ CDK4 gene effect
#     CRISPR_CDK6      expression ~ CDK6 gene effect
#     CRISPR_CDK4_6    expression ~ combined CDK4/6 dependency (cfg$combined)
#     joint            expression ~ AUC + CDK4 + CDK6
#                      -> drug_given_CRISPR: drug association NOT explained by
#                         CDK4/CDK6 genetic dependency ("pharmacology-specific")

CFG_PRISM_DRUG_VS_CRISPR <- modifyList(PRISM_CFG, list(
  restrict_to_clustered = TRUE,       # only lines with a cluster (same set as the ALL comparison)

  # --- CRISPR ---
  genes         = c("CDK4", "CDK6"),
  combined      = "min",              # CDK4/6 dependency: "min" (strongest of both) | "mean"
  dep_threshold = -0.5,               # Chronos value drawn as "dependent" in the scatter plots

  # --- AUC cutoffs (only for the reference lines in the scatter plots) ---
  q_sens       = 0.25,
  q_res        = 0.75,
  auc_sens_max = 0.80,

  # --- analysis / plots ---
  min_n_group   = 3,
  min_lineage_n = 3,
  fdr           = 0.05,
  label_n       = 15,
  volcano_label_genes = c("RB1", "CDKN2A", "CDKN2B", "CDKN2D", "CCNE1", "CCNE2",
                          "CCNG1", "CCNP", "ESD", "TACSTD2"),
  out_dir = file.path(DIRS$out, "04_palbociclib_drug_vs_CRISPR")
))

run_prism_drug_vs_crispr <- function(cfg) {
  cfg <- resolve_inputs(cfg, need = c("auc", "expr", "crispr"))
  odir <- file.path(cfg$out_dir, "pancancer")
  dir.create(odir, showWarnings = FALSE, recursive = TRUE)

  annot <- load_annotation(cfg)
  auc   <- load_official_auc(cfg)
  auc[, AUC := if (is.null(cfg$auc_cap)) auc_official else pmin(auc_official, cfg$auc_cap)]
  dep   <- load_crispr(cfg$crispr_file, cfg$genes)
  setnames(dep, cfg$genes, c("CDK4", "CDK6"))
  dep[, CDK4_6 := if (cfg$combined == "min") pmin(CDK4, CDK6) else (CDK4 + CDK6) / 2]
  msg("Loading expression ...")
  E   <- load_expr(cfg$expr_file, cfg$profile_map_file)
  cls <- load_clusters_any(cfg, annot)

  d <- merge(auc[, .(depmap_id, AUC)], dep, by = "depmap_id")
  d <- d[complete.cases(d) & depmap_id %in% colnames(E)]
  d <- merge(d, annot[, .(depmap_id = ModelID, lineage)], by = "depmap_id", all.x = TRUE)
  d <- merge(d, cls$pancancer[, .(depmap_id = ModelID, cluster)], by = "depmap_id", all.x = TRUE)
  if (cfg$restrict_to_clustered) d <- d[!is.na(cluster)]
  msg("\n%d cell lines with AUC + CRISPR CDK4/CDK6 + expression%s",
      nrow(d), if (cfg$restrict_to_clustered) " + cluster" else "")
  if (nrow(d) < 30) stop("Too few cell lines in common.")

  # Phenotype correlations
  cors <- rbindlist(lapply(c("CDK4", "CDK6", "CDK4_6"), function(g) {
    ct <- suppressWarnings(cor.test(d$AUC, d[[g]], method = "spearman"))
    data.table(comparison = paste("AUC vs", g), rho = unname(ct$estimate), p = ct$p.value)
  }))
  cors <- rbind(cors, {
    ct <- suppressWarnings(cor.test(d$CDK4, d$CDK6, method = "spearman"))
    data.table(comparison = "CDK4 vs CDK6", rho = unname(ct$estimate), p = ct$p.value) })
  msg("Spearman correlations between phenotypes:")
  print(cors[, .(comparison, rho = round(rho, 2), p = signif(p, 2))], row.names = FALSE)
  fwrite(cors, file.path(cfg$out_dir, "phenotype_correlations.csv"))

  # ---- Scatter: drug sensitivity vs genetic dependency (coloured by cluster) ----
  gr <- assign_groups(d$AUC, cfg)
  long <- melt(d, id.vars = c("depmap_id", "AUC", "cluster"), measure.vars = c("CDK4", "CDK6", "CDK4_6"),
               variable.name = "gene", value.name = "gene_effect")
  long[, gene := factor(gene, levels = c("CDK4", "CDK6", "CDK4_6"),
                        labels = c("CRISPR CDK4", "CRISPR CDK6",
                                   sprintf("CRISPR CDK4/6 (%s)", cfg$combined)))]
  quad <- long[, .(lab = sprintf("rho = %.2f\ndep. & sensitive: %d\ndep. & resistant: %d\nnon-dep. & sensitive: %d",
                                 cor(AUC, gene_effect, method = "spearman"),
                                 sum(gene_effect < cfg$dep_threshold & AUC <= gr$thr_sens),
                                 sum(gene_effect < cfg$dep_threshold & AUC >= gr$thr_res),
                                 sum(gene_effect >= cfg$dep_threshold & AUC <= gr$thr_sens))), by = gene]
  lv   <- cluster_order(long$cluster, cfg$cluster_colors)
  long[, cluster := factor(cluster, levels = lv)]
  g <- ggplot(long, aes(gene_effect, AUC, colour = cluster)) +
    geom_vline(xintercept = cfg$dep_threshold, linetype = 2, colour = "grey50", linewidth = 0.3) +
    geom_hline(yintercept = gr$thr_sens, linetype = 2, colour = GROUP_COLS[["sensitive"]], linewidth = 0.3) +
    geom_hline(yintercept = gr$thr_res,  linetype = 2, colour = GROUP_COLS[["resistant"]], linewidth = 0.3) +
    geom_point(size = 1.2, alpha = 0.8) +
    geom_text(data = quad, aes(x = -Inf, y = -Inf, label = lab), inherit.aes = FALSE,
              hjust = -0.05, vjust = -0.15, size = 2.6, lineheight = 0.9) +
    facet_wrap(~ gene, scales = "free_x") +
    scale_y_continuous(expand = expansion(mult = c(0.3, 0.04))) +   # room for the counts below the points
    scale_colour_manual(values = cluster_palette(lv, cfg$cluster_colors), name = "Cluster") +
    labs(title = sprintf("%s sensitivity vs genetic dependency (pan-cancer, n=%d)", cfg$drug_name, nrow(d)),
         subtitle = sprintf("Dashed: Chronos = %.1f | AUC sensitive <= %.2f, resistant >= %.2f",
                            cfg$dep_threshold, gr$thr_sens, gr$thr_res),
         x = "CRISPR gene effect (Chronos; lower = more dependent)",
         y = sprintf("AUC PRISM (%s)", cfg$screen_id)) +
    theme_bw(10) + theme(panel.grid = element_blank())
  save_plot(g, file.path(cfg$out_dir, "scatter_AUC_vs_CRISPR"), width = 12, height = 4.8)

  # ---- Differential expression models (same lines, z-scored predictors) ----
  Z  <- d[, .(AUC = as.numeric(scale(AUC)), CDK4 = as.numeric(scale(CDK4)),
              CDK6 = as.numeric(scale(CDK6)), CDK4_6 = as.numeric(scale(CDK4_6)))]
  Es <- E[, d$depmap_id, drop = FALSE]
  comb_name <- "CRISPR_CDK4_6"

  runs <- list(
    drug_AUC    = fit_multi(Es, Z[, .(AUC)],    d$lineage, "AUC",    cfg),
    CRISPR_CDK4 = fit_multi(Es, Z[, .(CDK4)],   d$lineage, "CDK4",   cfg),
    CRISPR_CDK6 = fit_multi(Es, Z[, .(CDK6)],   d$lineage, "CDK6",   cfg),
    CRISPR_CDK4_6 = fit_multi(Es, Z[, .(CDK4_6)], d$lineage, "CDK4_6", cfg),
    joint       = fit_multi(Es, Z[, .(AUC, CDK4, CDK6)], d$lineage, c("AUC", "CDK4", "CDK6"), cfg)
  )
  res <- list(drug_AUC          = runs$drug_AUC$res$AUC,
              CRISPR_CDK4       = runs$CRISPR_CDK4$res$CDK4,
              CRISPR_CDK6       = runs$CRISPR_CDK6$res$CDK6,
              CRISPR_CDK4_6     = runs$CRISPR_CDK4_6$res$CDK4_6,
              drug_given_CRISPR = runs$joint$res$AUC,
              CDK4_given_drug_CDK6 = runs$joint$res$CDK4,
              CDK6_given_drug_CDK4 = runs$joint$res$CDK6)
  for (nm in names(res)) fwrite(res[[nm]], file.path(odir, sprintf("DE_%s.csv", nm)))

  smry <- rbindlist(lapply(names(res), function(nm) {
    r <- res[[nm]]
    data.table(model = nm, n_lines = nrow(d), n_genes = nrow(r),
               n_up_resistant_or_nondependent = sum(r$adj.P.Val < cfg$fdr & r$logFC > 0),
               n_up_sensitive_or_dependent    = sum(r$adj.P.Val < cfg$fdr & r$logFC < 0))
  }))
  msg("\nDE genes (padj < %.2g) per model:", cfg$fdr)
  print(smry[, .(model, n_up_resistant_or_nondependent, n_up_sensitive_or_dependent)], row.names = FALSE)
  fwrite(smry, file.path(cfg$out_dir, "model_summary.csv"))
  if (!runs$drug_AUC$lineage_adjusted) msg("NOTE: lineage could not be added as covariate.")

  # ---- Drug vs CRISPR: which genes are shared and which differ ----
  comps <- list("CRISPR CDK4" = "CRISPR_CDK4", "CRISPR CDK6" = "CRISPR_CDK6",
                "CRISPR CDK4/6" = "CRISPR_CDK4_6")
  cc <- rbindlist(lapply(names(comps), function(lab) {
    m <- compare_models(res$drug_AUC, res[[comps[[lab]]]], cfg$fdr)
    cbind(comparison = paste("drug AUC vs", lab), m)
  }))
  fwrite(cc[order(comparison, category == "NS", padj_drug)], file.path(cfg$out_dir, "drug_vs_CRISPR_genes.csv"))
  tab <- dcast(cc[, .N, by = .(comparison, category)], comparison ~ category, value.var = "N", fill = 0)
  msg("\nGenes by category (drug vs CRISPR):"); print(tab, row.names = FALSE)
  fwrite(tab, file.path(cfg$out_dir, "drug_vs_CRISPR_counts.csv"))

  rho <- cc[, .(lab = sprintf("rho = %.2f", cor(t_drug, t_crispr, method = "spearman"))), by = comparison]
  comp_levels <- paste("drug AUC vs", names(comps))
  rho[, comparison := factor(comparison, levels = comp_levels)]
  cc_plot <- cc[order(category != "NS")]
  cc_plot[, comparison := factor(comparison, levels = comp_levels)]
  cc_plot[, category := factor(category, levels = names(COMP_COLS))]
  lab_pts <- cc_plot[category != "NS"][order(pmin(padj_drug, padj_crispr)), head(.SD, 8), by = .(comparison, category)]
  g <- ggplot(cc_plot, aes(t_drug, t_crispr, colour = category)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_point(size = 0.6, alpha = 0.7) +
    geom_text(data = rho, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
              hjust = -0.1, vjust = 1.5, size = 3) +
    facet_wrap(~ comparison) +
    scale_colour_manual(values = COMP_COLS, drop = FALSE) +
    labs(title = sprintf("Expression correlates of %s sensitivity vs CDK4/CDK6 dependency", cfg$drug_name),
         subtitle = "t statistics; positive = higher in resistant (drug) / non-dependent (CRISPR) lines",
         x = "t, drug AUC model", y = "t, CRISPR gene effect model", colour = NULL) +
    theme_bw(10) + theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size = 2.5)))
  if (requireNamespace("ggrepel", quietly = TRUE))
    g <- g + ggrepel::geom_text_repel(data = lab_pts, aes(label = gene), colour = "black", size = 2.3,
                                      max.overlaps = Inf, min.segment.length = 0, segment.size = 0.2,
                                      show.legend = FALSE)
  save_plot(g, file.path(cfg$out_dir, "concordance_drug_vs_CRISPR"), width = 12, height = 5)

  # ---- Volcanos: all models side by side + drug-specific (joint model) ----
  vol <- rbindlist(lapply(c("drug_AUC", "CRISPR_CDK4", "CRISPR_CDK6", "CRISPR_CDK4_6", "drug_given_CRISPR"),
                          function(nm) res[[nm]][, .(gene, logFC, P.Value, adj.P.Val, panel = nm)]))
  vol[, panel := factor(panel, levels = c("drug_AUC", "CRISPR_CDK4", "CRISPR_CDK6", "CRISPR_CDK4_6", "drug_given_CRISPR"),
                        labels = c("drug AUC", "CRISPR CDK4", "CRISPR CDK6",
                                   sprintf("CRISPR CDK4/6 (%s)", cfg$combined), "drug | CRISPR (joint model)"))]
  vol[, panel := as.character(panel)]
  plot_volcano(vol, file.path(cfg$out_dir, "volcano_models"), cfg$fdr, 0,
               sprintf("%s: expression vs drug sensitivity and CDK4/6 dependency (pan-cancer)", cfg$drug_name),
               lab_up = "Higher in resistant / non-dependent", lab_down = "Higher in sensitive / dependent",
               xlab = "Change in log2(TPM+1) per 1 SD of the phenotype",
               subtitle = sprintf("Same %d lines, lineage-adjusted; coloured: padj < %.2g", nrow(d), cfg$fdr))
  for (lab in c(FALSE, TRUE))
    plot_volcano(res$drug_given_CRISPR[, .(gene, logFC, P.Value, adj.P.Val, panel = "drug | CRISPR CDK4 + CDK6")],
                 file.path(cfg$out_dir, paste0("volcano_drug_given_CRISPR", if (lab) "_labeled" else "")),
                 cfg$fdr, 0,
                 sprintf("%s: drug association beyond CDK4/6 genetic dependency", cfg$drug_name),
                 lab_up = "Higher with resistance", lab_down = "Higher with sensitivity",
                 xlab = "Change in log2(TPM+1) per 1 SD AUC (adjusted for CDK4 and CDK6 gene effect)",
                 subtitle = "expression ~ AUC + CDK4 + CDK6 + lineage",
                 label = lab, label_n = cfg$label_n, label_genes = cfg$volcano_label_genes)

  msg("\nDone. Results in: %s", normalizePath(cfg$out_dir))
  invisible(res)
}


# Analysis 5. Palbociclib AUC vs CRISPR gene effect (CDK4, CDK6, CCND1) ----------
#  Figures (each saved with and without cell line labels, PNG + PDF):
#    1. pan-cancer           : one panel per gene
#    2. breast               : one panel per gene
#    3. pan-cancer by cluster: one figure per gene, one panel per cluster
#    4. summary              : correlation (95% CI) per cluster and gene
#  x = CRISPR gene effect (Chronos; lower = more dependent)
#  y = AUC PRISM MTS010 (lower = more sensitive; capped at 1)
#  Positive rho = more dependent lines are more sensitive to palbociclib.

CFG_PRISM_AUC_VS_CRISPR <- modifyList(PRISM_CFG, list(
  genes                 = c("CDK4", "CDK6", "CCND1"),
  restrict_to_clustered = TRUE,       # pan-cancer/breast: only lines with a cluster (as in the paper)

  # --- statistics / plots ---
  cor_method      = "spearman",       # "spearman" | "pearson"
  min_n_cluster   = 5,                # minimum lines to compute a correlation in a cluster
  dep_threshold   = -0.5,             # dashed vertical line (Chronos)
  label_n         = 12,               # labelled lines per panel (extremes + outliers)
  label_all_max   = 40,               # panels with <= this many lines: label every line
  label_n_cluster = 6,                # labelled lines per panel in the by-cluster figures
  label_all_max_cluster = 10,         # by-cluster panels with <= this many lines: label every line
  label_lines     = NULL,             # cell lines always labelled, e.g. c("MCF7", "T47D")
  gene_colors     = c(CDK4 = "#D95F02", CDK6 = "#7570B3", CCND1 = "#1B9E77"),
  out_dir = file.path(DIRS$out, "05_palbociclib_AUC_vs_CRISPR")
))

run_prism_auc_vs_crispr <- function(cfg) {
  cfg <- resolve_inputs(cfg, need = c("auc", "crispr"))
  dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)

  annot <- load_annotation(cfg)
  auc   <- load_official_auc(cfg)
  auc[, auc := if (is.null(cfg$auc_cap)) auc_official else pmin(auc_official, cfg$auc_cap)]
  dep   <- load_crispr(cfg$crispr_file, cfg$genes)
  cls   <- load_clusters_any(cfg, annot)

  d <- merge(auc[, .(depmap_id, auc)], dep, by = "depmap_id")
  d <- merge(d, annot[, .(depmap_id = ModelID, name = fcoalesce(CellLineName, StrippedCellLineName, ModelID),
                          lineage, is_breast)], by = "depmap_id", all.x = TRUE)
  d[is.na(name), name := depmap_id]
  d_pan <- merge(d, cls$pancancer[, .(depmap_id = ModelID, cluster)], by = "depmap_id", all.x = TRUE)
  d_br  <- merge(d[is_breast %in% TRUE], cls$breast[, .(depmap_id = ModelID, cluster)], by = "depmap_id", all.x = TRUE)
  if (cfg$restrict_to_clustered) { d_pan <- d_pan[!is.na(cluster)]; d_br <- d_br[!is.na(cluster)] }
  msg("\nLines with AUC + CRISPR%s: pan-cancer %d | breast %d",
      if (cfg$restrict_to_clustered) " + cluster" else "", nrow(d_pan), nrow(d_br))

  to_long <- function(x) {
    l <- melt(x, id.vars = c("depmap_id", "name", "cluster", "auc"), measure.vars = cfg$genes,
              variable.name = "gene", value.name = "gene_effect")
    l[, gene := as.character(gene)][!is.na(gene_effect)]
  }
  L_pan <- to_long(d_pan); L_br <- to_long(d_br)
  gene_panel <- function(g) factor(paste("CRISPR", g), levels = paste("CRISPR", cfg$genes))
  sub_txt <- sprintf("%s correlation with linear fit; positive = more dependent lines are more sensitive",
                     tools::toTitleCase(cfg$cor_method))

  # 1. Pan-cancer
  L_pan[, panel := gene_panel(gene)]
  st_pan <- plot_corr(L_pan, file.path(cfg$out_dir, "corr_AUC_vs_CRISPR_pancancer"),
                      sprintf("%s sensitivity vs CRISPR dependency - pan-cancer", cfg$drug_name), sub_txt, cfg)
  # 2. Breast
  st_br <- NULL
  if (nrow(L_br)) {
    L_br[, panel := gene_panel(gene)]
    st_br <- plot_corr(L_br, file.path(cfg$out_dir, "corr_AUC_vs_CRISPR_breast"),
                       sprintf("%s sensitivity vs CRISPR dependency - breast", cfg$drug_name), sub_txt, cfg)
  } else msg("No breast lines: breast figure skipped.")

  # 3. Pan-cancer by cluster (one figure per gene)
  ncl <- L_pan[, .(n = uniqueN(depmap_id)), by = cluster]
  keep_cl <- ncl[n >= cfg$min_n_cluster, cluster]
  if (length(setdiff(ncl$cluster, keep_cl)))
    msg("Clusters with < %d lines (no correlation): %s", cfg$min_n_cluster,
        paste(setdiff(ncl$cluster, keep_cl), collapse = ", "))
  lv <- cluster_order(keep_cl, cfg$cluster_colors)
  st_cl <- list()
  for (g0 in cfg$genes) {
    x <- L_pan[gene == g0 & cluster %in% keep_cl]
    x[, panel := factor(cluster, levels = lv)]
    st <- plot_corr(x, file.path(cfg$out_dir, sprintf("corr_AUC_vs_CRISPR_by_cluster_%s", g0)),
                    sprintf("%s sensitivity vs CRISPR %s, by cluster (pan-cancer)", cfg$drug_name, g0),
                    sub_txt, cfg, fill_by = "none", ncol = 4, label_n = cfg$label_n_cluster,
                    label_all_max = cfg$label_all_max_cluster,
                    width = 11, height = 2.7 * ceiling(length(lv) / 4) + 1.2)
    st_cl[[g0]] <- cbind(gene = g0, st[, .(cluster = as.character(panel), n, rho, lo, hi, p)])
  }

  # 4. Summary table and forest plot
  tab <- rbind(
    cbind(context = "Pan-cancer", st_pan[, .(gene = sub("^CRISPR ", "", panel), cluster = "All", n, rho, lo, hi, p)]),
    if (!is.null(st_br)) cbind(context = "Breast", st_br[, .(gene = sub("^CRISPR ", "", panel), cluster = "All", n, rho, lo, hi, p)]),
    cbind(context = "Pan-cancer by cluster", rbindlist(st_cl)[, .(gene, cluster, n, rho, lo, hi, p)])
  )
  tab[context == "Pan-cancer by cluster", padj := p.adjust(p, "BH"), by = gene]
  fwrite(tab, file.path(cfg$out_dir, "correlations_AUC_vs_CRISPR.csv"))
  msg("\n%s correlation AUC vs gene effect:", tools::toTitleCase(cfg$cor_method))
  print(dcast(tab, context + cluster ~ gene, value.var = "rho")[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 2) else v)],
        row.names = FALSE)

  fp <- copy(tab)
  fp[, row := fifelse(context == "Pan-cancer", "Pan-cancer", fifelse(context == "Breast", "Breast", cluster))]
  fp[, `:=`(row = factor(row, levels = rev(c("Pan-cancer", if (!is.null(st_br)) "Breast", lv))),
            gene = factor(gene, levels = cfg$genes), sig = p < 0.05)]
  n_ref <- 1 + !is.null(st_br)
  g <- ggplot(fp, aes(rho, row, colour = gene, group = factor(gene, levels = rev(cfg$genes)))) +   # CDK4 on top
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_hline(yintercept = length(lv) + 0.5, colour = "grey30", linewidth = 0.3) +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, orientation = "y",
                  position = position_dodge(width = 0.6), linewidth = 0.4) +
    geom_point(aes(shape = sig), size = 2.2, fill = "white", position = position_dodge(width = 0.6)) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), labels = c(`TRUE` = "p < 0.05", `FALSE` = "p >= 0.05"),
                       name = NULL) +
    scale_colour_manual(values = cfg$gene_colors[cfg$genes], name = "CRISPR") +
    scale_x_continuous(limits = c(-1, 1)) +
    labs(title = sprintf("%s AUC vs CRISPR gene effect", cfg$drug_name),
         subtitle = sprintf("%s rho (95%% CI); positive = dependent lines are more sensitive",
                            tools::toTitleCase(cfg$cor_method)),
         x = "rho", y = NULL) +
    theme_bw(10) + theme(panel.grid.minor = element_blank(), legend.position = "right")
  save_plot(g, file.path(cfg$out_dir, "corr_AUC_vs_CRISPR_summary"), width = 7.5, height = 0.42 * (length(lv) + n_ref) + 2)

  msg("\nDone. Results in: %s", normalizePath(cfg$out_dir))
  invisible(tab)
}


# Analysis 6. Expression correlates: palbociclib sensitivity vs dependency ------
#  For every expressed gene:
#     rho_drug   = Spearman(expression, palbociclib AUC)
#     rho_CRISPR = Spearman(expression, CRISPR gene effect of CDK4 / CDK6 / CCND1)
#  computed on the SAME cell lines (AUC + CRISPR + expression [+ cluster]).
#  Each point of the plots is a gene: x = rho_drug, y = rho_CRISPR.
#  Orientation (cfg$flip = FALSE):
#     rho_drug   > 0 : gene higher in RESISTANT lines (high AUC)
#     rho_CRISPR > 0 : gene higher in NON-DEPENDENT lines (high Chronos)
#  cfg$flip = TRUE multiplies both by -1 (positive = sensitive / dependent).
#  Contexts: pan-cancer, breast and pan-cancer by cluster. Every figure is
#  saved with and without gene labels (PNG + PDF).

CFG_PRISM_COR_EXPRESSION <- modifyList(PRISM_CFG, list(
  genes                 = c("CDK4", "CDK6", "CCND1"),
  restrict_to_clustered = TRUE,       # only lines with a cluster (same set as in the paper)
  min_n_cluster         = 10,         # minimum lines to correlate within a cluster

  # --- statistics ---
  flip           = FALSE,             # TRUE: positive = higher in sensitive / dependent lines
  adjust_lineage = FALSE,             # TRUE: partial Spearman adjusting for lineage (pan-cancer, clusters)
  min_lineage_n  = 3,
  fdr            = 0.05,
  rho_min        = 0.1,               # |rho| needed (with FDR) to call a gene significant

  # --- labels ---
  label_per_category         = 8,     # labelled genes per category and panel
  label_per_category_cluster = 3,
  label_genes = c("RB1", "CDKN2A", "CDKN2B", "CDKN2D", "CCNE1", "CCNE2", "CCNG1", "CCNP",
                  "ESD", "TACSTD2", "CDK4", "CDK6", "CCND1"),   # always labelled (ringed)
  out_dir = file.path(DIRS$out, "06_palbociclib_expression_drug_vs_CRISPR")
))

run_prism_cor_expression <- function(cfg) {
  cfg <- resolve_inputs(cfg, need = c("auc", "expr", "crispr"))
  dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)

  annot <- load_annotation(cfg)
  auc   <- load_official_auc(cfg)
  auc[, AUC := if (is.null(cfg$auc_cap)) auc_official else pmin(auc_official, cfg$auc_cap)]
  dep   <- load_crispr(cfg$crispr_file, cfg$genes)
  msg("Loading expression ...")
  E   <- load_expr(cfg$expr_file, cfg$profile_map_file)
  cls <- load_clusters_any(cfg, annot)

  d <- merge(auc[, .(depmap_id, AUC)], dep, by = "depmap_id")
  d <- d[complete.cases(d) & depmap_id %in% colnames(E)]
  d <- merge(d, annot[, .(depmap_id = ModelID, lineage, is_breast)], by = "depmap_id", all.x = TRUE)
  d_pan <- merge(d, cls$pancancer[, .(depmap_id = ModelID, cluster)], by = "depmap_id", all.x = TRUE)
  d_br  <- merge(d[is_breast %in% TRUE], cls$breast[, .(depmap_id = ModelID, cluster)], by = "depmap_id", all.x = TRUE)
  if (cfg$restrict_to_clustered) { d_pan <- d_pan[!is.na(cluster)]; d_br <- d_br[!is.na(cluster)] }
  msg("\nLines with AUC + CRISPR (%s) + expression%s: pan-cancer %d | breast %d",
      paste(cfg$genes, collapse = ", "), if (cfg$restrict_to_clustered) " + cluster" else "", nrow(d_pan), nrow(d_br))

  phen <- c("AUC", cfg$genes)
  run_context <- function(dd, adjust) {
    Y <- as.matrix(dd[, ..phen]); rownames(Y) <- dd$depmap_id
    spearman_genes(E[, dd$depmap_id, drop = FALSE], Y, if (adjust) dd$lineage else NULL, cfg)
  }
  panel_lab <- function(g) factor(paste("CRISPR", g), levels = paste("CRISPR", cfg$genes))
  sub_txt <- sprintf("Each point is a gene; Spearman%s; coloured: FDR < %.2g and |rho| >= %.2g; dotted: y = x",
                     if (cfg$adjust_lineage) " (lineage-adjusted)" else "", cfg$fdr, cfg$rho_min)
  counts <- list(); sims <- list()

  # 1-2. Pan-cancer and breast
  for (ctx in c("pancancer", "breast")) {
    dd <- if (ctx == "pancancer") d_pan else d_br
    if (nrow(dd) < cfg$min_n_cluster) { msg("[%s] only %d lines: skipped", ctx, nrow(dd)); next }
    cr <- run_context(dd, cfg$adjust_lineage && ctx == "pancancer")
    pp <- pair_up(cr, cfg)
    fwrite(pp, file.path(cfg$out_dir, sprintf("cor_expression_%s.csv", ctx)))
    pp[, panel := panel_lab(crispr_gene)]
    sim <- plot_rr(pp, file.path(cfg$out_dir, sprintf("cor_drug_vs_CRISPR_%s", ctx)),
                   sprintf("Expression correlates of %s sensitivity vs dependency - %s (n = %d lines)",
                           cfg$drug_name, if (ctx == "pancancer") "pan-cancer" else "breast", nrow(dd)), sub_txt, cfg)
    sims[[ctx]] <- cbind(context = ctx, cluster = "All", sim[, .(crispr_gene = sub("^CRISPR ", "", panel), r, n_shared)])
    counts[[ctx]] <- cbind(context = ctx, cluster = "All", pp[, .N, by = .(crispr_gene, category)])
    msg("[%s] genes by category:", ctx)
    print(dcast(pp[, .N, by = .(crispr_gene, category)], crispr_gene ~ category, value.var = "N", fill = 0), row.names = FALSE)
  }

  # 3. Pan-cancer by cluster
  ncl <- d_pan[!is.na(cluster), .N, by = cluster]
  keep_cl <- ncl[N >= cfg$min_n_cluster, cluster]
  if (length(setdiff(ncl$cluster, keep_cl)))
    msg("Clusters with < %d lines (skipped): %s", cfg$min_n_cluster, paste(setdiff(ncl$cluster, keep_cl), collapse = ", "))
  lv <- cluster_order(keep_cl, cfg$cluster_colors)
  pp_cl <- rbindlist(lapply(lv, function(k) {
    dd <- d_pan[cluster == k]
    cbind(cluster = k, pair_up(run_context(dd, cfg$adjust_lineage), cfg))
  }))
  if (nrow(pp_cl)) {
    fwrite(pp_cl, file.path(cfg$out_dir, "cor_expression_by_cluster.csv"))
    for (g0 in cfg$genes) {
      x <- pp_cl[crispr_gene == g0]
      nn <- d_pan[cluster %in% lv, .N, by = cluster]
      x[, panel := factor(sprintf("%s (n=%d)", cluster, nn$N[match(cluster, nn$cluster)]),
                          levels = sprintf("%s (n=%d)", lv, nn$N[match(lv, nn$cluster)]))]
      sim <- plot_rr(x, file.path(cfg$out_dir, sprintf("cor_drug_vs_CRISPR_%s_by_cluster", g0)),
                     sprintf("Expression correlates of %s sensitivity vs CRISPR %s, by cluster", cfg$drug_name, g0),
                     sub_txt, cfg, ncol = 4, label_per_cat = cfg$label_per_category_cluster,
                     width = 12, height = 3.3 * ceiling(length(lv) / 4) + 1.4)
      sims[[g0]] <- cbind(context = "pancancer_by_cluster", sim[, .(cluster = sub(" \\(n=.*$", "", panel),
                                                                    crispr_gene = g0, r, n_shared)])
    }
    counts[["clusters"]] <- cbind(context = "pancancer_by_cluster", pp_cl[, .N, by = .(cluster, crispr_gene, category)])
  }

  # 4. Summary: similarity between the drug and CRISPR expression programs
  smry <- rbindlist(sims, use.names = TRUE)
  fwrite(smry, file.path(cfg$out_dir, "similarity_drug_vs_CRISPR.csv"))
  fwrite(rbindlist(counts, use.names = TRUE, fill = TRUE), file.path(cfg$out_dir, "category_counts.csv"))
  smry[, column := fifelse(context == "pancancer", "Pan-cancer", fifelse(context == "breast", "Breast", cluster))]
  smry[, `:=`(column = factor(column, levels = c("Pan-cancer", "Breast", lv)),
              crispr_gene = factor(crispr_gene, levels = rev(cfg$genes)))]
  g <- ggplot(smry, aes(column, crispr_gene, fill = r)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.2f\n(%d)", r, n_shared)), size = 2.8, lineheight = 0.9) +
    geom_vline(xintercept = 2.5, colour = "black", linewidth = 0.5) +
    scale_fill_gradient2(low = "#2471A3", mid = "white", high = "#C0392B", midpoint = 0, limits = c(-1, 1),
                         name = "similarity r") +
    labs(title = sprintf("Similarity of the expression programs: %s sensitivity vs CRISPR dependency", cfg$drug_name),
         subtitle = "Spearman correlation between gene-wise rho_drug and rho_CRISPR (number of shared genes)",
         x = NULL, y = "CRISPR") +
    theme_minimal(10) + theme(panel.grid = element_blank())
  save_plot(g, file.path(cfg$out_dir, "similarity_drug_vs_CRISPR_summary"), width = 0.75 * nlevels(smry$column) + 3, height = 3.2)
  msg("\nSimilarity r (rho_drug vs rho_CRISPR):")
  print(dcast(smry, crispr_gene ~ column, value.var = "r")[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 2) else v)],
        row.names = FALSE)
  msg("\nDone. Results in: %s", normalizePath(cfg$out_dir))
  invisible(smry)
}


# Run -----------------------------------------------------------------------------
dir.create(DIRS$out, recursive = TRUE, showWarnings = FALSE)

if (RUN$adjusted_codependency || RUN$classifier_validation) depmap <- load_depmap(FILES)
if (RUN$adjusted_codependency) res_codependency <- run_adjusted_codependency(depmap)
if (RUN$classifier_validation) res_classifier   <- run_classifier_validation(depmap)
if (RUN$prism_de_by_cluster)   run_prism_de_by_cluster(CFG_PRISM_DE)
if (RUN$prism_drug_vs_crispr)  run_prism_drug_vs_crispr(CFG_PRISM_DRUG_VS_CRISPR)
if (RUN$prism_auc_vs_crispr)   run_prism_auc_vs_crispr(CFG_PRISM_AUC_VS_CRISPR)
if (RUN$prism_cor_expression)  run_prism_cor_expression(CFG_PRISM_COR_EXPRESSION)

writeLines(capture.output(sessionInfo()), file.path(DIRS$out, "sessionInfo.txt"))
