#!/usr/bin/env Rscript
# ==============================================================================
# Bulk RNA-seq DE Analysis Pipeline Runner
# Parses CLI arguments and runs the full DEG + enrichment + GSEA pipeline
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  params <- list()
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    if (i + 1 <= length(args) && !grepl("^--", args[i + 1])) {
      params[[key]] <- args[i + 1]
      i <- i + 2
    } else {
      params[[key]] <- TRUE
      i <- i + 1
    }
  }
  return(params)
}

params <- parse_args(args)

# --- Load pipeline R scripts ---
script_dir <- "/pipeline/scripts"
source(file.path(script_dir, "00_libraries.R"))
source(file.path(script_dir, "01_deg_based_enrichment_analysis.R"))
source(file.path(script_dir, "02_gsea.R"))
source(file.path(script_dir, "03_visualization.R"))
source(file.path(script_dir, "04_pipeline.R"))

# --- Parse parameters ---
metadata_path   <- params$metadata_path
tx2gene_path    <- params$tx2gene_path
project_name    <- params$project_name
group_col       <- params$group_col
output_dir      <- params$output_dir
sample_organism <- params$sample_organism
keytype         <- params$keytype
fc_threshold    <- as.numeric(params$fc_threshold)
pval_threshold  <- as.numeric(params$pval_threshold)
pval_method     <- params$pval_method
norm_method     <- params$norm_method

# Parse comb_set: "A,B;C,D" -> list(c("A","B"), c("C","D"))
parse_comb_set <- function(s) {
  pairs <- strsplit(s, ";")[[1]]
  lapply(pairs, function(p) trimws(strsplit(p, ",")[[1]]))
}
comb_set <- parse_comb_set(params$comb_set)

# Parse specific_genes
parse_genes <- function(s) {
  if (is.null(s) || s == "" || s == "None" || s == "[]") return(NULL)
  g <- trimws(strsplit(gsub("[\\[\\]'\"]", "", s), ",")[[1]])
  g[g != ""]
}
specific_genes <- parse_genes(params$specific_genes)

# Specific geneset GMT path
specific_geneset <- params$specific_geneset
if (is.null(specific_geneset) || specific_geneset == "" || !file.exists(specific_geneset)) {
  specific_geneset <- NULL
}

# Hallmark GMT path
hallmark_gmt <- params$hallmark_gmt
if (!is.null(hallmark_gmt) && hallmark_gmt != "" && file.exists(hallmark_gmt)) {
  assign("HALLMARK_OVERRIDE", hallmark_gmt, envir = .GlobalEnv)
}

# --- Set root path ---
root_path <- paste0(output_dir, "/")

# --- Build project_settings ---
project_settings <- list(
  Project_name  = paste0(project_name, "/"),
  group_col     = group_col,
  comb_set      = comb_set,
  tx2gene_path  = tx2gene_path,
  metadata_path = metadata_path
)

message("========================================")
message("Bulk RNA-seq DE Analysis Pipeline")
message("========================================")
message("Project: ", project_name)
message("Group column: ", group_col)
message("Organism: ", sample_organism)
message("Keytype: ", keytype)
message("FC threshold: ", fc_threshold)
message("P-value threshold: ", pval_threshold)
message("Method: ", pval_method)
message("Norm method: ", norm_method)
message("Comparisons: ", length(comb_set))
for (i in seq_along(comb_set)) {
  message("  ", paste(comb_set[[i]], collapse = " vs "))
}
message("========================================")

# --- Run DEG pipeline ---
pipeline_result <- run_deg_pipeline(
  project_settings = project_settings,
  root_path = root_path
)

result_path <- pipeline_result$result_path
DGE_result  <- pipeline_result$DGE_result

message("DEG pipeline complete. Running downstream analysis...")

# --- Run downstream analysis for each comparison ---
for (comparison_name in names(DGE_result)) {
  message("Processing comparison: ", comparison_name)
  analysis_data(
    path            = result_path,
    Data            = DGE_result[[comparison_name]],
    method          = pval_method,
    FC              = fc_threshold,
    pval            = pval_threshold,
    Comparison_name = comparison_name,
    Sample          = sample_organism,
    specific_geneset = specific_geneset,
    specific_genes  = specific_genes,
    Keytype         = keytype,
    norm_method     = norm_method
  )
}

message("========================================")
message("All analyses completed successfully!")
message("Results saved to: ", result_path)
message("========================================")