#!/usr/bin/env Rscript
# ============================================================
# Bulk RNA-seq DE Analysis Wrapper
# Parses CLI arguments and runs the edgeR + enrichment pipeline
# ============================================================

# --- Parse command-line arguments ---
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  params <- list()
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    if (i + 1 <= length(args) && !startsWith(args[i + 1], "--")) {
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

# --- Load R pipeline scripts ---
script_dir <- params$script_dir
source(file.path(script_dir, "00_libraries.R"))
source(file.path(script_dir, "01_deg_based_enrichment_analysis.R"))
source(file.path(script_dir, "02_gsea.R"))
source(file.path(script_dir, "03_visualization.R"))
source(file.path(script_dir, "04_pipeline.R"))

# --- Parse comb_set: "A,B;C,D" -> list(c("A","B"), c("C","D")) ---
parse_comb_set <- function(s) {
  pairs <- strsplit(s, ";")[[1]]
  lapply(pairs, function(p) {
    trimws(strsplit(p, ",")[[1]])
  })
}

# --- Parse specific_genes: "GENE1,GENE2" -> c("GENE1","GENE2") ---
parse_genes <- function(s) {
  if (is.null(s) || s == "" || s == "None" || s == "[]") return(NULL)
  trimws(strsplit(s, ",")[[1]])
}

# --- Build project_settings ---
root_path <- params$root_path

project_settings <- list(
  Project_name  = params$project_name,
  group_col     = params$group_col,
  comb_set      = parse_comb_set(params$comb_set),
  tx2gene_path  = params$tx2gene_path,
  metadata_path = params$metadata_path
)

message("========================================")
message("Bulk RNA-seq DE Analysis Pipeline")
message("========================================")
message("Project: ", project_settings$Project_name)
message("Group column: ", project_settings$group_col)
message("Comparisons: ", params$comb_set)
message("tx2gene: ", project_settings$tx2gene_path)
message("Metadata: ", project_settings$metadata_path)
message("Organism: ", params$organism)
message("Root path: ", root_path)
message("========================================")

# --- Run DEG pipeline ---
pipeline_result <- run_deg_pipeline(
  project_settings = project_settings,
  root_path = root_path
)

result_path <- pipeline_result$result_path
DGE_result  <- pipeline_result$DGE_result

# --- Parse optional parameters ---
specific_geneset <- if (!is.null(params$specific_geneset) && params$specific_geneset != "" && file.exists(params$specific_geneset)) {
  params$specific_geneset
} else {
  NULL
}

specific_genes <- parse_genes(params$specific_genes)

# --- Resolve Hallmark path ---
hallmark_path <- if (!is.null(params$hallmark_path) && params$hallmark_path != "" && file.exists(params$hallmark_path)) {
  params$hallmark_path
} else {
  NULL
}

message("Hallmark GMT: ", ifelse(is.null(hallmark_path), "Not provided or not found", hallmark_path))
message("Specific geneset: ", ifelse(is.null(specific_geneset), "Not provided or not found", specific_geneset))
message("Specific genes: ", ifelse(is.null(specific_genes), "None", paste(specific_genes, collapse = ", ")))

# --- Run downstream analysis for each comparison ---
for (comparison_name in names(DGE_result)) {
  message(">> Running analysis for: ", comparison_name)
  analysis_data(
    path            = result_path,
    Data            = DGE_result[[comparison_name]],
    method          = "padj",
    FC              = as.numeric(params$fc),
    pval            = as.numeric(params$pval),
    Comparison_name = comparison_name,
    Sample          = params$organism,
    specific_geneset = specific_geneset,
    specific_genes  = specific_genes,
    Keytype         = params$keytype,
    norm_method     = "edgeR",
    Hallmark_path   = hallmark_path
  )
}

message("========================================")
message("DE Analysis Complete!")
message("Results saved to: ", result_path)
message("========================================")