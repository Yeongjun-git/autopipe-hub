# Bulk RNA-seq DE Analysis Pipeline

Differential expression analysis using edgeR with GO/KEGG/Reactome/Hallmark enrichment, GSEA, PCA, and volcano plots.

## Required Inputs

Place in input directory:
- `salmon_metadata.csv` — CSV with columns: sample, file, Condition, Batch
- `tx2gene.tsv` — Transcript-to-gene mapping from preprocessing step
- `star_salmon/` — Salmon quant.sf files per sample (referenced in metadata)
- `hallmark/` (optional) — Hallmark GMT files
- `reference/` (optional) — Custom gene set GMT files

## Outputs

- `{project_name}/PCA_plot.png` — PCA visualization
- `{project_name}/0.DGE/` — DEG Excel files
- `{project_name}/1.volcano/` — Volcano plots
- `{project_name}/2.DEG_enrichment/` — GO/KEGG/Reactome/Hallmark enrichment
- `{project_name}/3.GSEA/` — GSEA results and figures

## How to Run

```bash
docker build -t autopipe-rnaseq-de .
docker run --rm -v /path/to/input:/input:ro -v /path/to/output:/output \
  autopipe-rnaseq-de snakemake --cores 4 --snakefile /pipeline/Snakefile
```