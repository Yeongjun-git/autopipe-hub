# Bulk RNA-seq DE Analysis Pipeline

Performs differential expression analysis on Salmon quantification results using edgeR, followed by comprehensive downstream analyses including GO/KEGG/Reactome enrichment, GSEA, Hallmark analysis, and visualization (PCA, volcano plots).

## Required Inputs

Mount these under `/input`:

- **salmon_metadata.csv**: CSV with columns `sample`, `file` (path to quant.sf), group column (e.g., `Condition`), optional `Batch`
- **tx2gene.tsv**: Transcript-to-gene mapping TSV from STAR-Salmon preprocessing (3 columns: tx_id, gene_id, gene_name)
- **Salmon quant.sf directories**: Per-sample Salmon output directories containing `quant.sf` files (paths referenced in metadata CSV)
- **hallmark/** (optional): Directory with Hallmark GMT files (`h.all.v2025.1.Hs.symbols.gmt` for human, `mh.all.v2025.1.Mm.symbols.gmt` for mouse)
- **specific_geneset.gmt** (optional): Custom gene set GMT file for enrichment

## Expected Outputs

Written to `/output/<project_name>/`:

- `TPM_gene.csv`, `TPM_gene_log2p1.csv` — TPM expression tables
- `TMM_length_normalized_count.csv`, `log_TMM_length_normalized_count.csv` — TMM normalized counts
- `TMM_length_normalized_count_filtered.csv` — Filtered TMM counts
- `PCA_plot.png`, `PCA_plot_loading.csv` — PCA visualization and loading matrix
- `0.DGE/` — DEG result Excel files per comparison
- `1.volcano/`, `1.volcano_annot/`, `1.volcano_specific_genes/` — Volcano plots (basic, annotated, gene-specific)
- `2.DEG_enrichment/` — DEG-based enrichment results (GO BP/CC/MF, KEGG, Reactome, Hallmark, Custom)
- `3.GSEA/` — GSEA results and summary plots

## How to Run

```bash
# Build the Docker image
docker build -t bulk_rnaseq_de .

# Run the pipeline
docker run --rm \
  -v /path/to/input_data:/input:ro \
  -v /path/to/output:/output \
  bulk_rnaseq_de \
  snakemake --cores 4 --snakefile /pipeline/Snakefile
```

### Input directory structure example

```
/input/
├── salmon_metadata.csv
├── tx2gene.tsv
├── ConditionA_1/quant.sf
├── ConditionA_2/quant.sf
├── ConditionA_3/quant.sf
├── ConditionB_1/quant.sf
├── ConditionB_2/quant.sf
├── ConditionB_3/quant.sf
├── hallmark/
│   └── h.all.v2025.1.Hs.symbols.gmt
└── specific_geneset.gmt  (optional)
```

**Note**: The `file` column in `salmon_metadata.csv` should use `/input/` paths (e.g., `/input/ConditionA_1/quant.sf`).

## Configuration

Edit `config.yaml` to adjust:

- `project_name`: Output subfolder name (default: DE_results)
- `group_col`: Metadata column for group comparison (default: Condition)
- `comb_set`: JSON array of comparison pairs (default: `[["ConditionA", "ConditionB"]]`)
- `sample_organism`: "human" or "mouse"
- `keytype`: Gene ID format — "SYMBOL", "ENSEMBL", or "ENTREZID"
- `fc_threshold`: Log2 fold change cutoff (default: 1)
- `pval_threshold`: Adjusted p-value cutoff (default: 0.05)
- `norm_method`: "edgeR" or "DESeq2"
- `specific_genes`: Genes to highlight in volcano plots (JSON array)
- `specific_geneset_path`: Path to custom GMT file
- `hallmark_gmt_path`: Override default Hallmark GMT path