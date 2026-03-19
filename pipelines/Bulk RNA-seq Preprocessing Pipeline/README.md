# Bulk RNA-seq Preprocessing Pipeline

Runs the nf-core/rnaseq pipeline (STAR-Salmon workflow) for bulk RNA-seq preprocessing including QC, trimming, alignment, and transcript quantification.

## Required Inputs

Mount these under `/input`:

- **fastq_metadata.csv**: nf-core samplesheet CSV with columns `sample`, `fastq_1`, `fastq_2`, `strandedness`
- **reference_fasta**: Reference genome FASTA file (e.g., `Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa`)
- **reference_gtf**: Gene annotation GTF file (e.g., `Homo_sapiens.GRCh38.112.gtf`)
- **FASTQ files**: Raw sequencing reads referenced in the metadata CSV

## Expected Outputs

Written to `/output`:

- `results/star_salmon/` — Per-sample Salmon quantification (`quant.sf`), `tx2gene.tsv`
- `results/star_salmon/*/` — Individual sample directories with BAM and Salmon outputs
- `results/multiqc/` — MultiQC quality report
- `preprocessing_complete.flag` — Completion indicator

## How to Run

```bash
# Build the Docker image
docker build -t bulk_rnaseq_preprocess .

# Run the pipeline
docker run --rm \
  -v /path/to/input_data:/input:ro \
  -v /path/to/output:/output \
  -v /var/run/docker.sock:/var/run/docker.sock \
  bulk_rnaseq_preprocess \
  snakemake --cores 16 --snakefile /pipeline/Snakefile
```

**Note**: This pipeline requires Docker socket access (`-v /var/run/docker.sock:/var/run/docker.sock`) because nf-core/rnaseq runs sub-containers with `-profile docker`.

## Configuration

Edit `config.yaml` to adjust:

- `nextflow_version`: Nextflow version (default: 25.04.0)
- `nfcore_rnaseq_version`: nf-core/rnaseq release (default: 3.22.2)
- `aligner`: Alignment method (default: star_salmon)
- `threads`: Number of CPU threads (default: 16)
- `extra_nextflow_args`: Additional Nextflow CLI arguments