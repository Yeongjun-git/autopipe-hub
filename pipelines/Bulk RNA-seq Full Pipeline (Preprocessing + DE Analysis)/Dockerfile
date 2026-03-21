FROM condaforge/miniforge3:latest

# ============================================================
# Bulk RNA-seq Pipeline Environment
# Step 1: Nextflow + nf-core/rnaseq (STAR-Salmon)
# Step 2: R + edgeR + clusterProfiler + enrichment/GSEA
# ============================================================

# --- System dependencies ---
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    curl \
    && rm -rf /var/lib/apt/lists/*

# --- Install Snakemake, Nextflow, and Java (for Nextflow) ---
RUN conda install -y -c bioconda -c conda-forge \
    snakemake-minimal \
    bash \
    nextflow=24.10.4 \
    openjdk=17.0.10 \
    && conda clean -afy

# Fix bash GLIBC mismatch
RUN ln -sf /opt/conda/bin/bash /usr/bin/bash && ln -sf /opt/conda/bin/bash /bin/sh

# --- Install R and Bioconductor packages ---
RUN conda install -y -c conda-forge -c bioconda \
    r-base=4.3.3 \
    r-ggplot2=3.5.0 \
    r-dplyr=1.1.4 \
    r-tibble=3.2.1 \
    r-stringr=1.5.1 \
    r-ggrepel=0.9.5 \
    r-writexl=1.5.0 \
    r-rcolorbrewer=1.1_3 \
    bioconductor-tximport=1.30.0 \
    bioconductor-edger=4.0.16 \
    bioconductor-clusterprofiler=4.10.1 \
    bioconductor-fgsea=1.28.0 \
    bioconductor-reactomepa=1.46.0 \
    bioconductor-org.hs.eg.db=3.18.0 \
    bioconductor-org.mm.eg.db=3.18.0 \
    bioconductor-enrichplot=1.22.0 \
    && conda clean -afy

# --- Copy pipeline files ---
COPY Snakefile /pipeline/Snakefile
COPY config.yaml /pipeline/config.yaml
COPY run_de_analysis.R /pipeline/run_de_analysis.R
COPY nextflow.config /pipeline/nextflow.config
COPY R_scripts/ /pipeline/R_scripts/

WORKDIR /pipeline