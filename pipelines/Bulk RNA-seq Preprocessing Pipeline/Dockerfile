FROM condaforge/miniforge3:latest

# Install system dependencies and Nextflow requirements
RUN conda install -c bioconda -c conda-forge -y \
    snakemake-minimal \
    bash \
    nextflow=24.10.4 \
    openjdk=17.0.10 \
    curl \
    wget \
    procps-ng \
    && conda clean -afy

# Fix bash for GLIBC compatibility
RUN ln -sf /opt/conda/bin/bash /usr/bin/bash && ln -sf /opt/conda/bin/bash /bin/sh

# Install Docker CLI (needed for -profile docker in nf-core)
RUN apt-get update && apt-get install -y docker.io && apt-get clean

# Copy pipeline files
COPY Snakefile /pipeline/Snakefile
COPY config.yaml /pipeline/config.yaml

WORKDIR /pipeline