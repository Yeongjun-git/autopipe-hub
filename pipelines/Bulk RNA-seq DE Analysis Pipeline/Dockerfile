FROM continuumio/miniconda3

# ============================================================
# Create conda env with R 4.5.1
# Referenced from: 0_env/2_DE_Analysis/Dockerfile
# ============================================================
RUN conda create -n rnaseq_env r-base=4.5.1 -c conda-forge -y

SHELL ["conda", "run", "-n", "rnaseq_env", "/bin/bash", "-c"]

# ============================================================
# Install Bioconductor + R packages via conda
# ============================================================
RUN conda install -y -c bioconda -c conda-forge \
    bioconductor-tximport \
    bioconductor-edger \
    bioconductor-clusterprofiler \
    bioconductor-fgsea \
    bioconductor-reactomepa \
    bioconductor-org.mm.eg.db \
    bioconductor-org.hs.eg.db \
    bioconductor-enrichplot \
    r-ggplot2 \
    r-ggrepel \
    r-dplyr \
    r-tibble \
    r-writexl \
    r-stringr \
    r-rcolorbrewer \
    && conda clean -afy

# ============================================================
# Install snakemake + bash
# Pin python to 3.12 to avoid 3.14 conflict with snakemake
# ============================================================
RUN conda install -y -c bioconda -c conda-forge \
    python=3.12 \
    snakemake-minimal \
    bash \
    && conda clean -afy

# Fix bash for GLIBC compatibility
RUN ln -sf /opt/conda/envs/rnaseq_env/bin/bash /usr/bin/bash && \
    ln -sf /opt/conda/envs/rnaseq_env/bin/bash /bin/sh

# ============================================================
# Verify critical R packages
# ============================================================
RUN conda run -n rnaseq_env R -e "\
must <- c('writexl','tximport','edgeR','clusterProfiler','fgsea', \
  'ReactomePA','org.Hs.eg.db','org.Mm.eg.db','enrichplot', \
  'ggplot2','ggrepel','dplyr','tibble','stringr','RColorBrewer'); \
for (p in must) { library(p, character.only=TRUE); message(p,' OK') }; \
message('=== ALL PACKAGES OK ===')"

# ============================================================
# Copy pipeline files
# ============================================================
WORKDIR /pipeline
COPY Snakefile .
COPY config.yaml .
COPY scripts/ scripts/

# Activate conda env in PATH
ENV PATH=/opt/conda/envs/rnaseq_env/bin:$PATH
ENV CONDA_DEFAULT_ENV=rnaseq_env