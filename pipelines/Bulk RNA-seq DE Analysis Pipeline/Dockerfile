FROM condaforge/miniforge3:latest

# Install snakemake, bash, and R base via conda
RUN conda install -c conda-forge -c bioconda -y \
    snakemake-minimal \
    bash \
    r-base=4.3 \
    r-dplyr \
    r-tibble \
    r-readr \
    r-readxl \
    r-stringr \
    r-ggplot2 \
    r-ggrepel \
    r-jsonlite \
    && conda clean -afy

# Install Bioconductor packages via R
RUN R -e " \
    if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org'); \
    BiocManager::install(c( \
        'edgeR', \
        'DESeq2', \
        'tximport', \
        'RUVSeq', \
        'sva', \
        'clusterProfiler', \
        'topGO', \
        'fgsea', \
        'ReactomePA', \
        'org.Hs.eg.db', \
        'org.Mm.eg.db', \
        'biomaRt', \
        'enrichplot' \
    ), ask = FALSE, update = FALSE) \
"

# Install additional R packages from CRAN
RUN R -e " \
    install.packages(c( \
        'pheatmap', \
        'ggridges', \
        'RColorBrewer', \
        'writexl', \
        'openxlsx' \
    ), repos = 'https://cloud.r-project.org') \
"

# Fix bash for GLIBC compatibility
RUN ln -sf /opt/conda/bin/bash /usr/bin/bash && ln -sf /opt/conda/bin/bash /bin/sh

# Copy pipeline files
COPY Snakefile /pipeline/Snakefile
COPY config.yaml /pipeline/config.yaml
COPY scripts/ /pipeline/scripts/

WORKDIR /pipeline