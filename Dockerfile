FROM anaconda/miniconda:latest

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# System dependencies often needed by bioinformatics tools
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    curl \
    git \
    unzip \
    libz-dev \
    && rm -rf /var/lib/apt/lists/*

# Configure conda channels (order matters: priority is top-down, conda-forge top of priority list)
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict

RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Copy environment file and create the environment
COPY environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml && conda clean -afy

# Activate environment by default
SHELL ["conda", "run", "-n", "bioinfo", "/bin/bash", "-c"]
ENV PATH /opt/conda/envs/bioinfo/bin:$PATH

WORKDIR /workspace

CMD ["/bin/bash"]
