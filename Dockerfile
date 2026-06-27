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
ENV PATH /opt/miniconda3/envs/bioinfo/bin:$PATH

WORKDIR /workspace

# Configure nbstripout in the container's git so notebook outputs are
# stripped on commit even when committing from inside the container.
RUN git config --global filter.nbstripout.clean nbstripout && \
    git config --global filter.nbstripout.smudge cat && \
    git config --global filter.nbstripout.required true && \
    git config --global diff.ipynb.textconv "nbstripout -t"

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "bioinfo"]
CMD ["/bin/bash"]
