# ---------------------------------------------------------------------------
# Dockerfile for the MKWS NH3-H2 flame project
#
# Builds a reproducible Linux environment with Cantera and the Python stack
# needed to run the simulations and generate the figures.
#
# Build:
#   docker build -t mkws-nh3h2 .
#
# Run (interactive shell, with the current folder mounted so outputs persist):
#   docker run -it --name mkws -v ${PWD}:/work mkws-nh3h2
#
# Inside the container, run the whole pipeline:
#   python flame_nh3_h2.py     # solve all flames -> results/results.csv
#   python plots.py            # generate the figures -> figures/
# ---------------------------------------------------------------------------
FROM mambaorg/micromamba:1.5.8

WORKDIR /work

# Install Cantera and the scientific Python stack from conda-forge.
RUN micromamba install -y -n base -c conda-forge \
        python=3.11 \
        cantera \
        numpy \
        pandas \
        matplotlib && \
    micromamba clean --all --yes

# Make the base environment active for all subsequent commands.
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Copy the project source into the image.
COPY . /work

# Default to an interactive shell so the user can run the scripts manually.
CMD ["/bin/bash"]
