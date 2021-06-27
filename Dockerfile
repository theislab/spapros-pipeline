FROM continuumio/miniconda:4.7.12

# Install utilities
RUN apt update

RUN apt install vim \
                hostname \
                procps -y

# Install the conda environment
COPY environment.yml .
RUN conda env create -f environment.yml && conda clean -a

# Activate the environment
RUN echo "source activate spapros" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name spapros > spapros_environment.yml
