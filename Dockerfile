# Base image
FROM continuumio/anaconda3

# setup
RUN apt-get update
RUN apt-get install -y build-essential
RUN apt-get install -y vim
RUN conda update -n base -c defaults conda

# Add conda-forge channel
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# Create myenv
ARG env_name=RIC-contacts
RUN conda create -yn ${env_name} python=3.8.8

# Activate environment
ENV CONDA_DEFAULT_ENV ${env_name}

# Switch default environment
RUN echo "conda activate ${env_name}" >> ~/.bashrc
ENV PATH /opt/conda/envs/${env_name}/bin:$PATH

# Install some packages from conda-forge
RUN conda install -y star
RUN conda install -y samtools
RUN conda install -y bedtools
RUN conda install -y bedops
RUN conda install -y snakemake
