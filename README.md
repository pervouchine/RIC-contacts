RIC-contacts
==========
Prediction of RNA-RNA contacts from RIC-seq data. Developed by Sergei Margasyuk (smargasyuk@gmail.com) and Dmitri Pervouchine (pervouchine@gmail.com).

This package contains a pipeline for prediction of RNA-RNA contacts from RIC-seq data (Cai et al, Nature 2020 582(7812):432-437).

Install
==========

Clone this repository by typing

```
git clone https://github.com/pervouchine/RIC-contacts
```

The pipeline requires Conda package manager. One can be found here https://docs.conda.io/en/latest/miniconda.html#linux-installers.

To create RIC-contacts environment, type

```
conda create -n RIC-contacts
conda activate RIC-contacts
```

Then, add channels and install modules by typing

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install star
conda install samtools
conda install bedopts
conda install snakemake
```

To make a test run, type

```make test```
