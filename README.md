RIC-contacts
============
Prediction of RNA-RNA contacts from RIC-seq data. Developed by Sergei Margasyuk (smargasyuk@gmail.com) and Dmitri Pervouchine (pervouchine@gmail.com).

Description
===========
This package contains a pipeline for prediction of RNA-RNA contacts from RIC-seq data (protocol from Cai et al, Nature 2020 582(7812):432-437). 
The input files are in fastq format, currently two bioreplicates of rRNA depleted RIC-seq samples and two bioreplicates of rRNA depleted total 
RNA-seq. The output is the list of contacts and their respective read counts in tsv format (columns 1-4 and 5-8 are the contacting coordinates, 
column 9 is read count). The description of input files must be provided in config.yaml. The output will be placed in data/hg19/contacts. 

Hardware/software requirements
==============================
  * x86-64 compatible processors
  * 64 bit Linux or Mac OS X

Installation
============
Clone this repository by typing
```
git clone https://github.com/pervouchine/RIC-contacts
```

The pipeline requires Conda package manager. One for Linux can be found here https://docs.conda.io/en/latest/miniconda.html#linux-installers.

Create a RIC-contacts environment by typing
```
conda create -n RIC-contacts
conda activate RIC-contacts
```

Add channels and install modules by typing
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install star
conda install samtools
conda install bedopts
conda install snakemake
```

Test run
============
To make a test run, type

```
cd RIC-contacts
make test
```

The script will download a toy dataset (truncated fastq files, genome, and genome annotation confined to the first 100MB of chr1), unpack, 
update the config file, and execute the pipeline. The output files in data/hg19/contacts will be compared to those provided in the archive.


Run on full RIC-seq data
========================

Download the RIC-seq files from GEO repository GSE127188 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127188). Download the control 
RNA-seq files from ENCODE consortium webpage (https://www.encodeproject.org/)

The files are as follows:
```
  RNASeq_HeLa_total_rep1:
    - fastq/ENCFF000FOM.fastq
    - fastq/ENCFF000FOV.fastq
  RNASeq_HeLa_total_rep2:
    - fastq/ENCFF000FOK.fastq
    - fastq/ENCFF000FOY.fastq
  RIC-seq_HeLa_rRNA_depleted_rep1:
    - fastq/SRR8632820_1.fastq
    - fastq/SRR8632820_2.fastq
  RIC-seq_HeLa_rRNA_depleted_rep2:
    - fastq/SRR8632821_1.fastq
    - fastq/SRR8632821_2.fastq
```

NOTE: You will also need a full version of the human genome (we recommend GRCh38 at https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz) 
and the full version of transcript annotation from GENCODE (https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz).
As soon as you get these files, you need to create STAR indicies by running
```
STAR --runMode genomeGenerate --genomeDir <path_to_your_local_index_dir> --genomeFastaFiles hg38.chromFa.fa --sjdbGTFfile gencode.v39.annotation.gtf
```

You will also need to update your config file accordingly.

To run on a single node, type
```
snakemake --cores 4 allContacts
```

To run on HPC cluster, type (for example)
```
snakemake --cluster "qsub -d . -l mem=40gb -l nodes=1:ppn=16 -l walltime=48:00:00" --jobs 70 allContacts
```

The output will be placed in data/hg19/contacts.






