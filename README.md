RIC-contacts
============
Prediction of RNA-RNA contacts from RIC-seq data. Developed by Sergei Margasyuk (smargasyuk@gmail.com) and Dmitri Pervouchine (pervouchine@gmail.com).

Description
===========
This package contains a pipeline for prediction of RNA-RNA contacts from RIC-seq data (protocol from Cai et al, Nature 2020 582(7812):432-437). 
The input are two bioreplicates of rRNA depleted RIC-seq samples and two bioreplicates of rRNA depleted total RNA-seq. The output is the list of 
contacts and their respective read counts. The description of data files is provided in config.yaml. The output is provided in data/hg19/contacts. 

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

The pipeline requires Conda package manager. One can be found here https://docs.conda.io/en/latest/miniconda.html#linux-installers.

Create RIC-contacts environment by typing
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

The script will download a toy dataset (fastq files, genome, and genome coontations confined to first 100MB of chr1), unpack, update the 
config file, and execute the pipeline. The output files in data/hg19/contacts will be compared to those provided in the archive.


Run on fill RIC-seq data
========================

Download the RIC-seq files from GEO repository GSE127188 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127188). Download the control 
RNA-seq files from ENCODE consortium webpage.

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
As soon as you get these files, you can run
```
STAR --runMode genomeGenerate --genomeDir <path to your local index dir> --genomeFastaFiles hg38.chromFa.fa --sjdbGTFfile gencode.v39.annotation.gtf
```

Finally, you will need to update your config.yaml file accordingly.





