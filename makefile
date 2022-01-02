.PHONY : test

test:   data/hg19/genome/index/SA config.yaml
	snakemake --cores 1 allContacts

config.yaml data/hg19/genome/chr1part.fa data/hg19/genome/v35lift37part.gtf:
	wget http://arkuda.skoltech.ru/~dp/RICseq_toy_data_31.12.2021.tar.gz
	tar -xf RICseq_toy_data_31.12.2021.tar.gz
	perl -e '$$star=`which STAR`;chomp $$star;$$config = `cat data/config.yaml`; $$config=~s/STAR/$$star/;print $$config' >config.yaml

data/hg19/genome/index/SA: data/hg19/genome/chr1part.fa data/hg19/genome/v35lift37part.gtf
	STAR --runMode genomeGenerate --genomeDir data/hg19/genome/index/ --genomeFastaFiles data/hg19/genome/chr1part.fa --sjdbGTFfile data/hg19/genome/v35lift37part.gtf
