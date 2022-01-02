configfile: "config.yaml"

rule alignControlPE:
    input:
        fastqs = lambda wildcards: config['samples'][wildcards["id"]],
        genome = lambda wildcards: config['star']['genome_path']
    params:
        STAR = lambda wildcards: config['star']['program'],
        out_prefix = lambda wildcards: f"data/hg19/bam/pass1/{wildcards['id']}/",
        STAR_params = '--runMode alignReads --outSAMtype BAM SortedByCoordinate',
	runThreadN = lambda wildcards: config['star']['runThreadN']
    output:
        sj = "data/hg19/bam/pass1/{id}/SJ.out.tab"
    shell: 
        """
mkdir -p $(dirname {output})
{params.STAR} {params.STAR_params} --runThreadN {params.runThreadN} --genomeDir {input.genome} --readFilesIn {input.fastqs} --outFileNamePrefix {params.out_prefix}
"""

rule filterPass1Junctions:
    input:
        "data/hg19/bam/pass1/{id}/SJ.out.tab"
    output:
        "data/hg19/bam/pass1/{id}/SJ.out.filtered.tab"
    shell:
        """
awk -v 'OFS="\t"' '$5 == 1' {input} > {output}
"""

rule alignRICSE:
    input:
        fastq = lambda wildcards: f'{config["samples"][wildcards["id"]][int(wildcards["mate"])]}',
        junctions = expand("data/hg19/bam/pass1/{id}/SJ.out.filtered.tab", id=config['samples_control']),
        genome = lambda wildcards: config['star']['genome_path']
    params:
        STAR = lambda wildcards: config['star']['program'],
        out_prefix = lambda wildcards: f"data/hg19/bam/pass2/{wildcards['id']}_{wildcards['mate']}/",
        STAR_params = '--runMode alignReads --outSAMtype BAM SortedByCoordinate --chimOutType Junctions --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreJunctionNonGTAG  -1 --scoreGapNoncan -1 --scoreGapATAC -1 --scoreGapGCAG -1 --chimSegmentReadGapMax 3 --outFilterMatchNminOverLread 0.5 --outFilterScoreMinOverLread 0.5',
        runThreadN = lambda wildcards: config['star']['runThreadN']
    output:
        bam = "data/hg19/bam/pass2/{id}_{mate}/Aligned.sortedByCoord.out.bam",
        chm = "data/hg19/bam/pass2/{id}_{mate}/Chimeric.out.junction"
    shell: 
        """
mkdir -p $(dirname {output})
{params.STAR} {params.STAR_params} --runThreadN {params.runThreadN} --genomeDir {input.genome} --sjdbFileChrStartEnd {input.junctions} --readFilesIn {input.fastq} --outFileNamePrefix {params.out_prefix}
"""

rule mergeKnownJunctions:
    input:
        expand("data/hg19/bam/pass1/{id}/SJ.out.tab", id=config['samples_control']),
	lambda wildcards: f"{config['star']['genome_path']}sjdbList.out.tab"
    output:
        "data/hg19/all_sj.tsv"
    shell:
        """
mkdir -p $(dirname {output})
cut -f1,2,3 {input} | sort -u > {output}
"""

rule extractChimericJunctions:
    input:
        mate0 = "data/hg19/bam/pass2/{id}_0/Chimeric.out.junction",
        mate1 = "data/hg19/bam/pass2/{id}_1/Chimeric.out.junction"
    output:
        "data/hg19/junctions/{id}/Chimeric.tsv"
    shell:
        """
mkdir -p $(dirname {output})
sort -k9,9 <(cut -f1-6,10 {input.mate0} | awk -v OFS="\t" 'BEGIN{{s["+"]="-";s["-"]="+";}}{{print $4, $5, $5+1, s[$6], $1, $2, $2+1, s[$3], $7, 0}}') \
           <(cut -f1-6,10 {input.mate1} | awk -v OFS="\t" '{{print $1, $2, $2+1, $3, $4, $5, $5+1, $6, $7, 1}}') | \
           grep -v GL > {output}
"""

rule extractNeoJunctions:
    input:
        mate0 = "data/hg19/bam/pass2/{id}_0/Aligned.sortedByCoord.out.bam",
	mate1 = "data/hg19/bam/pass2/{id}_1/Aligned.sortedByCoord.out.bam",
        junctions = "data/hg19/all_sj.tsv"
    output:
        "data/hg19/junctions/{id}/Neo.tsv"
    shell:
        """
mkdir -p $(dirname {output})
sort -k9,9 <(samtools view {input.mate0} | perl scripts/neo.pl {input.junctions} 0) \
           <(samtools view {input.mate1} | perl scripts/neo.pl {input.junctions} 1) | grep -v GL > {output}        
"""

rule clusterDonorsAcceptors:
    input:
        expand("data/hg19/junctions/{id}/{jtype}.tsv", id=config['samples'].keys(),jtype=["Neo", "Chimeric"])
    output:
        donors = "data/hg19/junctions/donors.bed",
        acceptors = "data/hg19/junctions/acceptors.bed"
    params:
        radius = lambda wildcards: config['junctions_merge_radius']
    shell:
        """
cut -f1-4 {input} | sort-bed - | bedops -m --range {params.radius} - > {output.donors}
cut -f5-8 {input} | sort-bed - | bedops -m --range {params.radius} - > {output.acceptors}
"""

rule clusterJunctions:
    input:
        junctions = "data/hg19/junctions/{id}/{jtype}.tsv",
	donors = "data/hg19/junctions/donors.bed",
	acceptors = "data/hg19/junctions/acceptors.bed"
    params:
        radius = lambda wildcards: config['junctions_merge_radius']
    output:
        "data/hg19/clusters/{id}/{jtype}.tsv"
    shell:
        """
mkdir -p $(dirname {output})
sort-bed {input.junctions} | \
intersectBed -a stdin -b {input.donors} -wa -wb | cut -f5- | sort-bed - | \
intersectBed -a stdin -b {input.acceptors} -wa -wb | cut -f5- | sort -k1,1 > {output}
"""

rule extractContacts:
    input:
        "data/hg19/clusters/{id}/{jtype}.tsv",
    output:
        "data/hg19/contacts/{id}/{jtype}.tsv",
    shell:
        """
mkdir -p $(dirname {output})
awk '{{n[$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10]++}}END{{for(j in n){{print j"\t"n[j]}}}}' {input} | sort > {output}
"""

rule allContacts:
    input:
        expand("data/hg19/contacts/{id}/{jtype}.tsv", id=config['samples'].keys(), jtype=["Neo", "Chimeric"]) 















rule getRICAligned:
    input:
        expand("data/hg19/bam/pass2/{id}_{mate}/Aligned.sortedByCoord.out.bam", id=config['samples'].keys(), mate=[0,1])

