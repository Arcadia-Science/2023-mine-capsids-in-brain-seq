# The goal of the snakefile is to bait out viral reads via mapping to see if:
# 1) the viral genome we think we detected is actually in the sample
# 2) if there is enough reads/information to assemble a genome or a capsid from the reads
#
# This file relies on the abundtrim and fastp outputs from 02a_assembly_graph_approach.snakefile.
# In the future however, we think it would be appropriate to directly map raw reads instead of working with quality controlled reads.
#
# these samples were selected from notebooks/20230428-blast-has-virus.ipynb.
# They are the samples from the capsid BLAST results whose reads also had:
# 1) viral matches against NCBI nt & no human or vector matches
# 2) no matches against NCBI nt at all

SAMPLES = ['SRR8750801', 'SRR1778915',  'SRR1779200', 'SRR8750456', 'SRR8750473', 'SRR8750734',
           'SRR14788345', 'SRR14862884', 'SRR14999724', 'SRR14862871']
SEQTYPE = ['nucleotides', 'proteins']

rule all:
    input:
        expand("outputs/mapping/{sample}_{seqtype}_pe.coverage", sample = SAMPLES, seqtype = SEQTYPE),
        expand("outputs/mapping/{sample}_proteins.tsv", sample = SAMPLES)

########################################################################
## Nucleotide
########################################################################

rule download_references:
    output: "inputs/references/nucleotides.fna"
    shell:'''
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/597/305/GCA_002597305.1_ASM259730v1/GCA_002597305.1_ASM259730v1_genomic.fna.gz # Pbunalikevirus phiFenriz phage
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz # COVID
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/904/055/GCF_000904055.1_ViralProj186434/GCF_000904055.1_ViralProj186434_genomic.fna.gz # ALREADY HAD: DELTAPOLYOMA VOG
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/864/845/GCF_000864845.1_ViralProj15492/GCF_000864845.1_ViralProj15492_genomic.fna.gz # ALREADY HAD: Gammapapillomavirus
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/903/895/GCF_000903895.1_ViralProj185188/GCF_000903895.1_ViralProj185188_genomic.fna.gz # NEW: ALPHAPOLYOMA VOG
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/839/645/GCF_000839645.1_ViralProj14090/GCF_000839645.1_ViralProj14090_genomic.fna.gz # NEW: PARVOVIRUS B19 
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/859/485/GCF_000859485.1_ViralProj15114/GCF_000859485.1_ViralProj15114_genomic.fna.gz # NEW: MASTADENO
    gunzip *fna.gz
    esearch -db nucleotide -query "MH777161.1" | efetch -format fasta > MH777161.1.fna
    cat *fna > nucleotides.fna
    '''

rule index_nucleotides:
    input: "inputs/references/nucleotides.fna"
    output: "inputs/references/nucleotides.1.bt2"
    conda: "envs/bowtie2.yml"
    shell:'''
    bowtie2-build {input} inputs/references/nucleotides
    '''

rule bowtie2:
    input: 
        fq="outputs/abundtrim/{sample}.abundtrim.fq.gz",
        index="inputs/references/nucleotides.1.bt2"
    output: "outputs/mapping/{sample}_nucleotides.bam"
    conda: "envs/bowtie2.yml"
    threads: 8
    shell:'''
    bowtie2 -p {threads} -x nucleotides --no-unal --very-sensitive-local --interleaved {input.fq} | samtools view -bSF4 - > {output}
    '''

###############################################################
## proteins
###############################################################

rule download_proteins:
    output: "inputs/references/proteins.faa"
    shell:'''
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/597/305/GCA_002597305.1_ASM259730v1/GCA_002597305.1_ASM259730v1_protein.faa.gz         # Pbunalikevirus phiFenriz phage
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_protein.faa.gz         # COVID
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/904/055/GCF_000904055.1_ViralProj186434/GCF_000904055.1_ViralProj186434_protein.faa.gz # DELTAPOLYOMA VOG
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/864/845/GCF_000864845.1_ViralProj15492/GCF_000864845.1_ViralProj15492_protein.faa.gz   # Gammapapillomavirus
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/903/895/GCF_000903895.1_ViralProj185188/GCF_000903895.1_ViralProj185188_protein.faa.gz # ALPHAPOLYOMA VOG
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/839/645/GCF_000839645.1_ViralProj14090/GCF_000839645.1_ViralProj14090_protein.faa.gz   # PARVOVIRUS B19 
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/859/485/GCF_000859485.1_ViralProj15114/GCF_000859485.1_ViralProj15114_protein.faa.gz   # MASTADENO
    gunzip *faa.gz
    esearch -db protein -query "AYA93430.1" | efetch -format fasta > AYA93430.1.faa
    esearch -db protein -query "AYA93431.1" | efetch -format fasta > AYA93431.1.faa
    esearch -db protein -query "AYA93432.1" | efetch -format fasta > AYA93432.1.faa
    esearch -db protein -query "AYA93433.1" | efetch -format fasta > AYA93433.1.faa
    esearch -db protein -query "AYA93434.1" | efetch -format fasta > AYA93434.1.faa
    esearch -db protein -query "AYA93435.1" | efetch -format fasta > AYA93435.1.faa
    cat *faa > {output}
    '''

rule index_proteins:
    input: "inputs/references/proteins.faa"
    output: "inputs/references/proteins.faa.bwt"
    conda: "envs/paladin.yml"
    shell: '''
    paladin index -r3 {input}
    '''

rule run_paladin:
    input:
        fq = "outputs/abundtrim/{sample}.abundtrim.fq.gz",
        index = "inputs/references/proteins.faa.bwt"
    output: "outputs/mapping/{sample}_proteins.bam"
    conda: "envs/paladin.yml"
    threads: 8
    shell:'''
    # filter to only mapped reads, but not only PE reads
    paladin align -t {threads} proteins.faa {input.fq} | samtools view -b -F 4 - > {output}
    '''

##############################################################
## process alignments
##############################################################

rule filter_bam:
    input: "outputs/mapping/{sample}_{seqtype}.bam"
    output: "outputs/mapping/{sample}_{seqtype}_pe.bam"
    conda: "envs/samtools.yml"
    shell:'''
    samtools view -b -f 0x02 {input} > {output}
    '''

rule index_bam:
    input: "outputs/mapping/{sample}_{seqtype}_pe.bam"
    output: "outputs/mapping/{sample}_{seqtype}_pe.sorted.bam"
    conda: "envs/samtools.yml"
    shell:'''
    samtools sort {input} > {output} && samtools index {output}
    '''

rule stat_bam:
    input: "outputs/mapping/{sample}_{seqtype}_pe.sorted.bam"
    output: 
        stats="outputs/mapping/{sample}_{seqtype}_pe.stats",
        flagstat="outputs/mapping/{sample}_{seqtype}_pe.flagstat",
        depth="outputs/mapping/{sample}_{seqtype}_pe.depth",
        coverage="outputs/mapping/{sample}_{seqtype}_pe.coverage",
    conda: "envs/samtools.yml"
    shell:'''
    samtools stats {input} > {output.stats}
    samtools flagstat {input} > {output.flagstat}
    samtools depth {input} > {output.depth}
    samtools coverage {input} > {output.coverage}
    '''

##############################################################
## Diamond protein map instead
##############################################################

rule make_diamond_db:
    input: "inputs/references/proteins.faa"
    output: "inputs/references/proteins.dmnd"
    conda: "envs/diamond.yml"
    params: outprefix = "inputs/references/proteins"
    shell:'''
    diamond makedb --in {input} -d {params.outprefix}
    '''

rule diamond_protein_map:
    input:
        r1 = "outputs/fastp/{sample}_1.fastp.fq.gz",
        r2 = "outputs/fastp/{sample}_2.fastp.fq.gz",
        db = "inputs/references/proteins.dmnd"
    output: 'outputs/mapping/{sample}_proteins.tsv'
    conda: "envs/diamond.yml"
    shell:'''
    diamond blastx -q {input.r1} {input.r2} -d {input.db} \
      --masking 0 --mid-sensitive -s 1 -c1 -p1 -k1 -b 0.75 \
      -f 6 qseqid qstart qend qlen qstrand sseqid sstart send \
           slen pident evalue cigar qseq_translated full_qseq full_qseq_mate > {output}
    '''
