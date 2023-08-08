# This goal of this pipeline was to iteratively remove human sequences from sequencing samples.
# We then analyzed the leftover fractions to determine what non-human were comprised of, with the goal of finding viral sequences.
# While the approach worked, it was computationally expensive to run.

import pandas as pd

metadata = pd.read_csv("inputs/concat_accessions_brain_03302023.csv", header = 0)
SAMPLES = metadata['Run'].unique().tolist()
LINEAGES = ['bacteria', 'contam', 'archaea', 'fungi', 'protozoa']

rule all:
    input:
        expand("outputs/bbduk_sourmash_gather/{sample}.csv", sample = SAMPLES),
	expand("outputs/megahit/{sample}.contigs.fa", sample = SAMPLES),
        expand('outputs/megahit_bwa/{sample}.stat', sample = SAMPLES),
        expand('outputs/capsid_blast_lca/{sample}_lca.tsv', sample = SAMPLES)

#########################################################
## Prepare human reference transcriptome
#########################################################

rule download_salmon_index:
    output: "inputs/references/salmon_sa_index.tgz"
    shell:'''
    curl -JLo {output} http://refgenomes.databio.org/v3/assets/archive/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/salmon_sa_index?tag=default
    '''

rule decompress_salmon_index:
    input: "inputs/references/salmon_sa_index.tgz"
    output: "inputs/references/salmon_index/default/info.json"
    params: outdir = "inputs/references/salmon_index/"
    shell:'''
    mkdir -p inputs/references/salmon_index
    tar xvzf {input} --directory {params.outdir} 
    '''

##########################################################
## Bait out human and low quality reads
##########################################################

rule download_runs:
    output:
        r1 = temp("inputs/raw/{sample}_pass_1.fastq.gz"),
        r2 = temp("inputs/raw/{sample}_pass_2.fastq.gz")
    conda: 'envs/sratools.yml'
    benchmark: "benchmarks/download_runs/{sample}.tsv"
    shell:'''
    fastq-dump --gzip --outdir inputs/raw --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip {wildcards.sample}
    '''

#rule download_sequencing_data:
# deal with this later -- didn't check whether sratools is still broken on osx64, just downloaded straight from ENA for now

rule salmon:
    input:
        r1="inputs/raw/{sample}_pass_1.fastq.gz",
        r2="inputs/raw/{sample}_pass_2.fastq.gz",
        index="inputs/references/salmon_index/default/info.json"
    output: 
        unmapped="outputs/salmon/{sample}_quant/aux_info/unmapped_names.txt"
    benchmark: "benchmarks/salmon/{sample}.tsv"
    params: 
        indexdir="inputs/references/salmon_index/default/",
        outdir=lambda wildcards: "outputs/salmon/" + wildcards.sample + "_quant" 
    conda: "envs/salmon.yml"
    shell:'''
    salmon quant -i {params.indexdir} -l A -1 {input.r1} -2 {input.r2} -o {params.outdir} --validateMappings --writeUnmappedNames
    '''

rule filter_unmapped_reads:
    input: "outputs/salmon/{sample}_quant/aux_info/unmapped_names.txt"
    output: "outputs/salmon/{sample}_quant/aux_info/unmapped_names_u.txt"
    benchmark: "benchmarks/filter_unmapped_reads/{sample}.tsv"
    conda: "envs/csvtk.yml"
    shell:'''
    csvtk filter2 -d " " -f '$2=="u"' {input} | csvtk cut -f 1 -o {output}
    '''

rule extract_unmapped_reads_salmon:
    input:
        fq = "inputs/raw/{sample}_pass_{pair}.fastq.gz",
        unmapped = "outputs/salmon/{sample}_quant/aux_info/unmapped_names_u.txt"
    output: temp("outputs/salmon_unmapped/{sample}_{pair}.fastq.gz")
    conda: "envs/seqtk.yml"
    benchmark: "benchmarks/extract_unmapped_reads_salmon/{sample}_{pair}.tsv"
    shell:'''
    seqtk subseq {input.fq} {input.unmapped} | gzip > {output}
    '''

rule fastp:
    input: 
        r1 = "outputs/salmon_unmapped/{sample}_1.fastq.gz",
        r2 = "outputs/salmon_unmapped/{sample}_2.fastq.gz"
    output: 
        r1 = temp('outputs/fastp/{sample}_1.fastp.fq.gz'),
        r2 = temp('outputs/fastp/{sample}_2.fastp.fq.gz'),
        json = 'outputs/fastp/{sample}.json',
        html = "outputs/fastp/{sample}.html"
    conda: 'envs/fastp.yml'
    benchmark: "benchmarks/fastp/{sample}.tsv"
    shell:'''
    fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -q 10 -h {output.html} -j {output.json} -l 21 -c --dedup --trim_poly_x --trim_poly_g
    '''

rule remove_human_host:
    """
    see make_human_extras.md for instructions on how to generate the input file inputs/references/human_extras.fasta.gz
    """
    input: 
        r1 = 'outputs/fastp/{sample}_1.fastp.fq.gz',
        r2 = 'outputs/fastp/{sample}_2.fastp.fq.gz',
        human='inputs/references/human_extras.fasta.gz'
    output:
        r1 = 'outputs/bbduk/{sample}_1.nohost.fq.gz',
        r2 = 'outputs/bbduk/{sample}_2.nohost.fq.gz',
        human_r1=temp('outputs/bbduk/{sample}_1.human.fq.gz'),
        human_r2=temp('outputs/bbduk/{sample}_2.human.fq.gz')
    conda: 'envs/bbmap.yml'
    benchmark: "benchmarks/bbduk/{sample}.tsv"
    shell:'''
    bbduk.sh -Xmx62g t={threads} in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outm={output.human_r1} outm2={output.human_r2} k=21 ref={input.human}
    '''

############################################################
## Assemble and annotate leftovers
############################################################

rule assemble:
    input:
        r1 = 'outputs/bbduk/{sample}_1.nohost.fq.gz',
        r2 = 'outputs/bbduk/{sample}_2.nohost.fq.gz',
    output: 'outputs/megahit/{sample}.contigs.fa'
    conda: 'envs/megahit.yml'
    benchmark: "benchmarks/megahit/{sample}.tsv"
    shell:'''
    mkdir -p tmp
    megahit -1 {input.r1} -2 {input.r2} --min-contig-len 100 \
        --out-dir tmp/{wildcards.sample}_megahit \
        --out-prefix {wildcards.sample} --continue
    # move the final assembly to a folder containing all assemblies
    mv tmp/{wildcards.sample}_megahit/{wildcards.sample}.contigs.fa {output}
    rm -rf tmp/{wildcards.sample}_megahit
    '''

############################################################
## tmp evaluate to see how much more cleaning we need to do
############################################################

rule download_sourmash_databases_genbank:
    input: "inputs/sourmash_databases/sourmash-database-info.csv"
    output: "inputs/sourmash_databases/genbank-2022.03-{lineage}-k21-scaled1k-cover.zip"
    run:
        sourmash_database_info = pd.read_csv(str(input[0]))
        ksize = 21 
        lineage_df = sourmash_database_info.loc[(sourmash_database_info['lineage'] == wildcards.lineage) & (sourmash_database_info['ksize'] == ksize)]
        if lineage_df is None:
            raise TypeError("'None' value provided for lineage_df. Are you sure the sourmash database info csv was not empty?")

        osf_hash = lineage_df['osf_hash'].values[0] 
        shell("curl -JLo {output} https://osf.io/{osf_hash}/download")


rule sourmash_sketch_bbduk:
    input:
        r1 = 'outputs/bbduk/{sample}_1.nohost.fq.gz',
        r2 = 'outputs/bbduk/{sample}_2.nohost.fq.gz'
    output: "outputs/bbduk_sourmash_sketch/{sample}.sig"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sketch dna -p k=21,k=31,k=51,abund,scaled=1000 -o {output} --merge {wildcards.sample}_bbduk {input}
    '''

rule sourmash_gather_bbduk:
    input:
        sig="outputs/bbduk_sourmash_sketch/{sample}.sig",
        db1=expand("inputs/sourmash_databases/genbank-2022.03-{lineage}-k21-scaled1k-cover.zip", lineage = LINEAGES),
        db2="outputs/sourmash_databases/genbank-2023.03-viral-k21.zip"
    output: "outputs/bbduk_sourmash_gather/{sample}.csv"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash gather --threshold-bp 0 -k 21 -o {output} {input.sig} {input.db1} {input.db2}
    '''

rule index_megahit:
    input: 'outputs/megahit/{sample}.contigs.fa'
    output: 'outputs/megahit/{sample}.contigs.fa.bwt'
    conda: "envs/bwa.yml"
    shell:'''
    bwa index {input}
    '''

rule bwa_mem:
    input:
        r1 = 'outputs/bbduk/{sample}_1.nohost.fq.gz',
        r2 = 'outputs/bbduk/{sample}_2.nohost.fq.gz',
        index = 'outputs/megahit/{sample}.contigs.fa.bwt',
        ref = 'outputs/megahit/{sample}.contigs.fa'
    output: 'outputs/megahit_bwa/{sample}.bam'
    conda: 'envs/bwa.yml'
    shell:'''
    bwa mem {input.ref} <(gunzip -c {input.r1}) <(gunzip -c {input.r2}) | samtools sort -o {output} -
    '''

rule bwa_evaluate:
    input:'outputs/megahit_bwa/{sample}.bam'
    output: 'outputs/megahit_bwa/{sample}.stat'
    conda: "envs/bwa.yml"
    shell:'''
    samtools view -h {input} | samtools stats > {output}
    '''

############################################################
## Search for capsid proteins
############################################################

rule search_for_capsid_proteins_in_leftover_reads:
    input:
        r1 = 'outputs/bbduk/{sample}_1.nohost.fq.gz',
        r2 = 'outputs/bbduk/{sample}_2.nohost.fq.gz',
        db = 'inputs/capsid_db/capsid_vogs_rep_seq.dmnd'
    output: temp('outputs/capsid_blast/{sample}.tsv')
    benchmark: 'benchmarks/capsid_blast/{sample}.tsv'
    conda: "envs/diamond.yml" 
    shell:'''
    diamond blastx -q {input.r1} {input.r2} -d {input.db} \
      --masking 0 --mid-sensitive -s 1 -c1 -p1 -k1 -b 0.75 \
      -f 6 qseqid qstart qend qlen qstrand sseqid sstart send \
           slen pident evalue cigar qseq_translated full_qseq full_qseq_mate > {output}
    '''

rule add_lineage_information_to_capsid_results:
    input:
        diamond='outputs/capsid_blast/{sample}.tsv',
        lca='inputs/capsid_db/capsid_vogs_rep_lca_final.tsv'
    output: lca='outputs/capsid_blast_lca/{sample}_lca.tsv'
    conda: "envs/tidyverse.yml"
    script: "scripts/add_lineage_information_to_capsid_results.R"
