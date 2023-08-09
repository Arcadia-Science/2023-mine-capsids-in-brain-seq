# The goal of this pipeline was to quickly discover sequencing samples that potentially contained viral reads.
# We tried two approaches.
# The first was an exact k-mer matching approach between a sample of interest and databases of genomes.
# The second was a diamond BLASTP search against a database of capsid amino acid sequences.
# The capsid BLAST search was fast and accurate after filtering.
# Sourmash was too slow to be used as a screening tool for this use case.

import pandas as pd

metadata = pd.read_csv("inputs/concat_accessions_brain_03302023_small.tsv", header = 0, sep = "\t")
SAMPLES = metadata['Run'].unique().tolist()
LINEAGES = ['bacteria', 'contam', 'archaea', 'fungi', 'protozoa']
KSIZES = [21] 

rule all:
    input:
        expand("outputs/sourmash_gather_raw/{sample}_k{ksize}.csv", sample = SAMPLES, ksize = KSIZES),
        expand('outputs/capsid_blast_raw_lca/{sample}_lca.tsv', sample = SAMPLES)

rule download_runs:
    output:
        r1 = "inputs/raw/{sample}_pass_1.fastq.gz",
        r2 = "inputs/raw/{sample}_pass_2.fastq.gz"
    conda: 'envs/sratools.yml'
    benchmark: "benchmarks/download_runs/{sample}.tsv"
    shell:'''
    fastq-dump --gzip --outdir inputs/raw --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip {wildcards.sample}
    '''

############################################################
## Assess taxonomic breakdown of input samples
############################################################

rule download_sourmash_databases_genbank:
    input: "inputs/sourmash_databases/sourmash-database-info.csv"
    output: "inputs/sourmash_databases/genbank-2022.03-{lineage}-k{ksize}-scaled1k-cover.zip"
    run:
        sourmash_database_info = pd.read_csv(str(input[0]))
        lineage_df = sourmash_database_info.loc[(sourmash_database_info['lineage'] == wildcards.lineage) & (sourmash_database_info['ksize'] == wildcards.ksize)]
        if lineage_df is None:
            raise TypeError("'None' value provided for lineage_df. Are you sure the sourmash database info csv was not empty?")

        osf_hash = lineage_df['osf_hash'].values[0] 
        shell("curl -JLo {output} https://osf.io/{osf_hash}/download")


rule sourmash_sketch:
    input:
        r1 = "inputs/raw/{sample}_pass_1.fastq.gz",
        r2 = "inputs/raw/{sample}_pass_2.fastq.gz"
    output: "outputs/sourmash_sketch_raw/{sample}.sig"
    benchmark: "benchmarks/sourmash_sketch_raw/{sample}.tsv"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sketch dna -p k=21,k=31,k=51,abund,scaled=1000 -o {output} --merge {wildcards.sample}_raw {input}
    '''

rule sourmash_gather:
    input:
        sig="outputs/sourmash_sketch_raw/{sample}.sig",
        db1=expand("inputs/sourmash_databases/genbank-2022.03-{lineage}-k{ksize}-scaled1k-cover.zip", lineage = LINEAGES),
        db2="outputs/sourmash_databases/genbank-2023.03-viral-k{ksize}.zip"
    output: "outputs/sourmash_gather_raw/{sample}_k{ksize}.csv"
    benchmark: "benchmarks/sourmash_gather_raw/{sample}.tsv"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash gather --threshold-bp 0 -k {wildcards.ksize} -o {output} {input.sig} {input.db1} {input.db2}
    '''

############################################################
## Search for capsid proteins
############################################################

rule search_for_capsid_proteins_in_leftover_reads:
    input:
        r1 = "inputs/raw/{sample}_pass_1.fastq.gz",
        r2 = "inputs/raw/{sample}_pass_2.fastq.gz",
        db = 'inputs/capsid_db/capsid_vogs_rep_seq.dmnd'
    output: temp('outputs/capsid_blast_raw/{sample}.tsv')
    benchmark: 'benchmarks/capsid_blast_raw/{sample}.tsv'
    conda: "envs/diamond.yml" 
    shell:'''
    diamond blastx -q {input.r1} {input.r2} -d {input.db} \
      --masking 0 --mid-sensitive -s 1 -c1 -p1 -k1 -b 0.75 \
      -f 6 qseqid qstart qend qlen qstrand sseqid sstart send \
           slen pident evalue cigar qseq_translated full_qseq full_qseq_mate > {output}
    '''

rule add_lineage_information_to_capsid_results:
    input:
        diamond='outputs/capsid_blast_raw/{sample}.tsv',
        lca='inputs/capsid_db/capsid_vogs_rep_lca_final.tsv'
    output: lca='outputs/capsid_blast_raw_lca/{sample}_lca.tsv'
    conda: "envs/tidyverse.yml"
    script: "scripts/add_lineage_information_to_capsid_results.R"
