# The goal of this Snakemake file was to try and extract the assembly graph neighborhood around the samples of interest such that we could assemble a viral genome.
# This approach didn't work, either because:
# 1) sample were too shallowly sequenced so the viral genome of interest didn't have enough coverage
# 2) samples were too deeply sequenced and the neighborhood around the capsid of interest was too complex to pull anything out of.

# The samples were selected from notebooks/20230428-blast-has-virus.ipynb.
# They are the samples from the capsid BLAST results whose reads also had:
# 1) viral matches against NCBI nt & no human or vector matches
# 2) no matches against NCBI nt at all

SAMPLES = ['SRR1779416', 'SRR1779420', 'SRR1779422', 'SRR1779424', 'SRR8750801', 'SRR1778915', 
           'SRR1779200', 'SRR8750453', 'SRR8750456', 'SRR8750473', 'SRR8750491', 'SRR8750734',
           'SRR14788345', 'SRR14862884', 'SRR14999724', 'SRR14862871']
LINEAGES = ['contam']
KSIZES = [21]

# hard-bind sample wildcard with VOGs id'd in those samples
# these were made by interpretting the capsid BLAST results.
SAMPLE_VOG = ['SRR14788345-1891718.YP_007346963.1', 'SRR14862871-1891718.YP_007346963.1', 
              'SRR14862884-1891718.YP_007346963.1', 'SRR14999724-1891718.YP_007346963.1', 
              'SRR8750801-1891718.YP_007346963.1', 'SRR1778915-1277649.YP_007354884.1', 
              'SRR1779200-1277649.YP_007354884.1', 'SRR14862871-10798.YP_004928146.1', 
              'SRR8750456-10617.NP_040895.1', 'SRR8750473-10617.NP_040895.1', 'SRR8750734-10519.AP_000545.1']

rule all:
    input:
        expand("outputs/sourmash_gather_abundtrim/{sample}_k{ksize}.csv", sample = SAMPLES, ksize = KSIZES),
        expand("outputs/capsid_blast_sgc/all_outputs/{sample_vog}.fasta.cdbg_ids.reads.gz", sample_vog = SAMPLE_VOG)

rule download_runs:
    output:
        r1 = "inputs/raw/{sample}_pass_1.fastq.gz",
        r2 = "inputs/raw/{sample}_pass_2.fastq.gz"
    conda: 'envs/sratools.yml'
    benchmark: "benchmarks/download_runs/{sample}.tsv"
    shell:'''
    fastq-dump --gzip --outdir inputs/raw --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip {wildcards.sample}
    '''

rule fastp:
    input: 
        r1 = "inputs/raw/{sample}_pass_1.fastq.gz",
        r2 = "inputs/raw/{sample}_pass_2.fastq.gz"
    output: 
        r1 = 'outputs/fastp/{sample}_1.fastp.fq.gz',
        r2 = 'outputs/fastp/{sample}_2.fastp.fq.gz',
        json = 'outputs/fastp/{sample}.json',
        html = "outputs/fastp/{sample}.html"
    conda: 'envs/fastp.yml'
    benchmark: "benchmarks/fastp/{sample}.tsv"
    shell:'''
    fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -q 10 -h {output.html} -j {output.json} -l 21 -c --dedup 
    '''

rule interleave_reads:
    input:
        r1 = 'outputs/fastp/{sample}_1.fastp.fq.gz',
        r2 = 'outputs/fastp/{sample}_2.fastp.fq.gz',
    output: 'outputs/interleave/{sample}.fq.gz'
    conda: "envs/bbmap.yml"
    benchmark: "benchmarks/interleave/{sample}.tsv"
    shell:'''
    reformat.sh in1={input.r1} in2={input.r2} out={output}
    '''

rule kmer_trim:
    input:  "outputs/interleave/{sample}.fq.gz"
    output: "outputs/abundtrim/{sample}.abundtrim.fq.gz"
    benchmark: "benchmarks/kmer_trim/{sample}.tsv"
    conda: 'envs/khmer.yml'
    shell:'''
    trim-low-abund.py --gzip -C 3 -Z 18 -M 30e9 -V {input} -o {output}
    '''

rule make_spacegraphcats_config:
    """
    for spacegraphcats, which orchestrates assembly graph queries, we use k-mer size 31 since this is the only validated k-size for this type of method.
    """
    input: query = "outputs/capsid_blast_pident90/{sample}-{vog}.fasta"
    output: conf = "outputs/capsid_blast_sgc_conf/{sample}-{vog}.config"
    run:
        config_template = """\
catlas_base: {wildcards.sample}
input_sequences:
- outputs/abundtrim/{wildcards.sample}.abundtrim.fq.gz
ksize: 31
radius: 1
paired_reads: false
search:
- {input.query}
"""
        with open(output.conf, 'wt') as fp:
            fp.write(config_template.format(sample=wildcards.sample, query=input.query))


rule spacegraphcats:
    input:
        fq = "outputs/interleave/{sample}.fq.gz",
        query = "outputs/capsid_blast_pident90/{sample}-{vog}.fasta",
        conf = "outputs/capsid_blast_sgc_conf/{sample}-{vog}.config"
    output: "outputs/capsid_blast_sgc/all_outputs/{sample}-{vog}.fasta.cdbg_ids.reads.gz"
    params: outdir = "outputs/capsid_blast_sgc/"
    conda: "envs/spacegraphcats.yml"
    shell:'''
    mkdir -p outputs/capsid_blast_sgc/all_outputs/
    python -m spacegraphcats run {input.conf} extract_reads --nolock --outdir={params.outdir} --rerun-incomplete 
    cp {params.outdir}{wildcards.sample}_k31_r1_search_oh0/{wildcards.sample}-{wildcards.vog}.fasta.cdbg_ids.reads.gz {output}
    '''

############################################################
## tmp evaluate to see how much more cleaning we need to do
#
# while with spacegraphcats we use a k-mer size of 31, with
# sourmash we use a k-mer size of 21 because it's more
# permissive. This means that even if the contaminant in our
# sample is similar to but different from what we have in our
# database, we'll still detect it. This is ideal for a contam
# screen.
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
    input: "outputs/abundtrim/{sample}.abundtrim.fq.gz"
    output: "outputs/sourmash_sketch_abundtrim/{sample}.sig"
    benchmark: "benchmarks/sourmash_sketch_abundtrim/{sample}.tsv"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sketch dna -p k=21,k=31,k=51,abund,scaled=1000 -o {output} --merge {wildcards.sample}_abundtrim {input}
    '''

rule sourmash_gather:
    input:
        sig="outputs/sourmash_sketch_abundtrim/{sample}.sig",
        db=expand("inputs/sourmash_databases/genbank-2022.03-{lineage}-k{ksize}-scaled1k-cover.zip", lineage = LINEAGES),
    output: "outputs/sourmash_gather_abundtrim/{sample}_k{ksize}.csv"
    benchmark: "benchmarks/sourmash_gather_abundtrim/{sample}.tsv"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash gather -k {wildcards.ksize} -o {output} {input.sig} {input.db}
    '''
