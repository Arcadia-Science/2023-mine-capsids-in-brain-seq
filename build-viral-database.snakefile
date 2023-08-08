GROUPS = ["viral"]
KSIZES = [21, 31, 51]

rule all:
    input:
        expand("outputs/sourmash_database_covers/genbank-2023.03-viral-k{ksize}-scaled10k-cover.zip", ksize = KSIZES),
        "outputs/sourmash_database_lineages/genbank-2023.03-viral.lineages.csv"

###################################################################
## Download, repeat mask, and sketch genomes
###################################################################

rule download_viral_assembly_summary:
    output: "inputs/genbank/viral_assembly_summary.txt"
    shell:'''
    curl -JLo {output} https://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/assembly_summary.txt
    '''

rule grab_viral_taxids:
    input:"inputs/genbank/viral_assembly_summary.txt"
    output:"inputs/genbank/viral_assembly_summary_taxids_uniq.txt"
    conda: "envs/csvtk.yml"
    shell:'''
    cut -f6 {input} | tail -n +2 | csvtk freq -f taxid | csvtk cut -f 1 | csvtk del-header -o {output}
    '''

checkpoint download_reference_genomes:
    input: "inputs/genbank/viral_assembly_summary_taxids_uniq.txt"
    output: directory("inputs/genbank/viral/")
    conda: "envs/ncbi-genome-download.yml"
    params: outdir= "inputs/genbank/viral/"
    benchmark: "benchmarks/download_reference_genomes/viral.tsv"
    threads: 7
    shell:'''
    while read inline
    do 
        ncbi-genome-download viral --taxids ${{inline}} -s genbank -o {params.outdir} --flat-output -F fasta -p {threads} -P 
    done < {input}
    '''

rule sketch_genome:
   input: "inputs/genbank/viral/{accession}_genomic.fna.gz"
   output: "outputs/sourmash_sketch/viral/{accession}_genomic.sig.gz"
   conda: "envs/sourmash.yml"
   shell:'''
   sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --name {wildcards.accession} -o {output} {input}
   '''

###################################################################
## Create databases
###################################################################

def checkpoint_download_reference_genomes(wildcards):
    # expand checkpoint to get accession values, and place them in the final file name that uses that wildcard
    # checkpoint_output encodes the output dir from the checkpoint rule. 
    checkpoint_output = checkpoints.download_reference_genomes.get(**wildcards).output[0]    
    file_names = expand("outputs/sourmash_sketch/viral/{accession}_genomic.sig.gz",
                        accession = glob_wildcards(os.path.join(checkpoint_output, "{accession}_genomic.fna.gz")).accession)
    return file_names

rule combine_to_filelist:
   input: checkpoint_download_reference_genomes
   output: "outputs/sourmash_sketch/viral_sig_paths.txt"
   run:
       with open(str(output[0]), 'w') as f:
           for line in input:
               f.write(f"{line}\n")


rule build_database:
   input: "outputs/sourmash_sketch/viral_sig_paths.txt"
   output: "outputs/sourmash_databases/genbank-2023.03-viral-k{ksize}.zip"
   conda: "envs/sourmash.yml"
   shell:'''
   sourmash sig cat -k {wildcards.ksize} --from-file {input} -o {output}
   '''

rule build_cover_databases:
    input: "outputs/sourmash_databases/genbank-2023.03-viral-k{ksize}.zip"
    output: "outputs/sourmash_database_covers/genbank-2023.03-viral-k{ksize}-scaled10k-cover.zip"
    conda: "envs/sourmash.yml"
    shell:'''
    scripts/make-db-cover.py {input} -o {output}
    '''

###################################################################
## Create taxonomy sheets to accompany databases
###################################################################

rule create_csv_describing_contents_of_database:
    input: expand("outputs/sourmash_databases/genbank-2023.03-viral-k{ksize}.zip", ksize = KSIZES)
    output: "outputs/sourmash_databases/genbank-2023.03-viral.csv"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sig describe --csv {output} {input[0]} 
    '''

rule download_assembly_to_taxid_information:
    output: "inputs/genbank/assembly_summary_genbank.txt"
    shell:'''
    curl -JLo {output} https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt
    '''

rule get_taxid_for_each_genome_accession:
    input:
        assembly_summary="inputs/genbank/assembly_summary_genbank.txt",
        sig_describe="outputs/sourmash_databases/genbank-2023.03-viral.csv"
    output: taxids="outputs/sourmash_database_lineages/genbank-2023.03-viral-taxids.tsv"
    conda: "envs/tidyverse.yml"
    script: "scripts/get_taxid_for_each_genome_accession.R"

rule download_taxdump:
    output: "inputs/genbank/taxdump.tar.gz"
    shell:'''
    curl -JLo {output} https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
    '''

rule decompress_taxdump:
    input: "inputs/genbank/taxdump.tar.gz"
    output: "inputs/genbank/taxdump/nodes.dmp"
    params: outdir="inputs/genbank/taxdump/"
    shell:'''
    tar xf {input} -C {params.outdir}
    '''

rule convert_taxid_to_lineage:
    input:
        taxids="outputs/sourmash_database_lineages/genbank-2023.03-viral-taxids.tsv",
        dmp="inputs/genbank/taxdump/nodes.dmp"
    output: "outputs/sourmash_database_lineages/genbank-2023.03-viral-lineages.tsv"
    params: taxdir="inputs/genbank/taxdump/"
    conda: "envs/taxonkit.yml"
    shell:'''
    taxonkit reformat -I 2 -f "{{k}};{{p}};{{c}};{{o}};{{f}};{{g}};{{s}};{{t}}" --data-dir {params.taxdir} -o {output} {input.taxids}
    '''

rule reformat_for_sourmash_taxonomy:
    input: tsv = "outputs/sourmash_database_lineages/genbank-2023.03-viral-lineages.tsv"
    output: csv = "outputs/sourmash_database_lineages/genbank-2023.03-viral.lineages.csv"
    conda: "envs/tidyverse.yml"
    script: "scripts/reformat_for_sourmash_taxonomy.R"
