import pandas as pd
import os


# read in VOG information (ID, description, and number of genomes/proteins)
# check if the file exists, and if so, read it. If not, download it.
fname = "inputs/capsid_db/allvogs.csv"
if os.path.isfile(fname):
    vogs = pd.read_csv(fname, header = 0)
else: 
    vogs = pd.read_csv("https://vogdb.org/reports/allvogs.csv?mingenomes=1&minproteins=1", header = 0)
    fdir = "inputs/capsid_db/"
    if not os.path.isdir(fdir):
        os.makedirs(fdir)
    vogs.to_csv(fname)

# filter to VOGs that have anything to do with capsid proteins
capsid_vogs = vogs[vogs.Description.str.contains("capsid", regex=False)]
# turn into a list
CAPSID_VOGS = capsid_vogs['ID'].unique().tolist() 

rule all:
    input: 
        "inputs/capsid_db/capsid_vogs_rep_seq.dmnd",
        "inputs/capsid_db/capsid_vogs_rep_lca_final.tsv"

rule download_vog_protein_sequences:
    output: "inputs/capsid_db/vog.faa.tar.gz"
    shell:'''
    curl -JLo {output} http://fileshare.csb.univie.ac.at/vog/latest/vog.faa.tar.gz
    '''

rule decompress_vog_protein_sequences:
    input: "inputs/capsid_db/vog.faa.tar.gz"
    output: expand("inputs/capsid_db/vogs/{capsid_vog}.faa", capsid_vog = CAPSID_VOGS)
    params: outdir="inputs/capsid_db/vogs/"
    shell:'''
    tar xf {input} -C {params.outdir}
    '''

rule combine_capsid_protein_sequences:
    input: expand("inputs/capsid_db/vogs/{capsid_vog}.faa", capsid_vog = CAPSID_VOGS)
    output: "inputs/capsid_db/capsid_vogs.faa"
    shell:'''
    cat {input} > {output}
    '''

rule cluster_capsid_protein_sequences:
    input: "inputs/capsid_db/capsid_vogs.faa"
    output:
        "inputs/capsid_db/capsid_vogs_rep_seq.fasta",
        "inputs/capsid_db/capsid_vogs_cluster.tsv"
    conda: "envs/mmseqs2.yml"
    params: outprefix="inputs/capsid_db/capsid_vogs"
    threads: 1
    shell:'''
    mkdir -p tmp_mmseqs2
    mmseqs easy-linclust {input} {params.outprefix} tmp_mmseqs2 --min-seq-id 0.9 -c 0.9 --similarity-type 2 --cov-mode 1 --threads {threads}
    '''

rule build_diamond_db:
    input: "inputs/capsid_db/capsid_vogs_rep_seq.fasta"
    output: "inputs/capsid_db/capsid_vogs_rep_seq.dmnd"
    params: outprefix="inputs/capsid_db/capsid_vogs_rep_seq"
    conda: "envs/diamond.yml"
    shell:'''
    diamond makedb --in {input} -d {params.outprefix}
    '''

rule download_taxdump:
    output: "inputs/taxdump.tar.gz"
    shell:'''
    curl -JLo {output} https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
    '''

rule decompress_taxdump:
    input: "inputs/taxdump.tar.gz"
    output: "inputs/taxdump/nodes.dmp"
    shell:'''
    tar xf {input} -C inputs/taxdump/
    '''

rule get_taxid_for_each_protein:
    input: "inputs/capsid_db/capsid_vogs_cluster.tsv"
    output: "inputs/capsid_db/capsid_vogs_rep_taxid.tsv"
    conda: "envs/csvtk.yml"
    shell:'''
    cat {input} | csvtk mutate -H -t -f 2 -n taxid -p "^(.+?)\." | csvtk fold -H -f 1 -v 3 -s ";" -t -o {output}
    '''

rule get_lca_taxid_for_each_cluster:
    input: 
        tsv="inputs/capsid_db/capsid_vogs_rep_taxid.tsv",
        dmp="inputs/taxdump/nodes.dmp"
    output: "inputs/capsid_db/capsid_vogs_rep_taxid_lca.tsv"
    params: datadir="inputs/taxdump/"
    conda: "envs/taxonkit.yml" # needs to be 0.14.2, which isn't on conda yet
    # in the meantime, the executable can be downloaded from this url: https://github.com/shenwei356/taxonkit/files/11073880/taxonkit_linux_amd64.tar.gz
    shell:'''
    taxonkit lca --data-dir {params.datadir} -i 2 -s ";" -o {output} {input.tsv}
    '''

rule reformat_lca_taxid_to_lineage:
    input: 
        tsv="inputs/capsid_db/capsid_vogs_rep_taxid_lca.tsv",
        dmp="inputs/taxdump/nodes.dmp"
    output: "inputs/capsid_db/capsid_vogs_rep_lca.tsv"
    params: datadir="inputs/taxdump/"
    conda: "envs/taxonkit.yml"
    shell:'''
    taxonkit reformat -I 3 -f "{{k}};{{K}};{{p}};{{c}};{{o}};{{f}};{{g}};{{s}};{{t}}" -F --data-dir {params.datadir} -t -o {output} {input.tsv}
    '''

rule add_header:
    input: "inputs/capsid_db/capsid_vogs_rep_lca.tsv" 
    output: "inputs/capsid_db/capsid_vogs_rep_lca_final.tsv"
    conda: "envs/csvtk.yml"
    shell:'''
    csvtk add-header -t -I -H -n rep,taxid,lca_taxid,lca_lineage_named,lca_lineage_taxid -o {output} {input}
    '''
