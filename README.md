# Mining viral sequences from human brain sequencing data via capsid searches

This repository is associated with the pub, ["A capsid-based search recovers viral sequences from human brain sequencing data"](ADD LINK DOI WHEN AVAILABLE).

## Explanation of the files in this repository

Below we provide an explanation of the files in this repository.

### Analysis code

The majority of the code in this repository is orchestrated by snakefiles.
Skip to the [Running the code in this repository](#running-the-code-in-this-repository) section of the README for a description of how to run these files.

* [`00_subtractive_approach.snakefile`](00_subtractive_approach.snakefile): This snakefile contains our original approach, which used increasingly expensive computational tools to remove human sequences from DNA and RNA seq data. Our goal was to be able to explore what was in the "leftover" fraction of reads, and to perform viral discovery on that fraction. This approach was too computationally expensive to run on hundreds of samples, so we moved to a new approach.
* [`01_capsid_blast_approach.snakefile`](01_capsid_blast_approach.snakefile): This snakefile implements a DIAMOND BLASTx search for viral capsid protein sequences. We used this to quickly discover potential viral sequences in DNA and RNA seq data. It also uses sourmash gather to detect background contamination in the raw sequencing data.
* [`02a_assembly_graph_approach.snakefile`](02a_assembly_graph_approach.snakefile): This snakefile uses spacegraphcats to perform assembly graph queries with viral capsid reads identified in `01_capsid_blast_approach.snakefile`. The goal of this snakefile was to recover reads that are nearby to capsid sequences in the assembly graph and to use those reads to assemble more complete viral genomes. However, this approach either failed because the reads were too low coverage in areas with viral capsids, leading to fragmented assembly graphs, or because the regions of the graph were too complex around the viral genomes and not resolvable to contiguous sequences. 
* [`02b_mapping_approach.snakefile`](02b_mapping_approach.snakefile): This snakefile orchestrates mapping of the reads against viral genomes to estimate coverage of the genome and whether that virus was truly present in the sequencing data.

Two snakefiles build input databases for the other workflows.
* [`build-viral-database.snakefile`](./build-viral-database.snakefile): builds a k-mer database of viral sequences (sourmash FracMinHash database).
* [`build-capsid-db.snakefile`](./build-capsid-db.snakefile): builds the viral capsid sequence database used to search sequencing data sets with diamond blastx.

The snakefiles make use of the [inputs](./inputs), [envs](./envs), and [scripts](./scripts) folders.
The [inputs](./inputs) folder contains tabular files with metadata about the samples that are analyzed by the workflows.
The [envs](./envs) folder contains conda environments for the rules in the snakefile.
The [scripts](./scripts) folder contains scripts that are run by different rules. 

#### [Notebooks](./notebooks)

This repository contains three notebooks:

* [`20230426-explore-capsid-blast-results.ipynb`](./notebooks/20230426-explore-capsid-blast-results.ipynb): explores the results of the capsid BLAST analysis undertaken in `01_capsid_blast_approach.snakefile`. It examines filtering cut offs for the capsid BLAST results and writes FASTA files of the BLAST hits.
* [`20230428-blast-has-virus.ipynb`](./notebooks/20230428-blast-has-virus.ipynb): looks at the results of BLASTing the capsid reads (FASTA) against the NCBI nucleotide (nt) database. It filters out human and vector/clone hits.
* [`20230512-interpret-mapping.ipynb`](./notebooks/20230512-interpret-mapping.ipynb): interprets the mapping undertaken in `02b_mapping_approach.snakefile`. For viral genomes that were high coverage, it produces visualizations of read depth across the reference genomes.

### [Outputs](./outputs)


## Running the code in this repository

This repository uses snakemake to run the pipeline and conda to manage software environments and installations.
You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html).
We ran this pipeline on an AWS EC2 instance with an Ubuntu image (ubuntu/images/hvm-ssd/ubuntu-jammy-22.04-amd64-server-20230208).
We used the following commands to set up conda and mamba.

```
curl -JLO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh # download the miniconda installation script
bash Miniconda3-latest-Linux-x86_64.sh # run the miniconda installation script. Accept the license and follow the defaults.
source ~/.bashrc # source the .bashrc for miniconda to be available in the environment
# configure miniconda channel order
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict # make channel priority strict so snakemake doesn't yell at you
conda install mamba # install mamba for faster software installation.
```

After installing conda and mamba, you can clone the repository and use the following command to create the conda environment.
```
conda env create -n brain -f environment.yml
conda activate brain
```

You can then use the following command to run one of the snakefiles:

```
snakemake -s 00_subtractive_approach.snakefile -j 1 --use-conda --rerun-incomplete -k -n
```

where:
* `-s` specifies which snakefile to run 
* `-j` specifies the number of threads to run with
* `--use-conda` uses conda to manage software environments
* `--rerun-incomplete` re-runs incomplete files
* `-k` tells the pipeline to continue with independent steps when one step fails
* `-n` signifies to run a dry run first
