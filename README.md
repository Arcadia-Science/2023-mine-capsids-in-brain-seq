# Mining viral sequences from human brain sequencing data via capsid searches

This repository is associated with the pub, ["A capsid-based search recovers viral sequences from human brain sequencing data"](ADD LINK DOI WHEN AVAILABLE).

## Explanation of the files in this repository

Below we provide an explanation of the files in this repository.


### Analysis code

#### Snakefiles

This repository 
[inputs](./inputs)
[dags](./dags)
[envs](./envs)
[scripts](./scripts)

#### [Notebooks](./notebooks)

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
