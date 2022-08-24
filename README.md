# Identifying sequences that are in a metatranscriptome but not in a metagenome

This repository curates a set of paired metagenomes and metatranscriptomes and provides a pipeline to rapidly identify the fraction of sequences in a metatranscriptome that are not in a metagenome.
The pipeline is shaped around a metadata file, `inputs/metadata-paired-mgx-mtx.tsv`, that contains a sample name (`sample_name`), metagenome SRA run accession (`mgx_run_accession`; `SRR*`, `ERR*`, `DRR*`), metatranscriptome SRA run accession (`mtx_run_accession`), and a sample type (`sample_type`).
Using the run accessions, it downloads the sequencing data from the SRA and generates a [FracMinHash sketch](https://www.biorxiv.org/content/10.1101/2022.01.11.475838v2.abstract) of each run.
Then, it uses the paired information encoded in the metadata table to subtract the metagenome sketch from the metatranscriptome sketch.
This produces an estimate of the fraction metatranscriptome sequences not found in the paired metagenome.
These estimates are also clustered by `sample_type` to generate biome-specific estimates.
The pipeline also analyzes the fraction of metatranscriptome-specific sequences that are shared between samples to discover what fraction of sequences we are systematically missing within and across biomes.

Some metagenome and metatranscripome pairs are true pairs that were extracted from the same sample while others are from separate samples taken from the same location at the same time.

## Getting started with this repository

This repository uses snakemake to run the pipeline and conda to manage software environments and installations.
You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html).
After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the pipeline run environment.

```
mamba env create -n mtx_mgx --file environment.yml
conda activate mtx_mgx
```

To start the pipeline, run:
```
snakemake --use-conda -j 2
```

### Running this repository on AWS

This repository was executed on an AWS EC2 instance (Ubuntu 22.04 LTS ami-085284d24fe829cd0, t2.large, 100 GiB root storage).
The instance was configured using the following commands:

```
curl -JLO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh # download the miniconda installation script
bash Miniconda3-latest-Linux-x86_64.sh # run the miniconda installation script. Accept the license and follow the defaults.
source ~/.bashrc # source the .bashrc for miniconda to be available in the environment

# configure miniconda channel order
conda config --add channels defaults 
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install mamba # install mamba for faster software installation.
```

Once miniconda is configured, clone the repository and `cd` into it, then follow the instructions in the above section.
```
git clone https://github.com/Arcadia-Science/2022-mtx-not-in-mgx-pairs.git
cd 2022-mtx-not-in-mgx-pairs
```

Note that there is a [known issue](https://github.com/ncbi/sra-tools/issues/645) with the conda package `sra-tools=2.11.0` that will prevent this workflow from running on Mac operating systems. 

## Next steps

TBD
