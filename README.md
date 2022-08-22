# Identifying sequences that are in a metatranscriptome but not in a metagenome

This repository curates a set of paired metagenomes and metatranscriptomes and provides a pipeline to rapidly identify the fraction of sequences in a metatranscriptome that are not in a metagenome.
The pipeline is shaped around a metadata file, 'inputs/metadata-paired-mgx-mtx.tsv', that contains a sample name (`sample_name`), metagenome SRA run accession (`mgx_run_accession`; `SRR*`, `ERR*`, `DRR*`), metatranscriptome SRA run accession (`mtxx_run_accession`), and a sample type (`sample_type`).
Using the run accessions, it downloads the sequencing data from the SRA and generates a [FracMinHash sketch](https://www.biorxiv.org/content/10.1101/2022.01.11.475838v2.abstract) of each run.
Then, it uses the paired information encoded in the metadata table to subtract the metagenome sketch from the metatranscriptome sketch.
This produces an estimate of the fraction metatranscriptome sequences not found in the paired metagenome.
These estimates are also clustered by `sample_type` to generate biome-specific estimates.
The pipeline also analyzes the fraction of metatranscriptome-specific sequences that are shared between samples to discover what fraction of sequences we are systematically missing within and across biomes.

Some metagenome and metatranscripome pairs are true pairs that were extracted from the same sample while others are from separate samples taken from the same location at the same time.

## Getting started with this repository

TBD

## Next steps

TBD
