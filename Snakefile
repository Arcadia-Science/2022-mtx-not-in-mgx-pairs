import pandas as pd
import sourmash
import os

# create wildcard variables for workflow
metadata = pd.read_csv("inputs/metadata-paired-mgx-mtx.tsv", sep = "\t") # read in metadata as a pandas dataframe
MTX = metadata['mtx_run_accession'].unique().tolist()                    # make run accession in mtx col into list
MGX = metadata['mgx_run_accession'].unique().tolist()                    # make run accession in mgx col into list
RUN_ACCESSIONS = MGX + MTX                                               # combine mtx and mgx into one list with all run accessions
SAMPLES = metadata['sample_name'].unique().tolist()                      # make a list of sample names
KSIZES = [21, 31, 51]                                                    # create a list of k-mer sizes for the workflow
MTX_MINUS_MGX = [x + '-minus-' + y for x, y in zip(MTX, MGX)]            # create list that binds mtx to mgx

rule all:
    input: 
        #expand("outputs/sourmash_sketch_subtract/{mtx_minus_mgx}-k{ksize}.sig", mtx_minus_mgx = MTX_MINUS_MGX, ksize = KSIZES)
        "outputs/sourmash_sig_describe/sourmash_sig_describe.csv"

rule sourmash_sketch:
    """
    Use fastq-dump to download sequencing data for a given accession.
    Pipe the sequencing data directly to a sketch without writing to disk.
    Using fastq-dump like this dumps both single end and paired end reads into std out, which can then be piped into sourmash sketch
    """
    output: "outputs/sourmash_sketch/{run_accession}.sig"
    conda: 'envs/sourmash.yml'
    shell:'''
    fastq-dump --disable-multithreading --fasta 0 --skip-technical --readids --read-filter pass --dumpbase --split-spot --clip -Z {wildcards.run_accession} | 
        sourmash sketch dna -p k=21,k=31,k=51,scaled=200,abund --name {wildcards.run_accession} -o {output} -
    '''

rule calculate_mtx_not_in_mgx:
    """
    Use the sourmash CLI to subtract a metagenome sketch from its paired metatranscriptome, retaining metatranscriptome abundances.
    Snakemake will auto-parse the wildcard mtx_minus_mgx on the "-minus-" string to back-propagate the correct wildcards to m*x_run_accession.
    Providing the accessions in this way limits the number of ways they can combine;
    If MTX solved for mtx_run_accession and MGX solved for mgx_run_accession, then snakemake would execute this rule for all combinations of MTX and MGX.
    Instead, by binding pairs together in MTX_MINUS_MGX, only pairs of metatranscriptomes and metagenomes are subtracted.
    """
    input:
        mtx_sig = 'outputs/sourmash_sketch/{mtx_run_accession}.sig', 
        mgx_sig = 'outputs/sourmash_sketch/{mgx_run_accession}.sig', 
    output: "outputs/sourmash_sketch_subtract/{mtx_run_accession}-minus-{mgx_run_accession}-k{ksize}.sig"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sig subtract -k {wildcards.ksize} -o {output} -A {input.mtx_sig} {input.mtx_sig} {input.mgx_sig}
    '''

rule sourmash_sig_describe_sketches:
    """
    Use the sourmash CLI to report detailed information about all sketches, including number of hashes.
    Output the information as a csv file. 
    """
    input: 
        expand("outputs/sourmash_sketch/{run_accession}.sig", run_accession = RUN_ACCESSIONS),
        expand("outputs/sourmash_sketch_subtract/{mtx_minus_mgx}-k{ksize}.sig", mtx_minus_mgx = MTX_MINUS_MGX, ksize = KSIZES)
    output: "outputs/sourmash_sig_describe/sourmash_sig_describe.csv"
    conda: 'envs/sourmash.yml'
    shell:'''
    sourmash sig describe --csv {output} {input}
    '''

