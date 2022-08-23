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

rule all:
    input: expand("outputs/sourmash_sketch_subtract/{sample}_k{ksize}.sig", sample = SAMPLES, ksize = KSIZES)

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
    Use the sourmash python API to subtract metagenome fracminhash sketches from metatranscriptome fracminhash sketches.
    This for loop re-implements a simplified veresion of the the CLI of sourmash sig subtract.
    https://github.com/sourmash-bio/sourmash/blob/latest/src/sourmash/sig/__main__.py
    Doing this in a python script is the easiest way to map pairs of samples to each other within the snakemake workflow.
    
    The run directive doesn't allow a run to have a conda environment; 
    therefore the run environment for this code is specified in environment.yml
    """
    input:
        metadata = 'inputs/metadata-paired-mgx-mtx.tsv',
        sigs = expand('outputs/sourmash_sketch/{run_accession}.sig', run_accession = RUN_ACCESSIONS) 
    output:
        sigs = expand("outputs/sourmash_sketch_subtract/{sample}_k{{ksize}}.sig", sample = SAMPLES)
    run:
        ksize = wildcards.ksize
        # read in metadata dataframe to derive sample pairs
        metadata = pd.read_csv(input.metadata, sep = "\t") # read in metadata as pandas df
        metadata = metadata.reset_index()  # make sure indexes pair with number of rows
        
        for index, row in metadata.iterrows():
            # grab the run accessions for a given paired metagenome and metatranscriptome
            mtx_run_accession = row['mtx_run_accession']
            mgx_run_accession = row['mgx_run_accession']
            # using the run assession, create paths of the signatures
            mtx_sigfile = os.path.join('outputs/sourmash_sketch', mtx_run_accession + '.sig')
            mgx_sigfile = os.path.join('outputs/sourmash_sketch', mgx_run_accession + '.sig')

            # read in mtx signature and grab both the minhash object and the hashes in that object
            mtx_sigobj = sourmash.load_one_signature(mtx_sigfile, ksize=ksize, select_moltype='DNA')
            mtx_mh = from_sigobj.minhash
            # turn the mtx hashes into a set
            mtx_subtract_mins = set(mtx_mh.hashes)

            # read in mgx signature
            mgx_sigobj = sourmash.load_one_signature(mgx_sig_path, ksize=ksize, select_moltype='DNA')

            # do the subtraction
            mtx_subtract_mins -= set(mgx_sigobj.minhash.hashes)

            # make a new signature that only contains the mtx hashes that were not in the mgx
            mtx_subtract_mh = mtx_sigobj.minhash.copy_and_clear().flatten()
            mtx_subtract_mh.add_many(mtx_subtract_mins)

            # re-establish hash abundances from mtx sig
            # re-read in mtx sig just in case abundances were removed in place
            abund_sig = sourmash.load_one_signature(mtx_sig_path, ksize=ksize, select_moltype='DNA')
            mtx_subtract_mh = mtx_subtract_mh.inflate(abund_sig.minhash)
  
            # create a new sig object that can be written to file
            mtx_subtract_sigobj = sourmash.SourmashSignature(mtx_subtract_mh)

            # create a name that reflects the signature contents
            # name = mtx_run_accession + "-minus-" + mgx_run_accession
            name = row['sample_name']
            mtx_subtract_sigobj.name = name
            
            # create output sig file name
            sig_filename = os.path.join('outputs/sourmash_sketch_subtract/' + name + "_k" + ksize + ".sig")
            # write sig to file
            with sourmash.sourmash_args.FileOutput(sig_filename, 'wt') as fp:
                sourmash.save_signatures([mtx_subtract_sigobj], fp=fp)



