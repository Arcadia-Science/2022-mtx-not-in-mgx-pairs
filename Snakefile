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
        expand("outputs/sourmash_sketch_subtract_describe/{mtx_minus_mgx}-k{ksize}.csv", mtx_minus_mgx = MTX_MINUS_MGX, ksize = KSIZES),
        expand("outputs/sourmash_sketch_describe/{run_accession}.csv", run_accession = RUN_ACCESSIONS)
        un = "outputs/sourmash_sketch_subtract_gather_unassigned/{mtx_minus_mgx}-vs-genbank-2022.03-k{ksize}-unassigned.sig"

####################################################
## Sketch metagenomes and metatranscriptomes 
####################################################

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

########################################################
## Subtract the metagenomes from the metatranscriptomes
########################################################

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

#########################################################
## Determine the number of k-mers in original signatures
## and in subtracted signatures
#########################################################

rule sourmash_sig_describe_sketches:
    """
    Use the sourmash CLI to report detailed information about all sketches, including number of hashes.
    Output the information as a csv file. 
    """
    input: "outputs/sourmash_sketch/{run_accession}.sig"
    output: "outputs/sourmash_sketch_describe/{run_accession}.csv"
    conda: 'envs/sourmash.yml'
    shell:'''
    sourmash sig describe --csv {output} {input}
    '''

rule sourmash_sig_describe_subtracted_sketches:
    """
    Use the sourmash CLI to report detailed information about all sketches, including number of hashes.
    Output the information as a csv file. 
    """
    input: "outputs/sourmash_sketch_subtract/{mtx_minus_mgx}-k{ksize}.sig"
    output: "outputs/sourmash_sketch_subtract_describe/{mtx_minus_mgx}-k{ksize}.csv"
    conda: 'envs/sourmash.yml'
    shell:'''
    sourmash sig describe --csv {output} {input}
    '''

#######################################################
## Profile the taxonomy of the sequences leftover in a
## metatranscriptome after subtracting a metagenome
#######################################################

rule sourmash_gather:
    input:
        sig = "outputs/sourmash_sketch_subtract/{mtx_minus_mgx}-k{ksize}.sig",
        databases=expand("inputs/sourmash_databases/genbank-2022.03-{lineage}-k{{ksize}}.zip", lineage = LINEAGES)
    output: 
        csv="outputs/sourmash_sketch_subtract_gather/{mtx_minus_mgx}-vs-genbank-2022.03-k{ksize}.csv",
        un = "outputs/sourmash_sketch_subtract_gather_unassigned/{mtx_minus_mgx}-vs-genbank-2022.03-k{ksize}-unassigned.sig"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash gather -k {wildcards.ksize} --scaled 1000 --output-unassigned {output.un} --threshold-bp 0 -o {output.csv} {input.sig} {input.databases}
    '''
   
rule gunzip_lineage_csvs:
    input: "inputs/sourmash_databases/genbank-2022.03-{lineage}.lineages.csv.gz"
    output: "inputs/sourmash_databases/genbank-2022.03-{lineage}.lineages.csv"
    shell:'''
    gunzip -c {input} > {output}
    '''

rule sourmash_taxonomy_prepare:
    input: expand("inputs/sourmash_databases/genbank-2022.03-{lineage}.lineages.csv", lineage = LINEAGES),
    output: "inputs/sourmash_databases/genbank-2022.03-prepared-lineages.sqldb"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash tax prepare --taxonomy-csv {input} -o {output}
    '''

rule sourmash_taxonomy_annotate:
   input:
       lin_prepared="inputs/sourmash_databases/genbank-2022.03-prepared-lineages.sqldb",
       gather="outputs/sourmash_sketch_subtract_gather/{mtx_minus_mgx}-vs-genbank-2022.03-k{ksize}.csv"
   output: "outputs/sourmash_sketch_subtract_taxonomy/{mtx_minus_mgx}-vs-genbank-2022.03-k{ksize}.with-lineages.csv"
   params: outdir = "outputs/sourmash_sketch_subtract_taxonomy/"
   conda: "envs/sourmash.yml"
   shell:'''
   sourmash tax annotate -g {input.gather} -t {input.lin_prepared} -o {params.outdir}
   '''

##########################################################
## Download sourmash databases & taxonomy files
##########################################################

rule download_genbank_bacteria_zip_k21:
    output: "inputs/sourmash_databases/genbank-2022.03-bacteria-k21.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/6qxfp/download
    '''

rule download_genbank_bacteria_zip_k31:
    output: "inputs/sourmash_databases/genbank-2022.03-bacteria-k31.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/9ue5g/download
    '''

rule download_genbank_bacteria_zip_k51:
    output: "inputs/sourmash_databases/genbank-2022.03-bacteria-k51.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/5gvbw/download
    '''

rule download_genbank_bacteria_lineage:
    output: "inputs/sourmash_databases/genbank-2022.03-bacteria.lineages.csv.gz"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/4agsp/download
    '''

rule download_genbank_fungi_zip_k21:
    output: "inputs/sourmash_databases/genbank-2022.03-fungi-k21.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/fy82q/download
    '''
rule download_genbank_fungi_zip_k31:
    output: "inputs/sourmash_databases/genbank-2022.03-fungi-k31.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/4pdbj/download
    '''
rule download_genbank_fungi_zip_k51:
    output: "inputs/sourmash_databases/genbank-2022.03-fungi-k51.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/b9a86/download
    '''

rule download_genbank_fungi_lineage:
    output: "inputs/sourmash_databases/genbank-2022.03-fungi.lineages.csv.gz"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/s4b85/download
    '''

rule download_genbank_archaea_zip_k21:
    output: "inputs/sourmash_databases/genbank-2022.03-archaea-k21.zip"
    conda: "envs/wget.yml"
    shell:''' 
    wget -O {output} https://osf.io/g94n5/download
    '''

rule download_genbank_archaea_zip_k31:
    output: "inputs/sourmash_databases/genbank-2022.03-archaea-k31.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/hfybv/download
    '''

rule download_genbank_archaea_zip_k51:
    output: "inputs/sourmash_databases/genbank-2022.03-archaea-k51.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/dehrc/download
    '''

rule download_genbank_archaea_lineage:
    output: "inputs/sourmash_databases/genbank-2022.03-archaea.lineages.csv.gz"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/kcbpn/download
    '''

rule download_genbank_viral_zip_k21:
    output: "inputs/sourmash_databases/genbank-2022.03-viral-k21.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/updvc/download
    '''

rule download_genbank_viral_zip_k31:
    output: "inputs/sourmash_databases/genbank-2022.03-viral-k31.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/egkt2/download
    '''

rule download_genbank_viral_zip_k51:
    output: "inputs/sourmash_databases/genbank-2022.03-viral-k51.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/z8scg/download
    '''

rule download_genbank_viral_lineage:
    output: "inputs/sourmash_databases/genbank-2022.03-viral.lineages.csv.gz"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/j4tsu/download
    '''

rule download_genbank_protozoa_zip_k21:
    output: "inputs/sourmash_databases/genbank-2022.03-protozoa-k21.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/m23r6/download 
    '''

rule download_genbank_protozoa_zip_k31:
    output: "inputs/sourmash_databases/genbank-2022.03-protozoa-k31.zip"
    conda: "envs/wget.yml"
    shell:''' 
    wget -O {output} https://osf.io/zm5vg/download
    '''

rule download_genbank_protozoa_zip_k51:
    output: "inputs/sourmash_databases/genbank-2022.03-protozoa-k51.zip"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/32y98/download
    '''

rule download_genbank_protist_lineage:
    output: "inputs/sourmash_databases/genbank-2022.03-protozoa.lineages.csv.gz"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://osf.io/2x8u4/download
    '''

rule download_human_signature:
    output: "inputs/sourmash_databases/GCF_000001405.39_GRCh38.p13_genomic.sig"
    shell:'''
    wget -O {output} https://osf.io/fxup3/download
    '''

rule download_euk_rna_db:

