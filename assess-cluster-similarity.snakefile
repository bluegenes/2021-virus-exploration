"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
"""
import os
import re
import pandas as pd

configfile: "conf-cluster-similarity.yml"
basename = "pigeon1.0-vc"
out_dir = config["output_dir"]
logs_dir = os.path.join(out_dir, "logs")
expt = config.get("experiment", "")
if expt:
    expt = "_" + expt
compare_dir = os.path.join(out_dir, "compare" + expt)

# load acc::fasta filenames info
fasta_fileinfo = pd.read_csv(config["fasta_info"]).set_index("accession")
# load compareInfo and turn ;-separated list into python list
compareInfo = pd.read_csv(config["comparison_info"]).set_index("cluster")
compareInfo["cluster_members"] = compareInfo["cluster_members"].str.split(";")
# get correct alphabets, ksizes for comparisons
alphabet_info = config["alphabet_info"]
genomic_alphaksizes, protein_alphaksizes = [],[]
for alpha in alphabet_info:
    ak = expand("{alpha}-k{k}", alpha=alpha, k=alphabet_info[alpha]["ksizes"])
    if alpha != "nucleotide":
        protein_alphaksizes+=ak
    else:
        #I didn't do translated sigs with pigeon -- nucl should only have nucl sigs
        genomic_alphaksizes+=ak

rule all:
    input: 
        # fastani
        os.path.join(compare_dir, "fastani", f"{basename}.fastani.csv.gz"),
        # compareM
        expand(os.path.join(compare_dir, "compareM", "{basename}.{input_type}.compareM.csv.gz"), input_type=["genomic", "protein"], basename=basename),
        # sourmash
        expand(os.path.join(compare_dir, "cluster-compare", "{basename}.{input_type}.clustercompare.csv.gz"), basename=basename, input_type=["genomic", "protein"]),


#####################
# fastANI comparisons
######################

def get_fastani_comparison_genome_files(w):
    compare_accs = compareInfo.at[w.cluster, "cluster_members"]
    anchor = compareInfo.at[w.cluster, "cluster_anchor"]
    genome_paths = []
    for acc in compare_accs:
        if acc != anchor:
            genome_paths += [os.path.abspath(fasta_fileinfo.at[acc, "genome"])]
    return genome_paths 


localrules: write_genomic_fastani_fastalist
rule write_genomic_fastani_fastalist:
    input: ancient(get_fastani_comparison_genome_files)
    output: os.path.join(compare_dir, "fastani", "{cluster}", "{cluster}.genomic.fastalist")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")

def get_fastani_comparison_info(w):
    compare_accs = compareInfo.at[w.cluster, "cluster_members"]
    anchor = compareInfo.at[w.cluster, "cluster_anchor"]
    genome_paths = []
    for acc in compare_accs:
        if acc != anchor:
            genome_paths += [os.path.abspath(fasta_fileinfo.at[acc, "genome"])]
    anchor_g = os.path.abspath(fasta_fileinfo.at[anchor, "genome"])
    c_filelist = os.path.join(compare_dir, "fastani", f"{w.cluster}", f"{w.cluster}.genomic.fastalist")
    return {"anchor_genome" : anchor_g, "comparison_filelist": c_filelist}

rule compare_via_fastANI:
    input:
        unpack(get_fastani_comparison_info)
    output: 
        os.path.join(compare_dir, "fastani", "{cluster}.fastani.tsv"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        runtime=1200,
    log: os.path.join(logs_dir, "fastani", "{cluster}.fastani.log")
    benchmark: os.path.join(logs_dir, "fastani", "{cluster}.fastani.benchmark")
    conda: "/home/ntpierce/2020-distance-compare/envs/fastani-env.yml"
    shell:
        """
        fastANI -q {input.anchor_genome:q} --rl {input.comparison_filelist:q} -o {output} > {log} 2>&1
        """

## aggreagate fastani results
def get_all_fastani(w):
    fastani_files =  expand(os.path.join(compare_dir, "fastani", "{cluster}.fastani.tsv"), cluster=compareInfo.index)
    return fastani_files

localrules: write_fastani_result_csv
rule write_fastani_result_csv:
    input: get_all_fastani 
    output: os.path.join(compare_dir, "fastani", "{basename}.fastani.filecsv"),
    run:
        with open(str(output), "w") as out:
            for inF in input:
                comparison_name = os.path.basename(str(inF)).rsplit(".fastani.tsv")[0]
                out.write(f"{comparison_name},{str(inF)}\n")

localrules: aggregate_fastani_results
rule aggregate_fastani_results:
    input:
        fastani=os.path.join(compare_dir, "fastani", "{basename}.fastani.filecsv"),
        comparison_info=config["comparison_info"],
    output: os.path.join(compare_dir, "fastani", "{basename}.fastani.csv.gz"),
    log: os.path.join(logs_dir, "fastani", "{basename}.fastani.aggregate.log")
    benchmark: os.path.join(logs_dir, "fastani", "{basename}.fastani.aggregate.benchmark")
    shell:
        """
        python aggregate-cluster-fastani.py --fastani-filecsv {input.fastani} \
                                            --comparison-info {input.comparison_info} \
                                            --output-csv {output} > {log} 2>&1
        """


#####################
# compareM comparisons
######################

def get_compareM_protein_fastas(w):
    compare_accs = compareInfo.at[w.cluster, "cluster_members"] # for pigeon, includes anchor genome
    protein_fastas = []
    for acc in compare_accs:
        protein_fastas += [os.path.abspath(fasta_fileinfo.at[acc, "protein"])]
    return protein_fastas

rule write_protein_compareM_fastalist:
    input: ancient(get_compareM_protein_fastas)
    output: os.path.join(compare_dir, "compareM", "{cluster}", "protein", "{cluster}.fastalist")
    params:
        outdir = lambda w: os.path.join(compare_dir, "compareM", f"{w.cluster}", "protein")
    run:
        with open(str(output), "w") as out:
            for fn in input:
                out.write(f"{str(fn)}\n")

rule protein_AAI_via_compareM:
    input:
        os.path.join(compare_dir, "compareM", "{cluster}", "protein", "{cluster}.fastalist")
    output:
        os.path.join(compare_dir, "compareM", "{cluster}/protein/aai/aai_summary.tsv"),
    params:
        proteins_cmd = "--proteins",
        file_ext = ".faa.gz",
        outdir = lambda w: os.path.join(compare_dir, "compareM", f"{w.cluster}/protein"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        runtime=60,
    log: os.path.join(logs_dir, "compareM/protein", "{cluster}.compareM.log")
    benchmark: os.path.join(logs_dir, "compareM/protein", "{cluster}.compareM.benchmark")
    conda: "/home/ntpierce/2020-distance-compare/envs/compareM-env.yml"
    shell:
        """
        comparem aai_wf --cpus {threads} {params.proteins_cmd} --file_ext {params.file_ext:q}  --sensitive {input} {params.outdir} > {log} 2>&1
        """

## nucleotide compareM ##
def get_compareM_genome_fastas(w):
    compare_accs = compareInfo.at[w.cluster, "cluster_members"] # includes anchor genome
    genome_paths = []
    for acc in compare_accs:
        genome_paths += [os.path.abspath(fasta_fileinfo.at[acc, "genome"])]
    return genome_paths
    
# note, prodigal cant use gzipped nucl files. not an issue here; see pseudomonas_compare.v2.snakefile for hacky workaround
rule write_genomic_compareM_fastalist:
    input: ancient(get_compareM_genome_fastas)
    output: os.path.join(compare_dir, "compareM", "{cluster}", "genomic", "{cluster}.fastalist")
    params:
        outdir = lambda w: os.path.join(compare_dir, "compareM", f"{w.cluster}", "genomic")
    group: "nuclcompareM"
    run:
        with open(str(output), "w") as out:
            for fn in input:
                out.write(f"{fn}\n")

rule nucl_AAI_via_compareM:
    input:
        os.path.join(compare_dir, "compareM", "{cluster}", "genomic", "{cluster}.fastalist")
    output:
        os.path.join(compare_dir, "compareM", "{cluster}/genomic/aai/aai_summary.tsv"),
    params:
        proteins_cmd = "",
        file_ext = ".fa",
        outdir = lambda w: os.path.join(compare_dir, "compareM", f"{w.cluster}/genomic"),
    group: "nuclcompareM"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        runtime=60,
    log: os.path.join(logs_dir, "compareM/genomic", "{cluster}.compareM.log")
    benchmark: os.path.join(logs_dir, "compareM/genomic", "{cluster}.compareM.benchmark")
    conda: "/home/ntpierce/2020-distance-compare/envs/compareM-env.yml"
    shell:
        """
        comparem aai_wf --cpus {threads} {params.proteins_cmd} --file_ext {params.file_ext:q}  --sensitive {input} {params.outdir} > {log} 2>&1
        """

## aggreagate compareM results
def get_all_compareM(w):
    compareM_files = expand(os.path.join(compare_dir, "compareM", "{cluster}/{inp}/aai/aai_summary.tsv"), cluster=compareInfo.index, inp=w.input_type)
    return compareM_files

rule compile_compareM_resultfiles:
    input: get_all_compareM
    output: 
        os.path.join(compare_dir, "compareM", "{basename}.{input_type}.compareM.filecsv",),
    run:
        with open(str(output), "w") as out:
            for inF in input:
                comparison_name = os.path.basename(str(inF).rsplit(f"/{wildcards.input_type}")[0])
                out.write(f"{comparison_name},{str(inF)}\n")

localrules: aggregate_compareM_results
rule aggregate_compareM_results:
    input:
        compareM=os.path.join(compare_dir, "compareM", "{basename}.{input_type}.compareM.filecsv"),
        comparison_info=config["comparison_info"],
    output: os.path.join(compare_dir, "compareM", "{basename}.{input_type}.compareM.csv.gz"),
    log: os.path.join(logs_dir, "compareM", "{basename}.{input_type}.compareM.aggregate.log")
    benchmark: os.path.join(logs_dir, "compareM", "{basename}.{input_type}.compareM.aggregate.benchmark")
    shell:
        """
        python aggregate-cluster-compareM.py --comparem-tsv-filecsv {input.compareM} \
                                             --comparison-info {input.comparison_info} \
                                             --output-csv {output} > {log} 2>&1
        """


#####################
# sourmash comparisons
######################
## compare to anchor sigs ##
alpha_to_moltype = {"nucleotide": "DNA", "protein": "protein", "dayhoff": "dayhoff", "hp": "hp"}
rule cluster_compare_genomic:
    input:
        comparison_csv=config["comparison_info"],
        sigfile=config["genome_sigfile"]
    output:
        csv=os.path.join(compare_dir, "cluster-compare/genomic", "{basename}.{alphabet}-k{ksize}.clustercompare.csv.gz"),
    params:
        sigdir = os.path.join(out_dir, "signatures"),
        moltype = lambda w: alpha_to_moltype[w.alphabet],
        #sigext = lambda w: f".{w.input_type}.sig"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=1200,
    log: os.path.join(logs_dir, "cluster-compare/genomic", "{basename}.{alphabet}-k{ksize}.clustercompare.log")
    benchmark: os.path.join(logs_dir, "cluster-compare/genomic", "{basename}.{alphabet}-k{ksize}.clustercompare.benchmark")
    conda: "/home/ntpierce/2020-distance-compare/envs/pathcompare.yml"
    shell:
        """
        python cluster-compare-singleton.py --comparison-csv {input.comparison_csv} \
        --alphabet {params.moltype} --ksize {wildcards.ksize} --sigdir {params.sigdir} \
        --sigfiles {input.sigfile} --output-csv {output.csv} > {log} 2>&1
        """

rule cluster_compare_protein:
    input:
        comparison_csv=config["comparison_info"],
        siglist=config["protein_siglist"]
    output:
        csv=os.path.join(compare_dir, "cluster-compare/protein", "{basename}.{alphabet}-k{ksize}.clustercompare.csv.gz"),
    params:
        sigdir = config.get("protein_sigdir", "./"),
        sigext= config.get("protein_sigext", ""),
        sigprefix = config.get("protein_sigprefix", ""),
        moltype = lambda w: alpha_to_moltype[w.alphabet],
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=1200,
    log: os.path.join(logs_dir, "cluster-compare/protein", "{basename}.{alphabet}-k{ksize}.clustercompare.log")
    benchmark: os.path.join(logs_dir, "cluster-compare/protein", "{basename}.{alphabet}-k{ksize}.clustercompare.benchmark")
    conda: "/home/ntpierce/2020-distance-compare/envs/pathcompare.yml"
    shell:
        """
        python cluster-compare.py --comparison-csv {input.comparison_csv} \
        --alphabet {params.moltype} --ksize {wildcards.ksize} --sigdir {params.sigdir} \
        --sig-extension {params.sigext:q} --sig-prefix {params.sigprefix:q} \
        --siglist {input.siglist} --output-csv {output.csv} > {log} 2>&1
        """

localrules: aggregate_genomic_clustercompare
rule aggregate_genomic_clustercompare:
    input:
        expand(os.path.join(compare_dir, "cluster-compare/genomic", "{basename}.{alphak}.clustercompare.csv.gz"), basename=basename, alphak=genomic_alphaksizes)
    output:
        os.path.join(compare_dir, "cluster-compare", "{basename}.genomic.clustercompare.csv.gz")
    run:
        # aggreate all csv.gzs --> single csv
        aggDF = pd.concat([pd.read_csv(str(csv), sep=",") for csv in input])
        aggDF["alpha-ksize"] = aggDF["alphabet"] + "-" + aggDF["ksize"].astype(str)
        aggDF.to_csv(str(output), index=False)

localrules: aggregate_protein_clustercompare
rule aggregate_protein_clustercompare:
    input:
        expand(os.path.join(compare_dir, "cluster-compare/protein", "{basename}.{alphak}.clustercompare.csv.gz"), basename=basename, alphak=protein_alphaksizes)
    output:
        os.path.join(compare_dir, "cluster-compare", "{basename}.protein.clustercompare.csv.gz")
    run:
        # aggreate all csv.gzs --> single csv
        aggDF = pd.concat([pd.read_csv(str(csv), sep=",") for csv in input])
        aggDF["alpha-ksize"] = aggDF["alphabet"] + "-" + aggDF["ksize"].astype(str)
        aggDF.to_csv(str(output), index=False)

