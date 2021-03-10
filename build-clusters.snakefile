#####
# Use sourmash uniqify + matches to cluster a collection of sigs 
# snakemake -s cluster-sigs.snakefile --profile default -n
# Author: N. Tessa Pierce Ward, ntpierce@gmail.com
#####

import os
import pandas as pd

configfile: "conf-cluster.yml"
out_dir = config["output_dir"]
logs_dir = os.path.join(out_dir,"logs")

#basename = config.get("basename", "pigeon1.0")
prefix = config.get("prefix", "pigeon1.0")

alphabet_info = config["alphabet_info"]
alphakmc_params = []
for alpha, info in alphabet_info.items():
    alphakmc_params += expand("{alphabet}-k{ksize}.mc{maxcontain}", alphabet=alpha, ksize=info["ksizes"], maxcontain=info["maxcontain_threshold"])

#{alphabet}-k{ksize}.mc{maxcontain}



# to do: edits for this workflow
class Checkpoint_MakePattern:
    def __init__(self, pattern):
        self.pattern = pattern

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'check_csv'; this will trigger an
        # exception until that rule has been run.
        #checkpoints.check_csv.get(**w)
        checkpoints.find_founders.get(**w)

        global clusterMembers
        clusterMembers = {}
        
        with open(f"{out_dir}/{w.prefix}.{w.alphabet}-k{w.ksize}.mc{w.maxcontain}.members.siglist.csv") as fp:
            for line in fp:
                name, sigfile = line.rstrip().split(',')
                clusterMembers[name] = sigfile
        
        names = clusterMembers.keys()
        
        pattern = expand(self.pattern, name=names, **w)
        return pattern


# snakemake rules
rule all:
    input:
        expand(os.path.join(out_dir, "{prefix}.{akm}.cluster-info.csv"), prefix=prefix, akm=alphakmc_params)


checkpoint find_founders:
    message:
        """
        Use sourmash-uniqify logic to nucleate and populate clusters based on maximum containment. 
        """
    input: 
        config["siglist"]
    output:
        founders = os.path.join(out_dir, "{prefix}.{alphabet}-k{ksize}.mc{maxcontain}.founders.siglist.csv"),
        founders_txt = os.path.join(out_dir, "{prefix}.{alphabet}-k{ksize}.mc{maxcontain}.founders.siglist.txt"),
        members = os.path.join(out_dir, "{prefix}.{alphabet}-k{ksize}.mc{maxcontain}.members.siglist.csv")
    conda:
        "envs/sourmash4.0.yml"
    log: os.path.join(logs_dir, "find_founders", "{prefix}.{alphabet}-k{ksize}.mc{maxcontain}.find-founders.log" )
    benchmark: os.path.join(logs_dir, "find_founders", "{prefix}.{alphabet}-k{ksize}.mc{maxcontain}.find-founders.benchmark" )
    shell:
        """
        python find-founders.py --siglist {input} --threshold {wildcards.maxcontain} \
                                --moltype {wildcards.alphabet} \
                                --ksize {wildcards.ksize} > {log} 2&>1
        """


## probably need to use a checkpoint here, like protein-pigeon.
# 1. run find_founders
# 2. checkpoint -- check that find_founders has produced relevant text file
# 3. as part of checking checkpoint, read in `members.siglist.txt` --> global python dictionary
# 4. use the members in cluster_sig_to_founders; write individual files. make these temporary?
# 5. aggregate the individual files / temp files --> summary csv of cluster info
#  -->  cluster_name/num, founder_name, founder_sigfile, membersig_filename, member sig name, max containment, num common hashes


#localrules: check_csv

#checkpoint check_csv:
#    input:
#        founders = os.path.join(out_dir, "{prefix}.{alphabet}-k{ksize}.mc{maxcontain}.founders.siglist.csv"),
#        members = os.path.join(out_dir, "{prefix}.{alphabet}-k{ksize}.mc{maxcontain}.members.siglist.csv")
#    output: touch(os.path.join(out_dir, ".{prefix}.{alphabet}-k{ksize}.mc{maxcontain}.make_spreadsheet.touch"))
   #output: touch(f"{out_dir}/.make_spreadsheet.touch")

rule sourmash_index:
    input:
        siglist = os.path.join(out_dir, "{prefix}.{alphabet}-k{ksize}.mc{maxcontain}.founders.siglist.txt"),
    output:
        db = os.path.join(out_dir, "{prefix}.{alphabet}-k{ksize}.mc{maxcontain}.founders.sbt.zip"),
    params:
        alpha_cmd= lambda w: config["alphabet_info"][w.alphabet]["alpha_cmd"],
    log: os.path.join(logs_dir, "index/{prefix}.{alphabet}-k{ksize}.mc{maxcontain}.founders.sbt.log"),
    benchmark: os.path.join(logs_dir, "index/{prefix}.{alphabet}-k{ksize}.mc{maxcontain}.founders.sbt.benchmark"),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=600000,
    conda: "envs/sourmash4.0.yml"
    shell:
        """
        sourmash index -k {wildcards.ksize} {params.alpha_cmd} --from-file {input.siglist} {output.db} 2> {log}
        """
        
rule cluster_sig_to_founders:
    message:
        """
        Find best cluster-founder match for a query signature
        """
    input:
        #namecheck=os.path.join(out_dir,".make_spreadsheet.touch"),
        #db = rules.find_founders.output.founders,
#        namecheck= lambda w: os.path.join(out_dir, f".{w.prefix}.{w.alphabet}-k{w.ksize}.mc{w.maxcontain}.make_spreadsheet.touch"),
        db = os.path.join(out_dir, "{prefix}.{alphabet}-k{ksize}.mc{maxcontain}.founders.sbt.zip"),
       # db = os.path.join(out_dir, "{prefix}.{alphabet}-k{ksize}.mc{maxcontain}.founders.siglist.csv"),
        query = lambda w: clusterMembers[w.name], # dictionary of member name :: sigfile
        #query = os.path.join(sigdir, "{name}.sig")
    output:
        # if doing a single query sig --> all founders, then write single file with info for this sig alone
        cluster_match = os.path.join(out_dir, "cluster_info", "{name}_x_{prefix}.{alphabet}-k{ksize}.mc{maxcontain}.best-founder.txt"),
    params:
        alpha_cmd= lambda w: config["alphabet_info"][w.alphabet]["alpha_cmd"],
    log: os.path.join(logs_dir, "find_founders", "{name}_x_{prefix}.{alphabet}-k{ksize}.mc{maxcontain}.cluster-to-founders.log" )
    benchmark: os.path.join(logs_dir, "find_founders", "{name}_x_{prefix}.{alphabet}-k{ksize}.mc{maxcontain}.cluster-to-founders.benchmark")
    conda: "envs/sourmash-maxcontain.yml"
    shell:
        """
        sourmash search --max-containment {input.query} {input.db} -o {output.cluster_match} \
                        --best-only --threshold 0.01 --ksize {wildcards.ksize} {params.alpha_cmd} 2> {log}
        """
        #clusters = os.path.join(out_dir, "clusters", "{prefix}.{alphabet}-{ksize}.mc{maxcontain}.cluster_info.csv")

rule gather_cluster_info:
    message:
        """
        gather cluster match info into a single csv
        """
    input:
        #namecheck=os.path.join(out_dir,".make_spreadsheet.touch"),
#        namecheck= lambda w: os.path.join(out_dir, f".{w.prefix}.{w.alphabet}-k{w.ksize}.mc{w.maxcontain}.make_spreadsheet.touch"),
        cluster_info=Checkpoint_MakePattern(os.path.join(out_dir, "cluster_info", "{name}_x_{prefix}.{alphabet}-k{ksize}.mc{maxcontain}.best-founder.txt")), # checkpoint makepattern should get member names
        founders=rules.find_founders.output.founders,
    output:
        os.path.join(out_dir, "{prefix}.{alphabet}-k{ksize}.mc{maxcontain}.cluster-info.csv")
    run:
        with open(str(output), "w") as outF:
        # edit this.  for each member sig, want line of:
        # cluster_name, founder_name, membersig_filename, member sig name, max containment / common hashes, etc
            for inF in input.cluster_info:
                outF.write(str(inF) + "\n")



