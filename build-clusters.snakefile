#####
# Use sourmash uniqify + matches to cluster a collection of sigs 
# snakemake -s cluster-sigs.snakefile --profile default -n
#####

import os
import pandas as pd

configfile: "conf-cluster.yml"
out_dir = config["output_dir"]
logs_dir = os.path.join(out_dir,"logs")

basename = config.get("basename", "pigeon1.0")



# to do: edits for this workflow
class Checkpoint_MakePattern:
    def __init__(self, pattern):
        self.pattern = pattern

    def get_member_names(self):
        #with open(f'{out_dir}/fastasplit/{basename}.names.txt', 'rt') as fp:
        #    names = [ x.rstrip() for x in fp ]
        #return names
        with open(f"os.path.join(out_dir, "{prefix}.{alphabet}-{ksize}.mc{maxcontain}.members.siglist.txt") as fp:

    #def get_member_siglist(self):
    #    # build global dictionary of member signame: sigfile
    #    global clusterMembers
    #    with open(f"os.path.join(out_dir, "{prefix}.{alphabet}-{ksize}.mc{maxcontain}.members.siglist.txt") as fp:
    #    with open(f'{out_dir}/fas/{basename}.lengths.txt', 'rt') as fp:
    #        genome_lengths = {}
    #        for line in fp:
    #            name, length = line.rstrip().split(',')
    #            genome_lengths[name] = length

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'check_csv'; this will trigger an
        # exception until that rule has been run.
        checkpoints.check_csv.get(**w)

        self.get_lengths()
        names = self.get_names()

        pattern = expand(self.pattern, genome=names, **w)
        return pattern


rule find_founders:
    input: 
        config["siglist"]
    output:
        founders = os.path.join(out_dir, "{prefix}.{alphabet}-{ksize}.mc{maxcontain}.founders.siglist.csv")
        members = os.path.join(out_dir, "{prefix}.{alphabet}-{ksize}.mc{maxcontain}.members.siglist.csv")
    conda:
        "envs/env-sourmash4.0.yml"
    log: os.path.join(logs_dir, "find_founders", "{prefix}.find-founders.log" )
    benchmark: os.path.join(logs_dir, "find_founders", "{prefix}.find-founders.benchmark" )
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


rule cluster_sig_to_founders:
    message:
        """
        Find best cluster-founder match for a query signature
        """
    input:
        namecheck=os.path.join(out_dir,".make_spreadsheet.touch"),
        db = rules.find_founders.output.founders,
        query = lambda w: sigD[w.name], # dictionary of name :: sigfile
        #query = os.path.join(sigdir, "{name}.sig")
    output:
        # if doing a single query sig --> all founders, then write single file with info for this sig alone
        cluster_match = os.path.join(out_dir, "cluster_info", "{name}.best-founder.txt"
    log: os.path.join(logs_dir, "find_founders", "{prefix}.{alphabet}-{ksize}.mc{maxcontain}.cluster-to-founders.log" )
    benchmark: os.path.join(logs_dir, "find_founders", "{prefix}.{alphabet}-{ksize}.mc{maxcontain}.cluster-to-founders.benchmark")
    shell:
        """
        sourmash search --max-containment --db {input.db} --query {input.query} -o {output.clusters} \
                                    --ksize {wildcards.ksize} --moltype {wildcards.alphabet} > {log} 2&>1
        """
        #clusters = os.path.join(out_dir, "clusters", "{prefix}.{alphabet}-{ksize}.mc{maxcontain}.cluster_info.csv")

rule gather_cluster_info:
    message:
        """
        gather cluster match info into a single csv
        """
    input:
        namecheck=os.path.join(out_dir,".make_spreadsheet.touch"),
        cluster_info=Checkpoint_MakePattern(os.path.join(out_dir, "cluster_info", "{name}.best-founder.txt")) # checkpoint makepattern should get member names
        founders=rules.find_founders.output.founders,
    output:
        os.path.join(out_dir, "{prefix}.cluster-info.csv")
    run:
        with open(str(output), "w") as outF:
        # edit this.  for each member sig, want line of:
        # cluster_name, founder_name, membersig_filename, member sig name, max containment / common hashes, etc
            for inF in input.cluster_info:
                outF.write(str(inF) + "\n")



