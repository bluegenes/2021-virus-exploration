import os
import pandas as pd
#####
#
#Compare distance estimation on pairs of unaligned nucleotide sequences
####

# Build pigeon reference seqs 
# 1. DNA sequence --> dna sigs
# 2. DNA sequence --> split fastas --> PRODIGAL translate
#  - https://github.com/hyattpd/prodigal/wiki/advice-by-input-type
#We have never extensively examined exactly how much sequence Prodigal requires in order to produce good results. It is our general feeling that 20kbp-100kbp is too short to gather enough statistics to predict genes well, but we have no actual data to support this claim. If you feel comfortable running 50kbp sequences in normal mode, and you're happy with the results you're getting, then, by all means, continue to do so. Start sites, especially, require a lot of data (100kbp may only be ~80 or so start sites in the training set, since Prodigal only trains on the starts of longer ORFs). We recommend 500kbp+ to get ideal 5' predictions, although this is likely on the conservative side. If your sequence is in this nebulous range (20kbp to 100kbp for 3' predictions, 100kbp to 500kbp for 5' predictions), you may try running both normal and anonymous mode and manually inspecting the results for any differences.

#Prodigal contains no special routines to deal with viruses. As such, it cannot handle certain phenomena that occur sometimes in viruses, such as translational frame shifts. Viruses should generally be analyzed as above, with short genomes analyzed in anonymous mode and longer ones in normal mode.


# 3. Check/report % coding for each genome
# 4. Prodigal sequence --> protein sigs

configfile: "conf-prot.yml"
out_dir = config["output_dir"]
logs_dir = out_dir + "/logs"

basename = config.get("basename", "pigeon1.0")
genomes_fasta = config["genomes_fasta"]
fasta_dir = config.get("fasta_dir", "")

# checkpoint code to read count and specify all the outputs
class Checkpoint_MakePattern:
    def __init__(self, pattern):
        self.pattern = pattern

    def get_names(self):
        with open(f'{out_dir}/fastasplit/{basename}.names.txt', 'rt') as fp:
            names = [ x.rstrip() for x in fp ]
        return names

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'check_csv'; this will trigger an
        # exception until that rule has been run.
        checkpoints.check_csv.get(**w)

        # the magic, such as it is, happens here: we create the
        # information used to expand the pattern, using arbitrary
        # Python code.
        names = self.get_names()

        pattern = expand(self.pattern, genome=names, **w)
        return pattern


rule all: 
    input: 
        expand(os.path.join(out_dir, "compare", "{name}.prodigal.fastalist.txt"), name=basename),
        expand(os.path.join(out_dir, "compare", "{name}.prodigal.siglist.txt"), name=basename),

rule split_fasta:
    input: config["genomes_fasta"]
    output: 
        csv=os.path.join(out_dir, "{name}.fastasplit.csv"),
        names=os.path.join(out_dir, "fastasplit", "{name}.names.txt"),
    params:
        outdir = os.path.join(out_dir, "fastasplit"),
    log: os.path.join(logs_dir, "{name}.fastasplit.log")
    benchmark: os.path.join(logs_dir, "{name}.fastasplit.benchmark")
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=120,
    shell:
        """
        python split-fasta-by-genome-dict.py {input} \
               --output-dir {params.outdir} \
               --output-csv {output.csv} \
               --prefix {wildcards.name} \
               --output-names {output.names} > {log} 2>&1
        """

checkpoint check_csv:
    input: 
        os.path.join(out_dir, "fastasplit", f"{basename}.names.txt"),
        #rules.split_fasta.output.names
    output: touch(f"{out_dir}/.make_spreadsheet.touch")


rule prodigal_translate:
    input:
        os.path.join(out_dir, "fastasplit", "{name}-{genome}.fa"),
    output: 
        gff=os.path.join(out_dir, "prodigal", "{name}-{genome}.genes.gff3"),
        proteins=os.path.join(out_dir, "prodigal", "{name}-{genome}.proteins.fasta"),
        #stats=os.path.join(out_dir, "prodigal", "{name}-{genome}.stats.txt")
    log:
        os.path.join(logs_dir, "prodigal","{name}-{genome}.prodigal.log")
    benchmark:
        os.path.join(logs_dir, "prodigal","{name}-{genome}.prodigal.benchmark")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *2000,
        runtime=120,
    conda: "envs/prodigal-env.yml"
    shell:
        """
        prodigal -i {input} -o {output.gff} -a {output.proteins} \
                  -f "gff" > {log} 2>&1
        """
        # --summ_file {output.stats}


def build_sketch_params(output_type):
    sketch_cmd = ""
    for alpha in ["protein", "dayhoff", "hp"]:
        if alpha in config["alphabet_info"].keys():
            ksizes = config["alphabet_info"][alpha]["ksizes"]
            scaled = min(config["alphabet_info"][alpha]["scaled"])
        sketch_cmd += " -p " + alpha + ",k=" + ",k=".join(map(str, ksizes)) + f",scaled={str(scaled)}" + ",abund"
    return sketch_cmd

rule sourmash_sketch_prodigal_input:
    input:
        os.path.join(out_dir, "prodigal", "{name}-{genome}.proteins.fasta")
    output:
        os.path.join(out_dir, "prodigal", "signatures", "{name}-{genome}.prodigal.sig"),
    params:
        sketch_params = build_sketch_params("protein"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *2000,
        runtime=120,
    log: os.path.join(logs_dir, "sourmash_sketch_prot_input", "{name}-{genome}.sketch.log")
    benchmark: os.path.join(logs_dir, "sourmash_sketch_prot_input", "{name}-{genome}.sketch.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash sketch protein {params.sketch_params} --name {wildcards.genome:q} -o {output} {input} 2> {log}
        """


localrules: write_prodigal_siglist, write_prodigal_fastalist
rule write_prodigal_siglist:
    input:
        namecheck=os.path.join(out_dir,".make_spreadsheet.touch"),
        sigs=Checkpoint_MakePattern(os.path.join(out_dir, "prodigal/signatures", "{name}-{genome}.prodigal.sig"))
        #lambda w: expand(os.path.join(out_dir, "prodigal", "signatures", "{name}.prodigal.sig"), name = dataset_samples[w.dataset]) 
    output: os.path.join(out_dir, "compare", "{name}.prodigal.siglist.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input.sigs:
                outF.write(str(inF) + "\n")

rule write_prodigal_fastalist:
    input:
        namecheck=os.path.join(out_dir,".make_spreadsheet.touch"),
        fasta=Checkpoint_MakePattern(os.path.join(out_dir, "prodigal", "{name}-{genome}.proteins.fasta"))
        #lambda w: ancient(expand(os.path.join(out_dir, "prodigal", "{name}.proteins.fasta"), name = dataset_samples[w.dataset])),
    output: os.path.join(out_dir, "compare", "{name}.prodigal.fastalist.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input.fasta:
                outF.write(str(inF) + "\n")

