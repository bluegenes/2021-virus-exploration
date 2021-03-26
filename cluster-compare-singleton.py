import os
import sys
import argparse
import glob
import pprint

import pandas as pd

import screed
import sourmash
from sourmash.logging import notify
from sourmash.sourmash_args import load_file_as_signatures

from collections import defaultdict, namedtuple

CompareResult = namedtuple('CompareResult',
                           'comparison_name, anchor_name, ref_name, cluster_name, alphabet, ksize, scaled, jaccard, max_containment, anchor_containment, anchor_hashes, query_hashes, num_common')


def load_sigs_from_list(siglistfiles, moltype, ksize, sigdir=None):
    # input lists of signatures instead
    #sigs = []
    sigs = {}
    for sl in siglistfiles:
        notify(f'loading from {sl}')
        sigfiles = sourmash.sourmash_args.load_file_list_of_signatures(sl)
        new_sigs = load_sigs(sigfiles, moltype, ksize, source_type= "input sigfile list", sigdir=sigdir)
        notify(f'...got {len(new_sigs.keys())} signatures from {sl} siglist file.')
        #sigs+=new_sigs
        sigs.update(new_sigs)
    return sigs

def load_sigs(sig_sources, moltype, ksize, source_type="input sigfiles", sigdir=None):
    siglist=[]
    sigD = {}
    for filename in sig_sources:
        if not os.path.exists(filename) and sigdir:
            filename = os.path.join(sigdir, filename)
        if source_type != "input sigfile list":
            notify(f'loading from {filename}')
        m = 0
        for sig in sourmash.sourmash_args.load_file_as_signatures(filename,
                                           select_moltype=moltype,
                                           ksize=ksize):
            m += 1
            #siglist.append((filename, sig))
            sigD[str(sig)] =  sig
        if source_type != "input sigfile list":
            notify(f'...got {m} signatures from {source_type}.')
    #return siglist
    return sigD

def compare_sigs(sigA, sigB, comparison_name, cluster_name, alpha, ksize, scaled):
    sigA_numhashes = len(sigA.minhash.hashes)
    sigB_numhashes = len(sigB.minhash.hashes)
    intersect_numhashes = sigA.minhash.count_common(sigB.minhash)
    jaccard = sigA.jaccard(sigB)
    containA = sigA.contained_by(sigB)
    max_contain = sigA.max_containment(sigB)
    #max_contain = max(containA,containB)
    return CompareResult(comparison_name, str(sigA).split(" ")[0], str(sigB).split(" ")[0], cluster_name, alpha, ksize, scaled, jaccard, max_contain, containA, sigA_numhashes, sigB_numhashes, intersect_numhashes)

def main(args):
    ksize=args.ksize
    scaled=args.scaled
    alphabet=args.alphabet
    sigext = args.sig_extension
    if alphabet == "nucleotide":
        moltype = "DNA"
    else:
        moltype = alphabet

    sigD = {}
    # load all sigs
    if args.sigfiles:
        #siglist = load_sigs(args.sigfiles, args.moltype, args.ksize, sigdir=args.sigdir)
        sigD.update(load_sigs(args.sigfiles, moltype, args.ksize, sigdir=args.sigdir))
    if args.siglist:
        #siglist+= load_sigs_from_list(args.siglist, args.moltype, args.ksize, sigdir=args.sigdir)
        sigD.update(load_sigs_from_list(args.siglist, moltype, args.ksize, sigdir=args.sigdir))

    cluster_comparisons = []
    compareInfo = pd.read_csv(args.comparison_csv).set_index("cluster")
    compareInfo["cluster_members"] = compareInfo["cluster_members"].str.split(";")
    # loop through comparisons
    for n, cluster in enumerate(compareInfo.index):
        anchor_acc = compareInfo.at[cluster, "cluster_anchor"]
        if n !=0 and n % 50 == 0:
            print(f"... assessing {n}th cluster comparison, cluster name: {cluster}, anchor: {anchor_acc}\n")

        # select and load anchor sig
        anchor_sig = sigD[anchor_acc]

        # iterate through comparison sigs
        compare_accs = compareInfo.at[cluster, "cluster_members"]
        for compare_acc in compare_accs:
            if compare_acc != anchor_acc:
                # select and load comparison sig
                compare_sig = sigD[compare_acc]
                comparison = compare_sigs(anchor_sig, compare_sig, f"{anchor_acc}_x_{compare_acc}", cluster, alphabet, ksize, scaled)
                cluster_comparisons.append(comparison)

    # convert path comparison info to pandas dataframe
    comparisonDF = pd.DataFrame.from_records(cluster_comparisons, columns = CompareResult._fields)

    # print to csv
    comparisonDF.to_csv(args.output_csv, index=False)
    print(f"done! taxon comparison info written to {args.output_csv}")

def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--comparison-csv")
    p.add_argument("--lineages-csv", default="gtdb-r95-reps.lineages.protein-filenames.reordered.csv")
    p.add_argument('--sigfiles', nargs='*',help='signature files, directories, and sourmash databases')
    p.add_argument("--siglist", action="append", help="provide list of signatures to assess")
    p.add_argument("--sigdir", default="gtdb95-evolpaths/signatures")
    p.add_argument("--sig-extension", default=".sig")
    p.add_argument("--alphabet", default="protein")
    p.add_argument("--ksize", default=10, type=int)
    p.add_argument("--scaled", default=100, type=int)
    p.add_argument("--output-csv", required=True)
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
