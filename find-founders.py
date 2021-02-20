#! /usr/bin/env python
"""
find-founders.py uses @ctb's sourmash-uniqify logic to perform iterative,
greedy clustering of a pile of sourmash signatures.

sourmash-uniqify is a practical alternative to a more principled clustering; see
https://github.com/ctb/2017-sourmash-cluster

The main goal here is just to obtain a set of cluster founders using iterative
`sourmash-uniqify` batches and similarity assessment. The list of non-founder ("member")
signatures can then be searched against all founder sigs to obtain the best cluster
placement for each member sig.

Authors:
  - N. Tessa Pierce Ward, github.com/bluegenes, ntpierce@gmail.com
  - C. Titus Brown, github.com/ctb/, titus@idyll.org

This code is under CC0.
"""
import sys
import argparse
import random
import csv
from collections import defaultdict, namedtuple

import sourmash
from sourmash.logging import notify

def load_sigs_from_list(siglistfiles, moltype, ksize):
    # input lists of signatures instead
    sigs = []
    for sl in siglistfiles:
        notify(f'loading from {sl}')
        sigfiles = sourmash.sourmash_args.load_file_list_of_signatures(sl)
        new_sigs = load_sigs(sigfiles, moltype, ksize, source_type= "input sigfile list")
        notify(f'...got {len(new_sigs)} signatures from {sl} siglist file.')
        sigs+=new_sigs
    return sigs

def load_sigs(sig_sources, moltype, ksize, source_type="input sigfiles"):
    siglist=[]
    for filename in sig_sources:
        if source_type != "input sigfile list":
            notify(f'loading from {filename}')
        m = 0
        for sig in sourmash.sourmash_args.load_file_as_signatures(filename,
                                           select_moltype=moltype,
                                           ksize=ksize):
            m += 1
            siglist.append((filename, sig))
        if source_type != "input sigfile list":
            notify(f'...got {m} signatures from {source_type}.')
    return siglist


def max_containment(sigA, sigB):
    c1 = sigA.contained_by(sigB)
    c2 = sigB.contained_by(sigA)
    return max(c1,c2)


def cluster_to_founders(founders, siglist, batch_n, pass_n, members):
    # splitting this out here may not be useful when single threaded, BUT
    # seems like it would be useful to split batches of remaining sigs, map to this batch of founders
    # then return lists that can be added to founders, members!

    # assign sigs to clusters via max containment to founder genomes
    notify(f'Attempting to cluster sigs to {len(founders)} new founders (pass {pass_n+1})')
    num_sigs = len(siglist)
    batch_size = len(founders)
    for n, (founder_from, founder) in enumerate(founders):
        leftover = []
        cluster_n = 0
        if n % 500 == 0:
            notify(f'batch {batch_n}: checking founder {n+1}/{str(batch_size)}. {str(len(siglist))}/{str(num_sigs)} sigs remaining')
        for (sig_from, sig) in siglist:
            maxcontain = max_containment(sig, founder)
            if not maxcontain >= args.threshold:
                leftover.append((sig_from, sig))
            else:
                members.append((sig_from, sig))
                cluster_n +=1
        if cluster_n:
            notify(f'    clustered {str(cluster_n)} signature(s) with founder sig {str(founder)[:30]}...')
        if leftover:
            siglist = leftover
            remaining = len(leftover)

            if remaining % 1000 == 0:
                print(f"{remaining} sigs left in this pass {pass_n+1}")
        else:
            break

    notify(f'{len(siglist)} signature(s) could not be assigned to existing clusters')
    return siglist, members


def get_new_founders_via_uniqify(siglist, batch_n, pass_n):
    '''
    use sourmash_uniqify code to build a new set of founders
    '''
    batch_size = len(siglist)
    notify(f'Finding new cluster founders from batch {batch_n} ({len(siglist)} sigs)')
    uniqify_pass_n = 0
    new_founders, new_members = [],[]
    while len(siglist):
        if uniqify_pass_n % 500 == 0:
            notify(f'batch {batch_n}: starting pass {uniqify_pass_n+1}. {str(len(siglist))}/{str(batch_size)} sigs remaining')
        # make the first one a founder; try to find matches; repeat.
        (founder_from, founder) = siglist.pop()
        new_founders.append((founder_from, founder))

        leftover = []
        cluster_n = 0
        for (sig_from, sig) in siglist:
            maxcontain = max_containment(sig, founder)
            if not maxcontain >= args.threshold:
                leftover.append((sig_from, sig))
            else:
                new_members.append((sig_from, sig))
                cluster_n +=1
        if cluster_n:
            notify(f'    clustered {str(cluster_n)} signature(s) with founder sig {str(founder)[:30]}...')

        siglist = leftover
        uniqify_pass_n += 1

    return new_founders, new_members


def main(args):
    batch_size = args.batch_size
    #load new sigs
    siglist=[]
    if args.signature_sources:
        siglist = load_sigs(args.signature_sources, args.moltype, args.ksize)
    if args.siglist:
        siglist+= load_sigs_from_list(args.siglist, args.moltype, args.ksize)

    notify(f'loaded {len(siglist)} new signatures total.')

    notify(f'setting random number seed to {args.seed} and shuffling input sigs')
    random.seed(args.seed)
    random.shuffle(siglist)

    founders, members = [],[]
    batch_n=0
    pass_n=0
    #if existing clusters, map to them first
    if args.seed_cluster_csv:
        notify(f'found existing input clusters.')
        founder_files = []
        # read in existing csv
        cluster_csv = csv.reader(open(args.seed_cluster_csv, "rt"))
        # grab 'filename' col from 'member_type' == "founder"
        for row in cluster_csv:
            if row["member_type"] == "founder":
                founder_files.append(row["filename"])
        # load in cluster sigs
        founders = load_sigs(founder_files, args.moltype, args.ksize, source_type="seed cluster founders")
        siglist, members  = cluster_to_founders(founders, siglist, batch_n, pass_n, members)
        pass_n+=1

    while siglist:
       # if unassigned sigs, uniqify to get new founders
       new_founders, new_members = get_new_founders_via_uniqify(siglist[:batch_size], batch_n, pass_n)
       founders += new_founders
       members += new_members
       siglist = siglist[batch_size:]
       batch_n+=1
       # cluster all sigs to list of new founders
       #siglist, members  = cluster_to_founders(founders, siglist, batch_n, pass_n, members)
       siglist, members  = cluster_to_founders(new_founders, siglist, batch_n, pass_n, members)
       pass_n +=1


    # write all founders, members
    prefix = args.prefix
    with open(f'{prefix}.founders.siglist.txt', 'wt') as fp:
        for (founder_from, founder) in founders:
            fp.write(founder_from + "\n")
    with open(f'{prefix}.members.siglist.txt', 'wt') as fp:
        for (member_from, member) in members:
            fp.write(member_from + "\n")


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--signature_sources', nargs='*',
                   help='signature files, directories, and sourmash databases')
    p.add_argument("--siglist",  action="append", help="provide list of signatures to assess")
    p.add_argument('-k', '--ksize', type=int, default=31)
    p.add_argument('--moltype', default='DNA')
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--threshold', type=float, default=0.05) # 0.2
    p.add_argument('--seed-cluster-csv')
    p.add_argument('--batch-size', type=int, default=5000)
    p.add_argument('--prefix', default='cluster',
                   help='output filename prefix (can include directories)')
    args = p.parse_args()
    if not any([args.signature_sources, args.siglist]):
        print("Please provide signatures via '--signature-sources' or '--siglist'")
        sys.exit(-1)
    sys.exit(main(args))
