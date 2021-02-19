#! /usr/bin/env python
"""
Uses @ctb's sourmash-uniqify + maximum containment assessement to perform
iterative, greedy clustering of a pile of sourmash signatures.

sourmash-uniqify is a practical alternative to a more principled clustering; see
https://github.com/ctb/2017-sourmash-cluster

Authors:
  - C. Titus Brown, github.com/ctb/, titus@idyll.org
  - N. Tessa Pierce Ward, github.com/bluegenes, ntpierce@gmail.com

This code is under CC0.
"""
import sys
import argparse
import random
import csv
from collections import defaultdict, namedtuple

import sourmash
from sourmash import load_file_as_signatures #, load_file_list_of_signatures
from sourmash.logging import notify


def load_sigs(sig_sources, moltype, ksize, source_type="input sigfiles"):
    siglist=[]
    for filename in sig_sources:
        # to do -- enable list of sigfiles here!
        notify(f'loading from {filename}')
        m = 0
        for sig in load_file_as_signatures(filename,
                                           select_moltype=moltype,
                                           ksize=args):
            m += 1
            siglist.append((filename, sig))
        notify(f'...got {m} signatures from {source_type}.')
    return siglist


def max_containment(sigA, sigB):
    c1 = sigA.contained_by(sigB)
    c2 = sigB.contained_by(sigA)
    return max(c1,c2)


#def cluster_to_founders(founders, query_sigs, clusterInfo, cluster_summary, batch_n, pass_n):
def cluster_to_founders(founders, query_sigs, batch_n, pass_n):
    # assign sigs to clusters via max containment to founder genomes
    notify(f'Attempting to cluster sigs to {len(founders)} current founders (pass {pass_n+1})')
    for (founder_from, founder) in founders:
        cluster = []
        leftover = []
        for (sig_from, sig) in siglist:
            maxcontain = max_containment(sig, founder)
            #if sig.similarity(founder) >= args.threshold:
            if not maxcontain >= args.threshold:
                #clusterInfo[founder].append((sig_from, sig))
                #cluster.append((sig_from, sig))
                #cluster_summary.append((sig_from, sig, batch_n, pass_n, 'member'))
            #else:
                leftover.append((sig_from, sig))
        if cluster:
            notify(f'clustered {len(cluster)} signature(s) with founder sig {str(founder)[:30]}...')
            #clusterInfo[str(founder)].extend(cluster)
        else:
            notify(f'No new members for cluster from founder sig {str(founder)[:30]}...')
        siglist = leftover

    notify(f'{len(siglist)} signature(s) could not be assigned to existing clusters')
    #return siglist, clusterInfo, cluster_summary
    return siglist


def get_new_founders_via_uniqify(siglist, clusterInfo, batch_n, pass_n):
    '''
    use sourmash_uniqify code to build a new set of founders
    '''
    notify(f'Finding new cluster founders from batch {batch_n} ({len(siglist)} sigs)')
    uniqify_pass_n = 0
    new_founders = []
    while len(siglist):
        notify(f'batch {batch_n}: starting pass {uniqify_pass_n+1}')
        # make the first one a founder; try to find matches; repeat.
        (founder_from, founder) = siglist.pop()
        new_founders.append((founder_from, founder))
        #cluster_summary.append((founder_from, founder, batch_n, uniqify_pass_n, 'founder'))

        cluster = []
        leftover = []
        for (sig_from, sig) in siglist:
            maxcontain = max_containment(sig, founder)
            if not maxcontain >= args.threshold:
            #if sig.similarity(founder) >= args.threshold:
                #cluster.append((sig_from, sig))
                #cluster_summary.append((sig_from, sig, batch_n, uniqify_pass_n, 'member'))
            #else:
                leftover.append((sig_from, sig))

        if cluster:
            notify(f'clustered {len(cluster)} signature(s) with founder sig {str(founder)[:30]}...')
            #clusterInfo[str(founder)].extend(cluster)

            #prefix = f'{args.prefix}.cluster.{pass_n}'
            #with open(f'{prefix}.founder.sig', 'wt') as fp:
            #    sourmash.save_signatures([founder], fp)
            #with open(f'{prefix}.cluster.sig', 'wt') as fp:
            #    cluster_sigs = [ x[1] for x in cluster ]
            #    sourmash.save_signatures(cluster_sigs, fp)

            #print(f'saved founder and {len(cluster)} signatures to {prefix}.*')
        else:
            notify(f'founder sig {str(founder)[:30]}... is a singleton.')

            #prefix = f'{args.prefix}.cluster.{pass_n}'
            #with open(f'{prefix}.founder.sig', 'wt') as fp:
            #    sourmash.save_signatures([founder], fp)
            #print(f'saved singleton signature to {prefix}.*')

        siglist = leftover
        #pass_n += 1
        uniqify_pass_n += 1

    return new_founders #, clusterInfo, cluster_summary


def main(args):
    batch_size = args.batch_size
    #load new sigs
    siglist = load_sigs(args.signature_sources, args.moltype, args.ksize)
    notify(f'loaded {len(siglist)} new signatures total.')

    notify(f'setting random number seed to {args.seed} and shuffling input sigs')
    random.seed(args.seed)
    random.shuffle(siglist)

    cluster_summary, founders = [],[]
    clusterInfo = defaultdict(list)
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
        #siglist, clusterInfo, cluster_summary = cluster_to_founders(founders, siglist, clusterInfo, cluster_summary, batch_n, pass_n)
        siglist  = cluster_to_founders(founders, siglist, batch_n, pass_n) #clusterInfo, cluster_summary, batch_n, pass_n)
        pass_n+=1

    while siglist:
       # if unassigned sigs, uniqify to get new founders
       #new_founders, clusterInfo, cluster_summary = get_new_founders_via_uniqify(siglist[:batch_size], clusterInfo, cluster_summary, batch_n, pass_n)
       new_founders  = get_new_founders_via_uniqify(siglist[:batch_size], batch_n, pass_n) #, clusterInfo, cluster_summary, batch_n, pass_n)
       founders += new_founders
       batch_n+=1
       # cluster all sigs to full list of founders
       #siglist, clusterInfo, cluster_summary = cluster_to_founders(founders, siglist, clusterInfo, cluster_summary, batch_n, pass_n)
       siglist  = cluster_to_founders(founders, siglist, batch_n, pass_n) #, clusterInfo, cluster_summary, batch_n, pass_n)
       pass_n +=1


    # this script is really about using a greedy alg to find a set of founders. We need to re-map all to founders after this
    # (to get best matches, not just first matches), so all we really need is the list of founder sigs

    # write all founders
    with open(f'{prefix}.founders.siglist', 'wt') as fp:
        for (founder_from, founder) in founders:
            fp.write(founder_from + "\n")
        #sourmash.save_signatures([founder], fp)

    # write output summary spreadsheet
    #headers = ['origin_path', 'name', 'filename', 'md5sum', 'cluster', 'member_type']
    #csv_name = f'{args.prefix}.summary.csv'

    #with open(csv_name, 'wt') as fp:
    #    w = csv.writer(fp)
    #    w.writerow(headers)

     #   for (origin_path, sig, cluster_n, member_type) in cluster_summary:
     #       name = str(sig)
     #       filename = sig.filename
     #       md5sum = sig.md5sum()

      #      w.writerow([origin_path, name, filename, md5sum, cluster_n, member_type])

    #notify(f"wrote {len(cluster_summary)} entries to clustering summary at '{csv_name}'")


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('signature_sources', nargs='+',
                   help='signature files, directories, and sourmash databases')
    p.add_argument('-k', '--ksize', type=int, default=31)
    p.add_argument('--moltype', default='DNA')
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--threshold', type=float, default=0.2)
    p.add_argument('--seed-cluster-csv')
    p.add_argument('--batch-size', type=int, default=5000)
    p.add_argument('--prefix', default='cluster',
                   help='output filename prefix (can include directories)')
    args = p.parse_args()
    sys.exit(main(args))
