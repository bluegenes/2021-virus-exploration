import os
import sys
import argparse
from collections import defaultdict, namedtuple

import screed
import sourmash
import numpy as np
import pandas as pd

SigInfo = namedtuple('SigInfo','name, ksize, scaled, num_hashes, genome_length')

def find_genome_lengths_single(fastafile):
    seqlens = defaultdict(int)

    with screed.open(fastafile) as records:
        for n, record in enumerate(records):
            if n % 10000 == 0:
                print(f"... processing {n}th record, {record.name}\n")
            record_len = len(record.sequence)
            seqlens[record.name] = record_len

    return seqlens


def main(args):
    scaled_vals = []
    if args.scaled:
        scaled_vals = args.scaled

    # find fasta lengths for each genome or proteome
    if args.length_csv:
        lenDF = pd.read_csv(args.length_csv, names=["name", "length"])
        genome_lengths = lenDF.set_index('name')['length'].to_dict()
    elif args.fastafile:
        genome_lengths = find_genome_lengths_single(args.fastafile)

    # load file list of sigs
    sigfiles = sourmash.sourmash_args.load_file_list_of_signatures(args.siglist)
    sigInfoList = []
    for sigF in sigfiles:
        # load sigs from each sigfile
        sigs = sourmash.sourmash_args.load_file_as_signatures(sigF)
        for n, sig in enumerate(sigs):
            # get signature information
            name = str(sig)
            genome_len = genome_lengths[name]
            if n % 10000 == 0:
                print(f"... processing {n}th sig, {name}\n")
            ksize = sig.minhash.ksize
            # first do existing scaled val
            scaled = sig.minhash.scaled
            num_hashes = len(sig.minhash.hashes)
            # store signature info
            sigInfoList.append(SigInfo(name=name, ksize=ksize, scaled=scaled, num_hashes=num_hashes, genome_length=genome_len))
            # now downsample to additional scaled vals, if desired
            for sc in scaled_vals:
                if sc < scaled:
                    print(f"Can't downsample: desired scaled {sc} is smaller than original scaled, {scaled}. Removing scaled {sc}...")
                    scaled_vals.remove(sc)
                    continue
                # store signature info
                num_hashes = len(sig.minhash.downsample(scaled=sc).hashes)
                sigInfoList.append(SigInfo(name=name, ksize=ksize, scaled=sc, num_hashes=num_hashes, genome_length=genome_len))
        print("total number sigs processed: {n}")

    # convert signature info to pandas dataframe
    sigInfoDF = pd.DataFrame.from_records(sigInfoList, columns = SigInfo._fields)

    # print to csv
    sigInfoDF.to_csv(args.output_csv, index=False)
    print("yay!")



def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--siglist", help="provide list of signatures to assess", required=True)
    p.add_argument("--output-csv", help="provide output filename for stats", required=True)
    p.add_argument("--length-csv", help="provide a csv of 'signame, fastalen' here")
    p.add_argument("--fastafile", help="alternatively, if using just one fasta file, calculate length per record here by providing the fasta")
    p.add_argument("-s", "--scaled", action="append", type=int, help= "provide additional scaled values for downsampling")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)

