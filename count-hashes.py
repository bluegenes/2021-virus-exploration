import os
import sys
import argparse
from collections import defaultdict, namedtuple

import screed
import sourmash
import numpy as np
import pandas as pd

SigInfo = namedtuple('SigInfo','name, ksize, scaled, num_hashes, genome_length')

def find_genome_lengths(fastafile):
    seqlens = defaultdict(int)

    with screed.open(fastafile) as records:
        for n, record in enumerate(records):
            if n % 10000 == 0:
                print(f"... processing {n}th record, {record.name}\n")
            record_len = len(record.sequence)
            seqlens[record.name] = record_len

    return seqlens


def main(args):
    # find genome length for each genome
    genome_lengths = find_genome_lengths(args.fastafile)

    # load each signature
    sigs = sourmash.sourmash_args.load_file_as_signatures(args.sigfile)
    sigInfoList = []
    for sig in sigs:
        # get signature information
        name = sig.name
        print(f"Working on {name}")
        ksize = sig.minhash.ksize
        scaled = sig.minhash.scaled
        num_hashes = len(list(sig.minhash.hashes))
        genome_len = genome_lengths[name]
        # store signature info
        sigInfoList.append(SigInfo(name=name, ksize=ksize, scaled=scaled, num_hashes=num_hashes, genome_length=genome_len))

    # convert signature info to pandas dataframe
    sigInfoDF = pd.DataFrame.from_records(sigInfoList, columns = SigInfo._fields)

    # print to csv
    sigInfoDF.to_csv(args.output_csv, index=False)
    print("yay!")



def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("sigfile")
    p.add_argument("output_csv")
    p.add_argument("--fastafile")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)

