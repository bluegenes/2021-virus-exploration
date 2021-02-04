import os
import sys
import argparse
import pprint
from collections import defaultdict

import screed
import numpy as np
import pandas as pd

def main(args):
    # get basename for these sequences
    records = screed.open(args.fastafile)
    seqlens = defaultdict(int)
    for record in records:
        import pdb;pdb.set_trace()
        record_len = len(record.sequence)
        seqlens[record_len]+=1
    lenDF = pd.Series(seqlens).to_frame().reset_index()
    lenDF.to_csv(args.output_csv, header=["length", "count"], index=False)


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("fastafile")
    p.add_argument("output_csv")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
