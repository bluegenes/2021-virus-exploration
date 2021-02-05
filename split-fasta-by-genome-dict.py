#! /usr/bin/env python
import os
import sys
import argparse
from collections import defaultdict
import screed

def make_outdir(output_dirname):
    if not os.path.exists(output_dirname):
        try:
            os.makedirs(output_dirname)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

def main(args):
    print(f"Splitting {args.fasta} by genome. Writing files to {args.output_dir} \n")
    groupD = defaultdict(list)
    if not args.prefix:
        prefix = ""
    else:
        prefix = args.prefix + "-"
    make_outdir(args.output_dir)
    # loop through; split fasta into groups
    for n, record in enumerate(screed.open(args.fasta)):
        if n > 0 and n % 100000 == 0:
            print(f"working on {str(n)}th contig\n")
        name = record.name.split("|")[0]
        groupD[name].append(record)

    if args.output_names:
        names = open(args.output_names, "w")
    with open(args.output_csv, "w") as outcsv:
        outcsv.write("accession,filename\n")
        filenum=0
        for name, records in groupD.items():
            filenum+=1
            outfile = os.path.join(args.output_dir, f"{prefix}{name}.fa")
            outcsv.write(f"{name},{outfile}\n")
            if args.output_names:
                names.write(f"{name}\n")
            with open(outfile, "w") as out:
                for record in records:
                    out.write(f">{record.name}\n{record.sequence}\n")

    print(f"{str(n)} contigs written to {str(filenum)} group fasta files\n")

    if args.output_names:
        names.close()


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("fasta")
    p.add_argument("--output-dir", default= "")
    p.add_argument("--output-csv", required=True)
    p.add_argument("--output-names")
    p.add_argument("--prefix")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
