import os
import sys
import argparse

import screed
import pandas as pd


def find_genome_lengths_by_record(fastafile, outfile):
    #seqlens = defaultdict(int)

    with open(outfile, "w") as outF:
        with screed.open(fastafile) as records:
            for n, record in enumerate(records):
                if n % 10000 == 0:
                    print(f"... processing {n}th record, {record.name}\n")
                record_len = len(record.sequence)
                outF.write(f"{record.name},{record_len}\n")
                #seqlens[record.name] = record_len

    #return seqlens

def find_lengths_by_file(acc2files, outfile):

    with open(outfile, "w") as outF:
        #seqlens = defaultdict(int)
        num_files = 0
        for accession, filename in acc2files.items():
            num_files+=1
            if num_files % 10000 == 0:
                print(f"... processing {num_files}th file, {accession}\n")
            record_len=0
            with screed.open(filename) as records:
                for record in records:
                    record_len += len(record.sequence)
            outF.write(f"{accession},{record_len}\n")
            #seqlens[accession] = record_len

    #return seqlens




def main(args):
    # find fasta lengths for each genome or proteome
    if args.acc2fastafilescsv:
        acc2file = pd.read_csv(args.acc2fastafilescsv, header=0)
        accD = pd.Series(acc2file.filename.values,index=acc2file.accession).to_dict()
        find_lengths_by_file(accD, args.lengths_csv)
    elif args.fastalist:
        #sigh, some bespoke stuff for the protein filelist
        fasta_files = [x.rstrip() for x in open(args.fastalist, 'r')]
        accD = {}
        for ff in fasta_files:
            acc = os.path.basename(ff).rsplit(".proteins.fasta")[0].split("pigeon1.0-")[1]
            accD[acc] = ff
        find_lengths_by_file(accD, args.lengths_csv)
    elif args.fastafile:
        find_genome_lengths_by_record(args.fastafile, args.lengths_csv)
    else:
        print("please provide fasta file information via --fastafile or --acc2fastafilecsv")

    print("done!")



def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--acc2fastafilescsv", help="provide csv of accession,filename to get length for each file")
    p.add_argument("--fastalist", help="alternatively, if just have fastalist, input here and use bespoke filtering to get accession")
    p.add_argument("--fastafile", help="alternatively, if using just one fasta file, calculate length per record here by providing the fasta")
    p.add_argument("--lengths-csv", help="output lengths csv", required=True)
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)

