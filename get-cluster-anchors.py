import os
import sys
import argparse
import glob
import pprint

import pandas as pd

# build anchor comparison csv, each line containing the cluster name, cluster anchor and ;-separated cluster members

def find_pigeon_name_based_on_original_name_substring(row):
    name = row["cluster_member"]
    convert_name = namemap[namemap["orig"].str.contains(name, case=False)]["pigeon"]
    print(f"{name}: {convert_name}")
    if not convert_name.empty:
        name_list=convert_name.tolist()
        name = max(name_list, key=len)
    return name

def main(args):
    infoDF = pd.read_csv(args.info_csv)
    # split cluster members into a list
    infoDF["cluster_member"] = infoDF["Members"].str.split(";")
    # explode the cluster member list into a long-form dataframes
    infoDF = infoDF.explode("cluster_member")
    # map all to pigeon names
    global namemap
    namemap = pd.read_csv(args.namemap_csv, names=["pigeon", "orig", "reference"], header=0)
    infoDF["cluster_member"] = infoDF["cluster_member"].str.replace("~", "_") # replace ~ which exist for some reason
    # drop environment info cols
    infoDF.drop(columns=["Freshwaterenvironment", "Marineenvironment", "SPRUCEpeatlandenvironment", "Otherpeatlandenvironment", "Othersoilenvironment"], inplace=True)
    infoDF["pigeon_names"] = infoDF.apply(find_pigeon_name_based_on_original_name_substring, axis=1)
    # groupby Clusternumber and choose n0th as cluster anchor (all clusters have at least two members)
    infoDF["cluster_anchor"] = infoDF.groupby("Clusternumber", as_index=False).nth(0)["pigeon_names"]
    # aggregate pigeon names into ;-separated list
    infoDF = infoDF.groupby(["Clusternumber", "cluster_anchor"], as_index=False).agg({'pigeon_names':lambda x: ";".join(list(x))})
    # reorder columns
    infoDF = infoDF[["Clusternumber", "cluster_anchor", "pigeon_names"]]
    # rename columns and drop duplicates
    infoDF = infoDF.rename(columns={"Clusternumber": "cluster", "pigeon_names": "cluster_members"})
    infoDF.drop_duplicates(inplace=True)
    # write csv:
    infoDF.to_csv(args.output_csv, index=False)

def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--info-csv", default = "pigeon1.0.vContact.clustering.csv")
    p.add_argument("--output-csv", default = "pigeon1.0.vContact.clustering.anchors.csv")
    p.add_argument("--namemap-csv", default = "pigeon1.0.pubseq.namemap.csv")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)

