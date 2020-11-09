#!/usr/bin/env python3

import json
import os, io
import csv
import glob
from argparse import ArgumentParser
import sys
import subprocess
from collections import OrderedDict
from utils import report, download, run_cmd, get_dict_value, annot_vcf, union_vcf, vcf2tsv
from evaluator import evaluate, countrecs
import copy
import numpy as np
import pandas as pd


def main():
    parser = ArgumentParser()
    parser.add_argument("-d", "--dump_path", dest="dump_path", type=str, default="data/rdpc-song.jsonl", help="path to song dump jsonl file")
    parser.add_argument("-m", "--metadata_url", dest="metadata_url", type=str, default="https://song.rdpc.cancercollaboratory.org")
    parser.add_argument("-s", "--storage_url", dest="storage_url", type=str, default="https://score.rdpc.cancercollaboratory.org")
    parser.add_argument("-l", "--tool", dest="tool", type=str, default='som')
    parser.add_argument("-n", "--cpu_number", dest="cpu_number", type=int, default=4)
    parser.add_argument("-c", "--conf", dest="conf", type=str, default="conf/af_only_gnomad.conf")
    parser.add_argument("-t", "--token", dest="token", type=str, required=True)
    parser.add_argument("-e", "--mode", dest="mode", type=str, default="argo_40donors")
    args = parser.parse_args()

    include = {}
    file_list = glob.glob('/'.join(['include', args.mode, '*'])) 
    for fl in file_list:    
        wf_name = os.path.splitext(os.path.basename(fl))[0]
        with open(fl, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                if not include.get(wf_name): include[wf_name] = set()
                include[wf_name].add(line.rstrip())


    #download data and annotate
    for wf in ['sanger', 'mutect2']:
        subfolder = args.mode + '/' + wf
        if not include.get(wf): continue 
        download(args.dump_path, 'snv', args.token, args.metadata_url, args.storage_url, include.get(wf), subfolder)
        download(args.dump_path, 'indel', args.token, args.metadata_url, args.storage_url, include.get(wf), subfolder)

        # annotate the vcf with gnomad AF, get human readable table
        data_dir = os.path.join("data", subfolder)
        annot_dir = os.path.join("data", subfolder+"_annot_vcf")
        #bed_dir = os.path.join("data", "beds")
        annot_vcf(args.cpu_number, args.conf, data_dir, annot_dir)

    # union the result from different callers by donor
    data_dir = os.path.join("data", args.mode)
    union_dir = os.path.join("data", args.mode, 'union')
    union_vcf(data_dir, union_dir)

    # generate tsv output for easy analysis
    vcf2tsv(union_dir)


if __name__ == "__main__":
    main()
