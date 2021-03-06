#!/usr/bin/env python3

import json
import os, io
import csv
import glob
from argparse import ArgumentParser
import sys
import subprocess
from collections import OrderedDict
from utils import report, download, run_cmd, get_dict_value, annot_vcf, union_vcf, vcf2tsv, region_query
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
    parser.add_argument("-e", "--mode", dest="mode", type=str, default="variant_calling")
    parser.add_argument("-f", "--force", dest="force", action='store_true')
    args = parser.parse_args()

    include = {}
    if not args.mode =='variant_calling':
        file_list = glob.glob('/'.join(['include', args.mode, '*'])) 
        for fl in file_list:    
            wf_name = os.path.splitext(os.path.basename(fl))[0]
            with open(fl, 'r') as f:
                for line in f:
                    if line.startswith('#'): continue
                    if not include.get(wf_name): include[wf_name] = set()
                    include[wf_name].add(line.rstrip())
    


    #download data and annotate
    process_flist = set() 
    for wf in ['sanger', 'mutect2']:
        subfolder = args.mode + '/' + wf
        if include and not include.get(wf): continue
        download_flist = download(args.dump_path, 'snv', wf, args.token, args.metadata_url, args.storage_url, include.get(wf, None), subfolder)
        process_flist.update(download_flist)
        download_flist = download(args.dump_path, 'indel', wf, args.token, args.metadata_url, args.storage_url, include.get(wf, None), subfolder)
        process_flist.update(download_flist)

        # annotate the vcf with gnomad AF
        data_dir = os.path.join("data", subfolder)
        annot_dir = os.path.join("data", subfolder+"_annot")
        
        annot_vcf(args.cpu_number, args.conf, data_dir, annot_dir, args.force)

        region_query(annot_dir)
    
    # union the result from different callers by donor
    data_dir = os.path.join("data", args.mode)
    union_dir = os.path.join("data", args.mode, 'union')
    union_vcf(data_dir, union_dir, process_flist)
    vcf2tsv(union_dir)

    bed_dir = os.path.join("data", "beds")
    regions = ['control_access', 'open_access', 'utr5', 'utr3', 'protein_coding_promoter', 'protein_coding_splice_site', 'lncRNA', 'lncRNA_promoter', 'lncRNA_splice_site', 'smallRNA', 'smallRNA_promoter', 'cds', 'exon']
    for region in regions:
        bed_file = os.path.join(bed_dir, '.'.join([region,'bed','gz']))
        region_dir = region_query(union_dir, region, args.force, bed_file)
        # generate tsv output for easy analysis
        vcf2tsv(region_dir)


if __name__ == "__main__":
    main()
