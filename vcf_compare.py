#!/usr/bin/env python3

import json
import os
import csv
import glob
from argparse import ArgumentParser
import sys
import subprocess
from collections import OrderedDict
from utils import report, download
from evaluator import evaluate, countrecs
import copy


def main():
    parser = ArgumentParser()
    parser.add_argument("-d", "--dump_path", dest="dump_path", type=str, default="data/rdpc-song.jsonl", help="path to song dump jsonl file")
    parser.add_argument("-i", "--include", dest="include", action='store_true', help="include list to work on")
    parser.add_argument("-m", "--metadata_url", dest="metadata_url", type=str, default="https://song.rdpc.cancercollaboratory.org")
    parser.add_argument("-s", "--storage_url", dest="storage_url", type=str, default="https://score.rdpc.cancercollaboratory.org")
    parser.add_argument("-n", "--cpu_number", dest="cpu_number", type=str, default=4)
    parser.add_argument("-t", "--token", dest="token", type=str, required=True)
    args = parser.parse_args()

    if args.include:
        file_list = glob.glob("include/*.*") 
        include = {}
        for fl in file_list:    
            wf_name = os.path.splitext(os.path.basename(fl))
            with open(fl, 'r') as f:
                for line in f:
                    if line.startswith('#'): continue
                    if not include.get(wf_name): include[wf_name] = set()
                    include[wf_name].add(line.rstrip())
    else:
        include = None

    #download data
    for wf in ['sanger', 'mutect2', 'mutect2-bqsr']:
        subfolder = 'evaluate/'+wf
        download(args.dump_path, 'snv', args.token, args.metadata_url, args.storage_url, include.get(wf, None), subfolder)
        download(args.dump_path, 'indel', args.token, args.metadata_url, args.storage_url, include.get(wf, None), subfolder)

    data_dir = 'data/evaluate'
    evaluate_result = []
    for truvcf in glob.glob(os.path.join(data_dir, "sanger", "*-*", "*.vcf.gz"), recursive=True):
        result_dict = OrderedDict()
        studyId, donorId, sampleId, library, date_string, workflow, variant_type, evtype = os.path.basename(truvcf).split(".")[0:8]
        for sub in ['mutect2', 'mutect2-bqsr']:
            fname = '.'.join([studyId, donorId, sampleId, library, '*', 'gatk-mutect2', variant_type, evtype, 'vcf', 'gz'])
            print(fname)
            if not len(glob.glob(os.path.join(data_dir, sub, studyId, fname), recursive=True)) == 1: continue
            subvcf = glob.glob(os.path.join(data_dir, sub, studyId, fname))[0]

            chromlist = None
            if not subvcf.endswith('.vcf') and not subvcf.endswith('.vcf.gz'):
                sys.stderr.write("submission VCF filename does not enc in .vcf or .vcf.gz\n")
                sys.exit(1)

            if not os.path.exists(truvcf + '.tbi'):
                sys.stderr.write("truth VCF does not appear to be indexed. bgzip + tabix index required.\n")
                sys.exit(1)

            if evtype not in ('SV', 'SNV', 'INDEL'):
                sys.stderr.write("last arg must be either SV, SNV, or INDEL\n")
                sys.exit(1)

            result = evaluate(subvcf, truvcf, vtype=evtype, ignorechroms=chromlist, truthmask=False)
            count  = countrecs(subvcf, truvcf, vtype=evtype, ignorechroms=chromlist, truthmask=False)
            
            result_dict = {
                'studyId': studyId,
                'donorId': donorId,
                'sampleId': sampleId,
                'library_strategy': library,
                'workflow': sub,
                'evtype': evtype,
                'sensitivity': result[0],
                'specificity': result[1],
                'accuracy': result[2],
                'mutation_count': count
            }
            evaluate_result.append(copy.deepcopy(result_dict)) 
    
    with open(os.path.join(report_dir, 'mutect2_evaluate_result.json'), 'w') as f:
        for s in evaluate_result:
            f.write(json.dumps(s)+'\n')

    # generate tsv file
    report(evaluate_result, os.path.join(report_dir, 'mutect2_evaluate_result.txt'))


if __name__ == "__main__":
    main()