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
    parser.add_argument("-a", "--pcawg_sample_sheet", dest="pcawg_sample_sheet", type=str, default="data/pcawg_sample_sheet.tsv", help="path to pcwag sample sheet file")
    parser.add_argument("-q", "--pcawg_sanger_qc", dest="pcawg_sanger_qc", type=str, default="data/pcawg_sanger_qc_metrics.jsonl", help="path to pcawg sanger qc file")

    parser.add_argument("-m", "--metadata_url", dest="metadata_url", type=str, default="https://song.rdpc.cancercollaboratory.org")
    parser.add_argument("-s", "--storage_url", dest="storage_url", type=str, default="https://score.rdpc.cancercollaboratory.org")
    parser.add_argument("-n", "--cpu_number", dest="cpu_number", type=str, default=4)
    parser.add_argument("-t", "--token", dest="token", type=str, required=True)
    args = parser.parse_args()

    song_dump = args.dump_path
    variant_calling_stats = {}

    #download data
    for subfolder in ['evaluate/sanger', 'evaluate/mutect2', 'evaluate/mutect2-bqsr']:
        download(song_dump, 'snv', args.token, args.metadata_url, args.storage_url, args.include, subfolder)
        download(song_dump, 'indel', args.token, args.metadata_url, args.storage_url, args.include, subfolder)

    data_dir = 'data/evaluate'
    evaluate_result = []
    for truvcf in glob.glob(os.path.join(data_dir, "sanger", "*-*", "*.vcf.gz"), recursive=True):
        result_dict = OrderedDict()
        studyId, donorId, sampleId, library, date_string, workflow, variant_type, evtype = os.path.basename(truvcf).split(".")[0:8]
        for sub in ['mutect2', 'mutect2-bqsr']:
            fname = '.'.join([studyId, donorId, sampleId, library, '*', 'gatk-mutect2', variant_type, evtype, 'vcf', 'gz']
            if not len(glob.glob(os.path.join(data_dir, sub, studyId, fname))) == 1: continue
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
