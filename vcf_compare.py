#!/usr/bin/env python3

import json
import os, io
import csv
import glob
from argparse import ArgumentParser
import sys
import subprocess
from collections import OrderedDict
from utils import report, download, run_cmd, get_dict_value, annot_vcf
from evaluator import evaluate, countrecs
import copy
import numpy as np
import pandas as pd

mutect2_report_fields = {
    'studyId': 'studyId',
    'donorId': 'donorId',
    'sampleId': 'sampleId',
    'library_strategy': 'library_strategy',
    'evtype': 'evtype',
    'mutect2_total_query': 'mutect2.total_query',
    'mutect2_tp': 'mutect2.tp',
    'mutect2_fp': 'mutect2.fp',
    'mutect2_recall': 'mutect2.recall',
    'mutect2_precision': 'mutect2.precision',
    'mutect2_specificity': 'mutect2.specificity',
    'mutect2_balaccuracy': 'mutect2.balaccuracy',
    'mutect2_f1_score': 'mutect2.f1_score',
    'mutect2-bqsr_total_query': 'mutect2-bqsr.total_query',
    'mutect2-bqsr_tp': 'mutect2-bqsr.tp',
    'mutect2-bqsr_fp': 'mutect2-bqsr.fp',
    'mutect2-bqsr_recall': 'mutect2-bqsr.recall',
    'mutect2-bqsr_precision': 'mutect2-bqsr.precision',
    'mutect2-bqsr_specificity': 'mutect2-bqsr.specificity',
    'mutect2-bqsr_balaccuracy': 'mutect2-bqsr.balaccuracy',
    'mutect2-bqsr_f1_score': 'mutect2-bqsr.f1_score'
}


def main():
    parser = ArgumentParser()
    parser.add_argument("-d", "--dump_path", dest="dump_path", type=str, default="data/rdpc-song.jsonl", help="path to song dump jsonl file")
    parser.add_argument("-m", "--metadata_url", dest="metadata_url", type=str, default="https://song.rdpc.cancercollaboratory.org")
    parser.add_argument("-s", "--storage_url", dest="storage_url", type=str, default="https://score.rdpc.cancercollaboratory.org")
    parser.add_argument("-l", "--tool", dest="tool", type=str, default='som')
    parser.add_argument("-n", "--cpu_number", dest="cpu_number", type=int, default=4)
    parser.add_argument("-c", "--conf", dest="conf", type=str, default="conf/af_only_gnomad.conf")
    parser.add_argument("-t", "--token", dest="token", type=str, required=True)
    parser.add_argument("-e", "--mode", dest="mode", type=str, default="validation")
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
    for wf in ['sanger', 'mutect2', 'mutect2-bqsr']:
        subfolder = args.mode + '/' + wf
        if not include.get(wf): continue 
        download(args.dump_path, 'snv', args.token, args.metadata_url, args.storage_url, include.get(wf), subfolder)
        download(args.dump_path, 'indel', args.token, args.metadata_url, args.storage_url, include.get(wf), subfolder)

        # annotate the vcf with gnomad AF, get human readable table
        data_dir = os.path.join("data", subfolder)
        annot_dir = os.path.join("data", subfolder+"_annot_vcf")
        bed_dir = os.path.join("data", "beds")
        annot_vcf(args.cpu_number, args.conf, data_dir, annot_dir, bed_dir)

    # # union the result from different callers by donor
    # data_dir = os.path.join("data", args.mode)
    # union_dir = os.path.join("data", args.mode, 'union')
    # union_vcf(data_dir, union_dir)
    
    

'''
    data_dir = "data/evaluate"

    evaluate_result = []
    for fn in glob.glob(os.path.join(data_dir, "*_report", "*-*", "*.*"), recursive=True):
        evaluate_result.append(os.path.basename(fn))

    for truvcf in glob.glob(os.path.join(data_dir, "sanger", "*-*", "*.vcf.gz"), recursive=True):
        result_dict = OrderedDict()
        studyId, donorId, sampleId, library, date_string, workflow, variant_type, evtype = os.path.basename(truvcf).split(".")[0:8]
        for sub in ['mutect2', 'mutect2-bqsr']:
            fname = '.'.join([studyId, donorId, sampleId, library, '*', 'gatk-mutect2', variant_type, evtype, 'vcf', 'gz'])
            if not len(glob.glob(os.path.join(data_dir, sub, studyId, fname), recursive=True)) == 1: continue
            
            subvcf = glob.glob(os.path.join(data_dir, sub, studyId, fname))[0]
            fname = os.path.basename(subvcf)

            report_dir = os.path.join(data_dir, sub+'_report', studyId)
            if not os.path.exists(report_dir):
                os.makedirs(report_dir)
            output = os.path.join(report_dir, fname)

            if args.tool == 'som':
                if fname+'.som.stats.csv' in evaluate_result: continue
                cmd = 'export HGREF=/data/reference/GRCh38/GRCh38_hla_decoy_ebv.fa \
                    && docker run -it --rm -e HGREF -v `pwd`:/data pkrusche/hap.py \
                    /opt/hap.py/bin/som.py /data/%s /data/%s -o /data/%s' % (truvcf, subvcf, output+'.som')
                run_cmd(cmd)

            elif args.tool == 'evaluator':
                if fname+'.evaluator.stats.csv' in evaluate_result: continue
                chromlist = None
                if not subvcf.endswith('.vcf') and not subvcf.endswith('.vcf.gz'):
                    sys.stderr.write("submission VCF filename does not enc in .vcf or .vcf.gz\n")
                    sys.exit(1)

                if not os.path.exists(truvcf + '.tbi'):
                    sys.stderr.write("truth VCF does not appear to be indexed. bgzip + tabix index required.\n")
                    sys.exit(1)

                if evtype not in ('sv', 'snv', 'indel'):
                    sys.stderr.write("last arg must be either SV, SNV, or INDEL\n")
                    sys.exit(1)

                result = evaluate(subvcf, truvcf, vtype=evtype.upper(), ignorechroms=chromlist, truthmask=False)
                # count  = countrecs(subvcf, truvcf, vtype=evtype.upper(), ignorechroms=chromlist, truthmask=False)
                
                with open(output+'.evaluator.stats.csv', 'w') as f:
                    f.write("total.truth, total.query, tp, fp, recall, precision, specificity, balaccuracy, f1_score")
                    f.write(','.join(map(str, result)))

            else:
                sys.exit(0)

    evaluate_dict = {} 
    evaluate_dict_list = []       
    for fn in glob.glob(os.path.join(data_dir, "*_report", "*-*", "*."+args.tool+'.stats.csv'), recursive=True):
        sub = os.path.dirname(fn).split('/')[-2].split('_')[0]
        studyId, donorId, sampleId, library, date_string, workflow, variant_type, evtype = os.path.basename(fn).split(".")[0:8]
        result_list = {
            'studyId': studyId,
            'donorId': donorId,
            'sampleId': sampleId,
            'library_strategy': library,
            'evtype': evtype,
            'workflow': sub
        }
        if not evaluate_dict.get(sampleId+'_'+evtype): 
            evaluate_dict[sampleId+'_'+evtype] = {
            'studyId': studyId,
            'donorId': donorId,
            'sampleId': sampleId,
            'library_strategy': library,
            'evtype': evtype
        }
           
        result_dict = {}
        with open(fn, 'r') as f:
            dict_reader = csv.DictReader(f, delimiter=",")
            for row in dict_reader:
                if not row.get('type') in ['indels', 'SNVs']: continue
                for k in ['total.truth', 'total.query', 'tp', 'fp', 'fn', 'recall', 'recall_lower', 'recall_upper', 'recall2', 'precision', 'precision_lower', 'precision_upper', 'specificity', 'balaccuracy', 'f1_score']:
                    result_dict.update({k.replace('.', '_'): row.get(k, None)})
                try:
                    result_dict['specificity'] = 1.0 - float(row['fp']) / float(row['total.query'])
                    result_dict['balaccuracy'] = (float(row['recall']) + result_dict['specificity']) / 2.0
                    result_dict['f1_score'] = 2*float(row['recall'])*float(row['precision'])/(float(row['recall'])+float(row['precision']))
                except ZeroDivisionError:
                    continue
        
        result_list.update(result_dict)
        evaluate_dict_list.append(copy.deepcopy(result_list))
        evaluate_dict[sampleId+'_'+evtype].update({sub: result_dict})
        
    
    report_dir = "report"
    if not os.path.exists(report_dir):
        os.makedirs(report_dir)
    
    # generate v1 reports
    with open(os.path.join(report_dir, 'mutect2_evaluate_result.v1.json'), 'w') as f:
        for s in evaluate_dict_list:
            f.write(json.dumps(s)+'\n')

    report(evaluate_dict_list, os.path.join(report_dir, 'mutect2_evaluate_result.v1.txt'))

    # generate v2 reports
    evaluate_list = []
    with open(os.path.join(report_dir, 'mutect2_evaluate_result.v2.json'), 'w') as f:
        for k, s in evaluate_dict.items():
            f.write(json.dumps(s)+'\n')
            evaluate_list.append(get_dict_value(None, s, mutect2_report_fields))

    # generate tsv file
    report(evaluate_list, os.path.join(report_dir, 'mutect2_evaluate_result.v2.txt'))

'''

if __name__ == "__main__":
    main()
