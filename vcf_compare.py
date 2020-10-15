#!/usr/bin/env python3

import json
import os
import csv
import glob
from argparse import ArgumentParser
import sys
import subprocess
from collections import OrderedDict
from utils import report, download, run_cmd
from evaluator import evaluate, countrecs
import copy


def main():
    parser = ArgumentParser()
    parser.add_argument("-d", "--dump_path", dest="dump_path", type=str, default="data/rdpc-song.jsonl", help="path to song dump jsonl file")
    parser.add_argument("-m", "--metadata_url", dest="metadata_url", type=str, default="https://song.rdpc.cancercollaboratory.org")
    parser.add_argument("-s", "--storage_url", dest="storage_url", type=str, default="https://score.rdpc.cancercollaboratory.org")
    parser.add_argument("-c", "--tool", dest="tool", type=str, default='som')
    parser.add_argument("-n", "--cpu_number", dest="cpu_number", type=int, default=4)
    parser.add_argument("-t", "--token", dest="token", type=str, required=True)
    args = parser.parse_args()

    include = {}
    file_list = glob.glob("include/*") 
    for fl in file_list:    
        wf_name = os.path.splitext(os.path.basename(fl))[0]
        with open(fl, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                if not include.get(wf_name): include[wf_name] = set()
                include[wf_name].add(line.rstrip())


    #download data
    for wf in ['sanger', 'mutect2', 'mutect2-bqsr']:
        subfolder = 'evaluate/'+wf
        if not include.get(wf): continue 
        download(args.dump_path, 'snv', args.token, args.metadata_url, args.storage_url, include.get(wf), subfolder)
        download(args.dump_path, 'indel', args.token, args.metadata_url, args.storage_url, include.get(wf), subfolder)

    data_dir = 'data/evaluate'
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

    evaluate_list = []        
    for fn in glob.glob(os.path.join(data_dir, "*_report", "*-*", "*."+args.tool+'.stats.csv'), recursive=True):
        sub = os.path.dirname(fn).split('/')[-2].split('_')[0]
        studyId, donorId, sampleId, library, date_string, workflow, variant_type, evtype = os.path.basename(fn).split(".")[0:8]
        result_dict = {
            'studyId': studyId,
            'donorId': donorId,
            'sampleId': sampleId,
            'library_strategy': library,
            'workflow': sub,
            'evtype': evtype
        }
        with open(fn, 'r') as f:
            dict_reader = csv.DictReader(f, delimiter=",")
            for row in dict_reader:
                if row.get('type')=='records': continue
                for k in ['total.truth', 'total.query', 'tp', 'fp', 'fn', 'recall', 'recall_lower', 'recall_upper', 'recall2', 'precision', 'precision_lower', 'precision_upper', 'specificity', 'balaccuracy', 'f1_score']:
                    result_dict.update({k: row.get(k, None)})
                if row['total.query'] > 0:
                    result_dict['specificity'] = 1.0 - float(row['fp']) / float(row['total.query'])
                    result_dict['balaccuracy'] = (float(row['recall']) + result_dict['specificity']) / 2.0
                result_dict['f1_score'] = 2*float(row['recall'])*float(row['precision'])/(float(row['recall'])+float(row['precision']))
                
        evaluate_list.append(copy.deepcopy(result_dict)) 
    
    with open(os.path.join(report_dir, 'mutect2_evaluate_result.json'), 'w') as f:
        for s in evaluate_list:
            f.write(json.dumps(s)+'\n')

    # generate tsv file
    report(evaluate_list, os.path.join(report_dir, 'mutect2_evaluate_result.txt'))


if __name__ == "__main__":
    main()
