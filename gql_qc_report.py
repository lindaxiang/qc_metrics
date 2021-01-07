#!/usr/bin/env python3
from argparse import ArgumentParser
import json
import functools
import os
import csv
import copy
from collections import OrderedDict
from utils import report, udf

wf_repo = {
    'alignment': 'https://github.com/icgc-argo/dna-seq-processing-wfs.git',
    'sanger-wgs': 'https://github.com/icgc-argo/sanger-wgs-variant-calling.git',
    'sanger-wxs': 'https://github.com/icgc-argo/sanger-wxs-variant-calling.git',
    'gatk-mutect2': 'https://github.com/icgc-argo/gatk-mutect2-variant-calling.git'
}

def process(gql_dump, analysisType, suppress):
    with open(gql_dump, 'r') as fp:
        for fline in fp:
            suppress_dict = OrderedDict()
            analysis = json.loads(fline)
            if not analysis.get('analysisType') == analysisType: continue
            if not analysis.get('inputForRuns'): continue
            # get donor info
            suppress_dict['studyId'] = analysis.get('studyId')
            suppress_dict['donorId'] = analysis['donors'][0]['donorId']
            suppress_dict['sampleId'] = analysis['donors'][0]['specimens'][0]['samples'][0]['sampleId']
            suppress_dict['tumourNormalDesignation'] = analysis['donors'][0]['specimens'][0]['tumourNormalDesignation']
            suppress_dict['experimental_strategy'] = analysis['experiment']['experimental_strategy']
            suppress_dict['run_input_analysisId'] = analysis.get('analysisId')
            suppress_dict['run_input_analysisType'] = analysis.get('analysisType')

            for wf in wf_repo:
                complete_count = 0
                complete_latest = 0
                inputForRuns = False
                for run in analysis['inputForRuns']:
                    if not wf_repo.get(wf) == run.get('repository'): continue
                    if run.get('state') == 'COMPLETE' and run.get('producedAnalyses'): 
                        complete_count = complete_count + 1
                        if int(run.get('completeTime')) > complete_latest:
                            complete_latest = int(run.get('completeTime'))
                        for out in run.get('producedAnalyses'):
                            if out.get('inputForRuns'): inputForRuns = True 
                        continue

                    if run.get('state') == 'RUNNING': continue
                    
                    # failed runs
                    failed_run = functools.reduce(udf, ['producedAnalyses', 'analysisId'], run)
                    if failed_run: 
                        suppress_dict['runId'] = run.get('runId')
                        suppress_dict['run_output_analysis_to_suppress'] = failed_run
                        suppress_dict['reason'] = 'failed_run_output_analysis'
                        suppress.append(copy.deepcopy(suppress_dict))
                
                if complete_count < 2: continue
                if suppress_dict['tumourNormalDesignation'] == 'Normal' and suppress_dict['run_input_analysisType'] == 'sequencing_alignment': continue
                # duplicate runs
                for run in analysis['inputForRuns']:
                    if not wf_repo.get(wf) == run.get('repository'): continue
                    if not run.get('state') == 'COMPLETE': continue
                    dup = True
                    out_analysis = []
                    for out in run.get('producedAnalyses'):
                        if not out.get('analysisState') == 'PUBLISHED': continue
                        out_analysis.append(out.get('analysisId'))
                        # keep the run which already being used as input for other workflows
                        if out.get('analysisType') == 'sequencing_alignment' and out.get('inputForRuns'): dup = False
                    # keep the latest run if no run results having been used
                    if not inputForRuns and int(run.get('completeTime')) == complete_latest: dup = False
                    # if no output analysis are found (already suppressed)
                    if not out_analysis: dup = False

                    if dup:
                        suppress_dict['runId'] = run.get('runId')
                        suppress_dict['run_output_analysis_to_suppress'] = out_analysis
                        suppress_dict['reason'] = 'duplicated_run_output_analysis'
                        suppress.append(copy.deepcopy(suppress_dict))      

    return suppress



def main():
    parser = ArgumentParser()
    parser.add_argument("-d", "--gql_dump", dest="gql_dump", type=str, default="data/rdpc-gql-analysis.jsonl", help="gql dump jsonl file")
    args = parser.parse_args()

    gql_dump = args.gql_dump
    suppress = []
    for analysisType in ['sequencing_experiment', 'sequencing_alignment']:
        suppress = process(gql_dump, analysisType, suppress)
    

    report_dir = "report"
    if suppress:
        with open(os.path.join(report_dir, 'suppress.json'), 'w') as f:
            for s in suppress:
                f.write(json.dumps(s)+'\n')

        # generate tsv file
        report(suppress, os.path.join(report_dir, 'suppress.txt'))
    else:
        print("No items to suppress!")

if __name__ == "__main__":
    main()