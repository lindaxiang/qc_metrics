#!/usr/bin/env python3
from argparse import ArgumentParser
import json
import functools
import os
import csv
import copy
from collections import OrderedDict

def report(donor, report_name):
    keys = donor[0].keys()
    with open(report_name, 'w') as output_file:
        dict_writer = csv.DictWriter(output_file, keys, delimiter="\t")
        dict_writer.writeheader()
        dict_writer.writerows(donor)

def udf(x, y):
    if not x:
        return
    if isinstance(x, dict):
        value = x.get(y)
    elif isinstance(x, list):
        value = [udf(a,y) for a in x]
    if value: 
        return value
    else:
        return

def process(gql_dump, suppress):
    with open(gql_dump, 'r') as fp:
        for fline in fp:
            suppress_dict = OrderedDict()
            analysis = json.loads(fline)
            if not analysis.get('inputForRuns'): continue
            # get donor info
            suppress_dict['studyId'] = analysis.get('studyId')
            suppress_dict['donorId'] = analysis['donors'][0]['donorId']
            suppress_dict['sampleId'] = analysis['donors'][0]['specimens'][0]['samples'][0]['sampleId']
            suppress_dict['tumourNormalDesignation'] = analysis['donors'][0]['specimens'][0]['tumourNormalDesignation']
            suppress_dict['run_input_analysisId'] = analysis.get('analysisId')
            suppress_dict['run_input_analysisType'] = analysis.get('analysisType')

            complete_count = 0
            complete_latest = 0
            inputForRuns = False
            for run in analysis['inputForRuns']:
                if run.get('state') == 'COMPLETE' and run.get('producedAnalyses'): 
                    complete_count = complete_count + 1
                    if int(run.get('completeTime')) > complete_latest:
                        complete_latest = int(run.get('completeTime'))
                    for out in run.get('producedAnalyses'):
                        if out.get('inputForRuns'): inputForRuns = True 
                    continue
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
                if not run.get('state') == 'COMPLETE': continue
                dup = True
                out_analysis = []
                for out in run.get('producedAnalyses'):
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
    parser.add_argument("-d", "--dump_path", dest="dump_path", type=str, default="gql", help="path to gql dump jsonl file")
    args = parser.parse_args()

    gql_dump_path = args.dump_path
    suppress = []
    suppress = process(os.path.join(gql_dump_path, 'sequencing_experiment.json'), suppress)
    suppress = process(os.path.join(gql_dump_path, 'sequencing_alignment.json'), suppress)
    #print(suppress_dict)

    report_dir = "report"
    with open(os.path.join(report_dir, 'suppress.json'), 'w') as f:
        for s in suppress:
            f.write(json.dumps(s)+'\n')

    # generate tsv file
    report(suppress, os.path.join(report_dir, 'suppress.txt'))

if __name__ == "__main__":
    main()