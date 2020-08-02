#!/usr/bin/env python3

import requests
import json
import os
import csv
from argparse import ArgumentParser
import glob

def get_analysis(file_dump):
    report_analysis = {}
    with open(file_dump, 'r') as fp:
        for fline in fp:
            analysis = json.loads(fline)            
            experimental_strategy = analysis['experiment']['experimental_strategy'] if analysis['experiment'].get('experimental_strategy') else analysis['experiment']['library_strategy']

            record = {
                'analysis_id': analysis['analysisId'],
                'analysis_type': analysis['analysisType']['name'],
                'analysis_state': analysis['analysisState'],
                'study_id': analysis['studyId'],
                'donor_id': analysis['samples'][0]['donor']['donorId'],
                'experimental_strategy': experimental_strategy,
                'sample_id': analysis['samples'][0]['sampleId'],
                'tumourNormalDesignation': analysis['samples'][0]['specimen']['tumourNormalDesignation'],
                'matched_normal_submitter_sample_id': analysis['samples'][0]['matchedNormalSubmitterSampleId'],
                'run_id': analysis['workflow']['run_id']
            }

            if analysis.get('workflow'):
                workflow_name = analysis['workflow']['workflow_name'] if analysis['workflow'].get('workflow_name') else analysis['workflow']['name']
                workflow_version = analysis['workflow']['workflow_version'] if analysis['workflow'].get('workflow_version') else analysis['workflow']['version']
                if workflow_name in ['Sanger WGS Variant Calling', 'Sanger WXS Variant Calling']:
                    if analysis['workflow']['inputs'][0].get('tumour_analysis_id'):
                        tumour_analysis_id = analysis['workflow']['inputs'][0]['tumour_analysis_id']
                        normal_analysis_id = analysis['workflow']['inputs'][1]['normal_analysis_id']
                    else:
                        tumour_analysis_id = analysis['workflow']['inputs'][1]['tumour_analysis_id']
                        normal_analysis_id = analysis['workflow']['inputs'][0]['normal_analysis_id']
                    workflow_info = {
                        'workflow_tumour_analysis_id': tumour_analysis_id,
                        'workflow_normal_analysis_id': normal_analysis_id,
                        'workflow_name': workflow_name,
                        'workflow_version': workflow_version
                    }
                    record_type = 'sanger_' + experimental_strategy.lower() + '_variant_calling'
                elif workflow_name in ['dna-seq-alignment', 'DNA Seq Alignment']:
                    workflow_info = {
                        'workflow_input_analysis_id': analysis['workflow']['inputs'][0]['input_analysis_id'],
                        'workflow_name': workflow_name,
                        'workflow_version': workflow_version
                    }
                    record_type = 'dna_seq_alignment_' + experimental_strategy.lower()
                else:
                    pass
            else:
                record_type = 'submitted_reads'

            if workflow_info: 
                record.update(workflow_info)

            if not report_analysis.get(record_type): report_analysis[record_type] = []
            report_analysis[record_type].append(record)
    return report_analysis

def process_sequencing_alignment(file_dump, exclude_alignment_analysis):
    job = {}
    sample_map = {}
    with open(file_dump, 'r') as fp:
        for fline in fp:
            analysis = json.loads(fline)
            if not analysis.get('analysisState') == 'PUBLISHED': continue
            if not analysis['analysisType']['name'] == 'sequencing_alignment': continue
            if analysis['workflow'].get('version') and analysis['workflow']['version'] == '1.0.4': continue
            if analysis['analysisId'] in exclude_alignment_analysis: continue
            if not analysis['samples'][0]['specimen']['tumourNormalDesignation'] == 'Tumour': continue
            sampleId = analysis['samples'][0]['sampleId']
            matchedNormal = analysis['samples'][0]['matchedNormalSubmitterSampleId']
            if not sample_map.get(analysis['studyId']+'::'+matchedNormal): 
                sample_map[analysis['studyId']+'::'+matchedNormal] = []
            sample_map[analysis['studyId']+'::'+matchedNormal].append(sampleId)
            if job.get(sampleId): print('sampleId: %s has duplicated analysis!' % sampleId)
            job[sampleId] = {
                'study_id': analysis['studyId'],
                'submitter_donor_id': analysis['samples'][0]['donor']['submitterDonorId'],
                'submitter_sample_id_n': None,
                'submitter_sample_id_t': analysis['samples'][0]['submitterSampleId'],
                'sequencing_stategy': analysis['experiment']['experimental_strategy'] if analysis['experiment'].get('experimental_strategy') else analysis['experiment']['library_strategy'],
                'donor_id': analysis['samples'][0]['donor']['donorId'],
                'sample_id_n': None,
                'sample_id_t': sampleId,
                'input_n_run_id': None,
                'input_t_run_id': analysis['workflow']['run_id'],
                'normal_aln_analysis_id': None,
                'tumour_aln_analysis_id': analysis['analysisId']
            }

    with open(file_dump, 'r') as fp:
        for fline in fp:
            analysis = json.loads(fline)
            if not analysis.get('analysisState') == 'PUBLISHED': continue
            if not analysis['analysisType']['name'] == 'sequencing_alignment': continue
            if analysis['workflow'].get('version') and analysis['workflow']['version'] == '1.0.4': continue
            if analysis['analysisId'] in exclude_alignment_analysis: continue
            if not analysis['samples'][0]['specimen']['tumourNormalDesignation'] == 'Normal': continue
            normal_sample_id = analysis['studyId']+'::'+analysis['samples'][0]['submitterSampleId']
            if not normal_sample_id in sample_map: continue
            for sa in sample_map[normal_sample_id]:          
                job[sa]['submitter_sample_id_n'] = analysis['samples'][0]['submitterSampleId']
                job[sa]['sample_id_n'] = analysis['samples'][0]['sampleId']
                job[sa]['input_n_run_id'] = analysis['workflow']['run_id']
                job[sa]['normal_aln_analysis_id'] = analysis['analysisId']
            
    return job

def report(donor, report_name):
    keys = donor[0].keys()
    with open(report_name, 'w') as output_file:
        dict_writer = csv.DictWriter(output_file, keys, delimiter="\t")
        dict_writer.writeheader()
        dict_writer.writerows(donor)

def main():
    parser = ArgumentParser()
    parser.add_argument("-d", "--dump_path", dest="dump_path", type=str, required=True, help="path to song dump jsonl file")
    parser.add_argument("-x", "--exclude", dest="exclude", action='store_true')
    
    args = parser.parse_args()

    song_dump = args.dump_path

    annotation = {
        'exclude_list': set(),
        'schedule_list': set()
    }
    file_list = glob.glob("exclude_*.*") if args.exclude else []
    for fl in file_list:
        with open(fl, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                annotation['exclude_list'].add(line.rstrip())
    
    file_list = glob.glob("schedule_*.*") if args.exclude else []
    for fl in file_list:
        with open(fl, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                annotation['schedule_list'].add(line.rstrip())
    

    report_analysis = get_analysis(song_dump)
    for k, v in report_analysis.items():
        report(v, 'report_'+k+'.tsv')

    job = process_sequencing_alignment(song_dump, annotation['exclude_list'])
    report_job = []
    for d, v in job.items():
        if not v.get('normal_aln_analysis_id') or not v.get('tumour_aln_analysis_id'): continue
        if v.get('normal_aln_analysis_id') in annotation['schedule_list'] and v.get('tumour_aln_analysis_id') in annotation['schedule_list']: continue
        report_job.append(v)
    report(report_job, 'sanger_jobs_to_stage.tsv')


if __name__ == "__main__":
    main()
