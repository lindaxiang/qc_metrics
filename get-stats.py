#!/usr/bin/env python3

import json
import os
import csv
import glob
from argparse import ArgumentParser
import sys
import subprocess
from collections import OrderedDict
import functools
import re
import pandas as pd
import numpy as np
from utils import report, run_cmd, get_dict_value, download, annot_vcf, get_extra_metrics

pd.options.mode.chained_assignment = None  # default='warn'

variant_calling_stats_fields = {
    'donor_id': 'donor_id',
    'study_id': 'study_id',
    # 'gender': 'gender',
    'experimental_strategy': 'experimental_strategy',
    # 'geno_infer_gender': 'geno_infer_gender',
    'normal_aligned': 'flags.normal_aligned',
    'tumour_aligned': 'flags.tumour_aligned',
    'sanger_called': 'flags.sanger_called',
    'mutect2_called': 'flags.mutect2_called',
    'is_pcawg': 'flags.is_pcawg',
    'normal_sample_id': 'normal.sample_id',
    #'normal_submitter_sample_id': 'normal.submitterSampleId',
    # 'normal_file_size_gb': 'normal.alignment.file_size',
    'normal_error_rate': 'normal.alignment.error_rate',
    'normal_duplicate_rate': 'normal.alignment.duplicate_rate',
    # 'normal_pairs_on_different_chromosomes': 'normal.alignment.pairs_on_different_chromosomes',
    'normal_pairs_on_different_chromosomes_rate': 'normal.alignment.pairs_on_different_chromosomes_rate',
    'normal_oxoQ_score': 'normal.alignment.oxoQ_score',
    'normal_insert_size_mean': 'normal.alignment.average_insert_size',
    'normal_avg_depth': 'normal.sanger.contamination.avg_depth',
    'normal_estimated_coverage': 'normal.alignment.estimated_coverage',
    'normal_sanger_contamination': 'normal.sanger.contamination.contamination',
    'normal_mutect2_contamination': 'normal.mutect2.contamination.contamination',
    'normal_properly_paired_reads': 'normal.alignment.properly_paired_reads',
    # 'normal_total_reads': 'normal.alignment.total_reads',
    'tumour_sample_id': 'tumour.sample_id',
    #'tumour_submitter_sample_id': 'tumour.submitterSampleId',
    # 'tumour_file_size_gb': 'tumour.alignment.file_size',
    'tumour_error_rate': 'tumour.alignment.error_rate',
    'tumour_duplicate_rate': 'tumour.alignment.duplicate_rate',
    # 'tumour_pairs_on_different_chromosomes': 'tumour.alignment.pairs_on_different_chromosomes',
    'tumour_pairs_on_different_chromosomes_rate': 'tumour.alignment.pairs_on_different_chromosomes_rate',
    'tumour_oxoQ_score': 'tumour.alignment.oxoQ_score',
    'tumour_insert_size_mean': 'tumour.alignment.average_insert_size',
    'tumour_avg_depth': 'tumour.sanger.contamination.avg_depth',
    'tumour_estimated_coverage': 'tumour.alignment.estimated_coverage',
    'tumour_sanger_contamination': 'tumour.sanger.contamination.contamination',
    'tumour_mutect2_contamination': 'tumour.mutect2.contamination.contamination',
    'tumour_properly_paired_reads': 'tumour.alignment.properly_paired_reads',
    # 'tumour_total_reads': 'tumour.alignment.total_reads',
    # 'ascat_normal_contamination': 'tumour.sanger.ascat_metrics.NormalContamination',
    'ascat_ploidy': 'tumour.sanger.ascat_metrics.Ploidy',
    # 'ascat_goodnessOfFit': 'tumour.sanger.ascat_metrics.goodnessOfFit',
    # 'ascat_psi': 'tumour.sanger.ascat_metrics.psi',
    'ascat_purity': 'tumour.sanger.ascat_metrics.rho',
    'mutect2_callable': 'tumour.mutect2.callable'

    # 'cgpPindel_cpu_hours': 'tumour.sanger.timing.cgpPindel.cpu_hours',
    # 'cgpPindel_max_memory_usage_per_core': 'tumour.sanger.timing.cgpPindel.maximum_memory_usage_per_core', 
    # 'CaVEMan_cpu_hours': 'tumour.sanger.timing.CaVEMan.cpu_hours',
    # 'CaVEMan_max_memory_usage_per_core': 'tumour.sanger.timing.CaVEMan.maximum_memory_usage_per_core',
    # 'BRASS_cpu_hours': 'tumour.sanger.timing.BRASS.cpu_hours',
    # 'BRASS_max_memory_usage_per_core': 'tumour.sanger.timing.BRASS.maximum_memory_usage_per_core',
    # 'ascat_cpu_hours': 'tumour.sanger.timing.ascat.cpu_hours',
    # 'ascat_max_memory_usage_per_core': 'tumour.sanger.timing.ascat.maximum_memory_usage_per_core',
    # 'snv_somatic_pass_total': 'gnomad_overlap.snv.somatic_pass_total',
    # 'indel_somatic_pass_total': 'gnomad_overlap.indel.somatic_pass_total',
    # 'gnomad_overlap_snv_t_0': 'gnomad_overlap.snv.t_0',
    # 'gnomad_overlap_snv_t_0.001': 'gnomad_overlap.snv.t_0_001',
    # 'gnomad_overlap_snv_t_0.01': 'gnomad_overlap.snv.t_0_01',
    # 'gnomad_overlap_snv_t_0.1': 'gnomad_overlap.snv.t_0_1',
    # 'gnomad_overlap_indel_t_0': 'gnomad_overlap.indel.t_0',
    # 'gnomad_overlap_indel_t_0.001': 'gnomad_overlap.indel.t_0_001',
    # 'gnomad_overlap_indel_t_0.01': 'gnomad_overlap.indel.t_0_01',
    # 'gnomad_overlap_indel_t_0.1': 'gnomad_overlap.indel.t_0_1',
    # 'gnomad_overlap_snv_t_0_count': 'gnomad_overlap.snv.t_0_count',
    # 'gnomad_overlap_snv_t_0.001_count': 'gnomad_overlap.snv.t_0_001_count',
    # 'gnomad_overlap_snv_t_0.01_count': 'gnomad_overlap.snv.t_0_01_count',
    # 'gnomad_overlap_snv_t_0.1_count': 'gnomad_overlap.snv.t_0_1_count',
    # 'gnomad_overlap_indel_t_0_count': 'gnomad_overlap.indel.t_0_count',
    # 'gnomad_overlap_indel_t_0.001_count': 'gnomad_overlap.indel.t_0_001_count',
    # 'gnomad_overlap_indel_t_0.01_count': 'gnomad_overlap.indel.t_0_01_count',
    # 'gnomad_overlap_indel_t_0.1_count': 'gnomad_overlap.indel.t_0_1_count'
}

pcawg_qc_fields = {
    'donor_id': 'donor_id',
    'study_id': 'study_id',
    'gender': 'gender',
    'experimental_strategy': 'experimental_strategy',
    'sanger_called': 'flags.sanger_called',
    'mutect2_called': 'flags.mutect2_called',
    'is_pcawg': 'flags.is_pcawg',
    'normal_sample_id': 'normal.sample_id',
    'tumour_sample_id': 'tumour.sample_id',
    'ARGO_alignment_normal_insert_size_mean': 'normal.alignment.average_insert_size',
    'ARGO_alignment_normal_read_length_mean': 'normal.alignment.average_length',
    'ARGO_alignment_tumour_insert_size_mean': 'tumour.alignment.average_insert_size',
    'ARGO_alignment_tumour_read_length_mean': 'tumour.alignment.average_length',
    'ARGO_alignment_normal_insert_size_sd': 'normal.alignment.insert_size_sd',
    'ARGO_alignment_tumour_insert_size_sd': 'tumour.alignment.insert_size_sd',
    'ARGO_alignment_normal_pairs_on_different_chromosomes_rate': 'normal.alignment.pairs_on_different_chromosomes_rate',
    'ARGO_alignment_tumour_pairs_on_different_chromosomes_rate': 'tumour.alignment.pairs_on_different_chromosomes_rate',
    'ARGO_sanger_normal_contamination': 'normal.sanger.contamination.contamination',
    'PCAWG_sanger_normal_contamination': 'pcawg.normal.sanger.contamination.contamination',
    'ARGO_mutect2_normal_contamination': 'normal.mutect2.contamination.contamination',
    'ARGO_sanger_tumour_contamination': 'tumour.sanger.contamination.contamination',
    'PCAWG_sanger_tumour_contamination': 'pcawg.tumour.sanger.contamination.contamination',
    'ARGO_mutect2_tumour_contamination': 'tumour.mutect2.contamination.contamination',
    'PCAWG_broad_tumour_contamination': 'pcawg.tumour.broad.contamination',
    'ARGO_sanger_ascat_normal_contamination': 'tumour.sanger.ascat_metrics.NormalContamination',
    'PCAWG_sanger_ascat_normal_contamination': 'pcawg.tumour.sanger.ascat_metrics.NormalContamination',
    'ARGO_sanger_ascat_ploidy': 'tumour.sanger.ascat_metrics.Ploidy',
    'PCAWG_sanger_ascat_ploidy': 'pcawg.tumour.sanger.ascat_metrics.Ploidy',
    'ARGO_sanger_normal_avg_depth': 'normal.sanger.contamination.avg_depth',
    'PCAWG_sanger_normal_avg_depth': 'pcawg.normal.sanger.contamination.avg_depth',
    'ARGO_sanger_tumour_avg_depth': 'tumour.sanger.contamination.avg_depth',
    'PCAWG_sanger_tumour_avg_depth': 'pcawg.tumour.sanger.contamination.avg_depth',
    'ARGO_sanger_ascat_goodnessOfFit': 'tumour.sanger.ascat_metrics.goodnessOfFit',
    'PCAWG_sanger_ascat_goodnessOfFit': 'pcawg.tumour.sanger.ascat_metrics.goodnessOfFit',
    'ARGO_sanger_ascat_psi': 'tumour.sanger.ascat_metrics.psi',
    'PCAWG_sanger_ascat_psi': 'pcawg.tumour.sanger.ascat_metrics.psi',
    'ARGO_sanger_ascat_purity': 'tumour.sanger.ascat_metrics.rho',
    'PCAWG_sanger_ascat_purity': 'pcawg.tumour.sanger.ascat_metrics.rho' 
}


total_size = {
    'wgs': 3200000000,
    'wxs': 99810084
}


def add_pcawg_info(variant_calling_stats, pcawg_sample_sheet, pcawg_sanger_qc, pcawg_broad_qc):
    pcawg_sample_info = {}
    with open(pcawg_sample_sheet, 'r') as fp:
        reader = csv.DictReader(fp, delimiter='\t')
        for row in reader:
            if not row.get('library_strategy') == "WGS": continue
            if pcawg_sample_info.get(row.get('aliquot_id')): continue
            pcawg_sample_info[row.get('aliquot_id')] = row.get('icgc_sample_id')

    pcawg_sample_qc = {}
    with open(pcawg_sanger_qc, 'r') as fp:
        for fline in fp:
            sanger_qc = json.loads(fline)
            for key, value in sanger_qc.items():
                if not key in pcawg_sample_info: continue
                if pcawg_sample_info[key] in pcawg_sample_qc: continue
                pcawg_sample_qc['WGS_'+pcawg_sample_info[key]] = {}
                pcawg_sample_qc['WGS_'+pcawg_sample_info[key]]['sanger'] = {}
                value['contamination'][key].pop('by_readgroup')
                pcawg_sample_qc['WGS_'+pcawg_sample_info[key]]['sanger'].update({'contamination': value['contamination'][key]})
                if 'cnv' in value: 
                    pcawg_sample_qc['WGS_'+pcawg_sample_info[key]]['sanger'].update({'ascat_metrics': value.get('cnv')})
    
    with open(pcawg_broad_qc, 'r') as fp:
        reader = csv.DictReader(fp, delimiter='\t')
        for row in reader:
            if not row.get('aliquot_GUUID') in pcawg_sample_info: continue
            sample_id = pcawg_sample_info[row.get('aliquot_GUUID')]
            if not 'WGS_'+sample_id in pcawg_sample_qc: 
                pcawg_sample_qc['WGS_'+sample_id] = {}
            pcawg_sample_qc['WGS_'+sample_id]['broad'] = {
                'contamination': float(row.get('contamination_percentage_whole_genome_no_array_value'))/100 if row.get('contamination_percentage_whole_genome_no_array_value') else None,
                'oxoQ': row.get('picard_oxoQ') if row.get('picard_oxoQ') else None,
                'callable': row.get('somatic_mutation_covered_bases_wgs') if row.get('somatic_mutation_covered_bases_wgs') else None
            }

    for key, value in variant_calling_stats.items():
        value['flags']['is_pcawg'] = False
        if not key in pcawg_sample_qc: continue
        value['flags']['is_pcawg'] = True
        value['pcawg'] = {}
        value['pcawg']['tumour'] = pcawg_sample_qc.get(key)
        if not 'sample_id' in value['normal']: continue
        if 'WGS_'+value['normal']['sample_id'] in pcawg_sample_qc:
            value['pcawg']['normal'] = pcawg_sample_qc.get('WGS_'+value['normal']['sample_id'])
        
    return variant_calling_stats



def process_qc_metrics(song_dump, variant_calling_stats):
    sample_map = {}
    with open(song_dump, 'r') as fp:
        for fline in fp:
            analysis = json.loads(fline)
            if not analysis.get('analysisState') == 'PUBLISHED': continue
            if not analysis['samples'][0]['specimen']['tumourNormalDesignation'] == 'Tumour': continue
            studyId = analysis['studyId']
            sampleId = analysis['samples'][0]['sampleId']
            submitterSampleId = analysis['samples'][0]['submitterSampleId']
            matchedNormal = analysis['samples'][0]['matchedNormalSubmitterSampleId']
            experimental_strategy = analysis['experiment']['experimental_strategy'] if analysis['experiment'].get('experimental_strategy') else analysis['experiment']['library_strategy']
            normal_sample_id = '_'.join([studyId, experimental_strategy, matchedNormal])
            if not sample_map.get(normal_sample_id): 
                sample_map[normal_sample_id] = []
            sample_map[normal_sample_id].append(experimental_strategy+"_"+sampleId)

            donorId = analysis['samples'][0]['donor']['donorId']
            gender = analysis['samples'][0]['donor']['gender']
            
            unique_sampleId = experimental_strategy+"_"+sampleId
            if not variant_calling_stats.get(unique_sampleId): variant_calling_stats[unique_sampleId] = {
                'study_id': studyId,
                'donor_id': donorId,
                'gender': gender,
                'experimental_strategy': experimental_strategy,
                'flags': {
                    'normal_aligned': False,
                    'tumour_aligned': False,
                    'sanger_called': False,
                    'mutect2_called': False
                },
                'normal': {
                    'alignment': {},
                    'sanger': {
                        'contamination': {}
                    },
                    'mutect2': {
                        'contamination': {}
                    }
                },
                'tumour': {
                    'sample_id': sampleId,
                    'submitterSampleId': submitterSampleId,
                    'alignment': {},
                    'sanger': {
                        'contamination': {},
                        'ascat_metrics': {},
                        'genotype_inference': {}
                    },
                    'mutect2': {
                        'contamination': {}
                    }
                }
            }

            if analysis['analysisType']['name'] == 'variant_calling': 
                if analysis['workflow']['workflow_short_name'] in ['sanger-wgs', 'sanger-wxs']:
                    variant_calling_stats[unique_sampleId]['flags']['sanger_called'] = True
                if analysis['workflow']['workflow_short_name'] == 'gatk-mutect2':
                    variant_calling_stats[unique_sampleId]['flags']['mutect2_called'] = True
            elif analysis['analysisType']['name'] == 'sequencing_alignment':
                variant_calling_stats[unique_sampleId]['flags']['tumour_aligned'] = True 
            elif not analysis['analysisType']['name'] == 'qc_metrics': 
                continue

            for fl in analysis['files']:
                if fl['dataType'] == 'Cross Sample Contamination':
                    if fl['info']['metrics']['sample_id'] == sampleId:
                        if 'sanger' in fl['fileName']:
                            variant_calling_stats[unique_sampleId]['tumour']['sanger']['contamination'].update(fl['info']['metrics'])
                        elif 'gatk-mutect2' in fl['fileName']:
                            variant_calling_stats[unique_sampleId]['tumour']['mutect2']['contamination'].update(fl['info']['metrics'])
                        else:
                            pass
                    else:
                        if 'sanger' in fl['fileName']:
                            variant_calling_stats[unique_sampleId]['normal']['sanger']['contamination'].update(fl['info']['metrics'])
                        elif 'gatk-mutect2' in fl['fileName']:
                            variant_calling_stats[unique_sampleId]['normal']['mutect2']['contamination'].update(fl['info']['metrics'])
                        else:
                            pass
                elif fl['dataType'] == 'Ploidy and Purity Estimation':
                    variant_calling_stats[unique_sampleId]['tumour']['sanger']['ascat_metrics'].update(fl['info']['metrics'])
                elif fl['dataType'] == 'Genotyping Inferred Gender':
                    variant_calling_stats[unique_sampleId]['tumour']['sanger']['genotype_inference'].update(fl['info']['metrics']['tumours'][0]['gender'])
                elif fl['dataType'] == 'Alignment QC' and 'qc_metrics' in fl['fileName']:
                    metrics = {}
                    for fn in ['error_rate', 'properly_paired_reads', 'total_reads', 'average_insert_size', 'average_length', 'pairs_on_different_chromosomes']:
                        metrics.update({fn: fl['info']['metrics'][fn]})
                    metrics.update({
                        'duplicate_rate': round(fl['info']['metrics']['duplicated_bases']/(fl['info']['metrics']['total_reads']*fl['info']['metrics']['average_length']), 3)
                        })
                    metrics.update({
                        'pairs_on_different_chromosomes_rate': round(fl['info']['metrics']['pairs_on_different_chromosomes']*2/(fl['info']['metrics']['paired_reads']), 3)
                    })
                    metrics.update({
                        'estimated_coverage': round(fl['info']['metrics']['mapped_bases_cigar']/total_size.get(experimental_strategy.lower()), 3)
                    })
                    fname = os.path.join("data", 'qc_metrics', analysis['studyId'], fl['fileName'])
                    extra_metrics = ['insert_size_sd']
                    metrics = get_extra_metrics(fname, extra_metrics, metrics)
                    variant_calling_stats[unique_sampleId]['tumour']['alignment'].update(metrics)
                elif fl['dataType'] == 'OxoG Metrics':
                    variant_calling_stats[unique_sampleId]['tumour']['alignment'].update({'oxoQ_score': fl['info']['metrics']['oxoQ_score']})
                    
                elif fl['dataType'] == 'Mutect2 Callable Stats':
                    variant_calling_stats[unique_sampleId]['tumour']['mutect2'].update(fl['info']['metrics'])

                elif fl['dataType'] == 'Aligned Reads':
                    variant_calling_stats[unique_sampleId]['tumour']['alignment'].update({"file_size": round(fl['fileSize']/(1024*1024*1024), 3)})

                else:
                    continue

    with open(song_dump, 'r') as fp:
        for fline in fp:
            analysis = json.loads(fline)
            if not analysis.get('analysisState') == 'PUBLISHED': continue
            if not analysis['analysisType']['name'] in ['qc_metrics', 'sequencing_alignment']: continue
            if not analysis['samples'][0]['specimen']['tumourNormalDesignation'] == 'Normal': continue

            experimental_strategy = analysis['experiment']['experimental_strategy'] if analysis['experiment'].get('experimental_strategy') else analysis['experiment']['library_strategy']           
            studyId = analysis['studyId']
            submitSampleId = analysis['samples'][0]['submitterSampleId']
            normal_sample_id = '_'.join([studyId, experimental_strategy, submitSampleId])
            
            if not normal_sample_id in sample_map: continue
            
            for fl in analysis['files']:
                if fl['dataType'] == 'Alignment QC' and 'qc_metrics' in fl['fileName']:
                    metrics = {}
                    for fn in ['error_rate', 'properly_paired_reads', 'total_reads', 'average_insert_size', 'average_length', 'pairs_on_different_chromosomes']:
                        metrics.update({fn: fl['info']['metrics'][fn]})
                    metrics.update({
                        'duplicate_rate': round(fl['info']['metrics']['duplicated_bases']/(fl['info']['metrics']['total_reads']*fl['info']['metrics']['average_length']), 3)
                        })
                    if fl['info']['metrics']['paired_reads']>0:
                        metrics.update({
                            'pairs_on_different_chromosomes_rate': round(fl['info']['metrics']['pairs_on_different_chromosomes']*2/(fl['info']['metrics']['paired_reads']), 3)
                        })
                    metrics.update({
                        'estimated_coverage': round(fl['info']['metrics']['mapped_bases_cigar']/total_size.get(experimental_strategy.lower()), 3)
                    })
                    fname = os.path.join("data", 'qc_metrics', analysis['studyId'], fl['fileName'])
                    extra_metrics = ['insert_size_sd']
                    metrics = get_extra_metrics(fname, extra_metrics, metrics)

                    for sa in sample_map[normal_sample_id]:
                        variant_calling_stats[sa]['normal']['sample_id'] = analysis['samples'][0]['sampleId']
                        variant_calling_stats[sa]['normal']['submitterSampleId'] = analysis['samples'][0]['submitterSampleId']  
                        variant_calling_stats[sa]['normal']['alignment'].update(metrics)
                        variant_calling_stats[sa]['flags']['normal_aligned'] = True 
                elif fl['dataType'] == 'OxoG Metrics':
                    for sa in sample_map[normal_sample_id]:  
                        variant_calling_stats[sa]['normal']['alignment'].update({'oxoQ_score': fl['info']['metrics'].get('oxoQ_score', None)})                 
                elif fl['dataType'] == 'Aligned Reads':
                    for sa in sample_map[normal_sample_id]:  
                        variant_calling_stats[sa]['normal']['alignment'].update({"file_size": round(fl['fileSize']/(1024*1024*1024), 3)})                    
                else:
                    continue                   

    return variant_calling_stats


def get_gnomad_overlap(vcf, af_threshold, annotated):
    
    basename = re.sub(r'.vcf.gz$', '', vcf)
    gnomad_file = f'{basename}.gnomad_af.txt'

    if not os.path.basename(gnomad_file) in annotated:
        bcftools = f"bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%gnomad_af\t%gnomad_filter\n' {vcf} > {basename}.gnomad_af.txt" 
        run_cmd(bcftools)

    df = pd.read_table(gnomad_file, sep='\t', \
         names=["CHROM", "POS", "REF", "ALT", "FILTER", "gnomad_af", "gnomad_filter"], \
         dtype={"CHROM": str, "POS": int, "REF": str, "ALT": str, "FILTER": str, "gnomad_af": str, "gnomad_filter": str}, \
         na_values=".")
    
    df_somatic_pass = df.loc[df['FILTER']=="PASS", :]
    somatic_pass_total = df_somatic_pass.shape[0]
    gnomad_af = {
        "somatic_pass_total": somatic_pass_total
    }
    # convert to numeric dtype for all pass variants
    df_somatic_pass['gnomad_af'] = pd.to_numeric(df_somatic_pass['gnomad_af'])
    for t in af_threshold:
        gnomad_af['t_'+str(t).replace('.', '_')+'_count'] = df_somatic_pass.loc[df_somatic_pass['gnomad_af'] > t, :].shape[0]
        if somatic_pass_total > 0:
            gnomad_af['t_'+str(t).replace('.', '_')] = round(gnomad_af['t_'+str(t).replace('.', '_')+'_count']/somatic_pass_total, 3)
        else:
            gnomad_af['t_'+str(t).replace('.', '_')] = None

    return gnomad_af


def process_annot_vcf(variant_calling_stats, af_threshold):
    annot_dir = os.path.join("data", "annot_vcf")
    annotated = []
    for fn in glob.glob(os.path.join(annot_dir, "*-*", "*.*"), recursive=True):
        annotated.append(os.path.basename(fn))

    for vcf in glob.glob(os.path.join(annot_dir, "*-*", '*.vcf.gz'), recursive=True):
        tumour_sample_id = os.path.basename(vcf).split('.')[2]
        data_type = os.path.basename(vcf).split('.')[7] 
        
        if not 'gnomad_overlap' in variant_calling_stats[tumour_sample_id]:
            variant_calling_stats[tumour_sample_id]['gnomad_overlap'] = {}
        variant_calling_stats[tumour_sample_id]['gnomad_overlap'][data_type] = {}
        gnomad = get_gnomad_overlap(vcf, af_threshold, annotated)
        variant_calling_stats[tumour_sample_id]['gnomad_overlap'][data_type].update(gnomad)
    
    return variant_calling_stats


def get_timing(fname, experimental_strategy):
    sanger_timing = {}
    start = False
    with open(fname, 'r') as f:
        sanger_timing[os.path.basename(fname).split('.')[-1]] = {} 
        if experimental_strategy == 'wgs':                            
            for row in f:
                if row.strip().startswith('Command being timed'): 
                    start = True
                    continue
                if not start == True: 
                    continue
                cols = row.strip().split(': ')
                if not cols[0] in ['Percent of CPU this job got', 
                                    'User time (seconds)', 
                                    'System time (seconds)',
                                'Elapsed (wall clock) time (h:mm:ss or m:ss)', 
                                'Maximum resident set size (kbytes)',
                                'File system inputs', 'File system outputs']: 
                    continue
                
                if cols[0] == 'Percent of CPU this job got':
                    value = round(float(cols[1][:-1])/100.0, 1)
                    fieldname = 'cpu_usage'
                elif cols[0] == 'Elapsed (wall clock) time (h:mm:ss or m:ss)':
                    hms = cols[1].split(':')
                    if len(hms) == 3:
                        value = int(hms[0])+int(hms[1])/60+float(hms[2])/3600
                    else:
                        value = int(hms[0])/60+float(hms[1])/3600
                    value = round(value, 3)
                    fieldname = 'real_running_hours'
                elif cols[0] == 'User time (seconds)':
                    cpu_time = float(cols[1])
                    continue
                elif cols[0] == 'System time (seconds)':
                    cpu_time = cpu_time + float(cols[1])
                    value = round(cpu_time/3600, 3)
                    fieldname = 'cpu_hours'
                elif cols[0] == 'Maximum resident set size (kbytes)':
                    fieldname = 'maximum_memory_usage_per_core'
                    value = round(int(cols[1])/(1024*1024), 3)
                elif cols[0] == 'File system inputs':
                    fieldname = 'file_inputs_size_gb'
                    value= round(int(cols[1])*512/(1024*1024*1024), 3)
                elif cols[0] == 'File system outputs':
                    fieldname = 'file_outputs_size_gb'
                    value= round(int(cols[1])*512/(1024*1024*1024), 3)
                    #total_output_size = total_output_size + value
                else:
                    continue
                sanger_timing[os.path.basename(fname).split('.')[-1]].update({fieldname: value})

        elif experimental_strategy == 'wxs':
            for row in f:
                if row.strip().startswith('command'): 
                    start = True
                    continue
                if not start == True: 
                    continue

                cols = row.strip().split(':')

                if not cols[0] in ['real', 'user', 'sys', 'pctCpu', 'max']: 
                    continue
                
                if cols[0] == 'pctCpu':
                    value = round(float(cols[1][:-1])/100.0, 1)
                    fieldname = 'cpu_usage'
                elif cols[0] == 'real':
                    value = round(float(cols[1])/3600, 3)
                    fieldname = 'real_running_hours'
                elif cols[0] == 'user':
                    cpu_time = float(cols[1])
                    continue
                elif cols[0] == 'sys':
                    cpu_time = cpu_time + float(cols[1])
                    value = round(cpu_time/3600, 3)
                    fieldname = 'cpu_hours'
                elif cols[0] == 'max':
                    fieldname = 'maximum_memory_usage_per_core'
                    value = round(int(cols[1].strip('k'))/(1024*1024), 3)
                else:
                    continue
                sanger_timing[os.path.basename(fname).split('.')[-1]].update({fieldname: value})
        else:
            pass
    return sanger_timing


def process_timing_supplement(variant_calling_stats):
    supplement_dir = os.path.join("data", 'variant_calling_supplement')
    for tgz in glob.glob(os.path.join(supplement_dir, '*-*', '*.timings-supplement.tgz'), recursive=True):
        tumour_sample_id = os.path.basename(tgz).split('.')[2]
        experimental_strategy = os.path.basename(tgz).split('.')[3]

        if os.path.isdir(os.path.join(supplement_dir, 'unzip')): 
            cmd = 'rm -rf %s && mkdir %s && tar -C %s -xzf %s' % (os.path.join(supplement_dir, 'unzip'), os.path.join(supplement_dir, 'unzip'), os.path.join(supplement_dir, 'unzip'), tgz)
        else:
            cmd = 'mkdir %s && tar -C %s -xzf %s' % (os.path.join(supplement_dir, 'unzip'), os.path.join(supplement_dir, 'unzip'), tgz)
        run_cmd(cmd)

        variant_calling_stats[tumour_sample_id]['tumour']['sanger']['timing'] = {}
        # total_output_size = 0 
        for fname in glob.glob(os.path.join(supplement_dir, 'unzip', 'timings', "*.time*")):
            if 'annot' in fname: continue
            sanger_timing = get_timing(fname, experimental_strategy)
            variant_calling_stats[tumour_sample_id]['tumour']['sanger']['timing'].update(sanger_timing)
        # variant_calling_stats[tumour_sample_id]['timing'].update({"total_output_size_gb": total_output_size})           
    return variant_calling_stats    


def main():
    parser = ArgumentParser()
    parser.add_argument("-d", "--dump_path", dest="dump_path", type=str, default="data/rdpc-song.jsonl", help="path to song dump jsonl file")
    parser.add_argument("-a", "--pcawg_sample_sheet", dest="pcawg_sample_sheet", type=str, default="data/pcawg_sample_sheet.tsv", help="path to pcwag sample sheet file")
    parser.add_argument("-q", "--pcawg_sanger_qc", dest="pcawg_sanger_qc", type=str, default="data/pcawg_sanger_qc_metrics.jsonl", help="path to pcawg sanger qc file")
    parser.add_argument("-b", "--pcawg_broad_qc", dest="pcawg_broad_qc", type=str, default="data/pcawg_broad_qc_metrics.tsv", help="path to pcawg broad qc file")
    parser.add_argument("-m", "--metadata_url", dest="metadata_url", type=str, default="https://song.rdpc.cancercollaboratory.org")
    parser.add_argument("-s", "--storage_url", dest="storage_url", type=str, default="https://score.rdpc.cancercollaboratory.org")
    parser.add_argument("-n", "--cpu_number", dest="cpu_number", type=str, default=4)
    parser.add_argument("-c", "--conf", dest="conf", type=str, default="conf/af_only_gnomad.conf")
    parser.add_argument("-f", "--af_threshold", dest="af_threshold", type=list, default=[0, 0.001, 0.01, 0.1])
    parser.add_argument("-t", "--token", dest="token", type=str, required=True)
    args = parser.parse_args()

    song_dump = args.dump_path
    variant_calling_stats = {}

    #download qc_metrics
    #download(song_dump, 'qc_metrics', args.token, args.metadata_url, args.storage_url)

    variant_calling_stats = process_qc_metrics(song_dump, variant_calling_stats)

    pcawg_sample_sheet = args.pcawg_sample_sheet
    pcawg_sanger_qc = args.pcawg_sanger_qc
    pcawg_broad_qc = args.pcawg_broad_qc
    variant_calling_stats = add_pcawg_info(variant_calling_stats, pcawg_sample_sheet, pcawg_sanger_qc, pcawg_broad_qc)


    # download timing-supplement
    #download(song_dump, 'timing_metrics', args.token, args.metadata_url, args.storage_url)

    # process the timing-supplement
    #variant_calling_stats = process_timing_supplement(variant_calling_stats)
    
    # download snv vcf
    #download(song_dump, 'snv', args.token, args.metadata_url, args.storage_url)

    # download indel vcf
    #download(song_dump, 'indel', args.token, args.metadata_url, args.storage_url)

    # annotate the vcf with gnomad AF
    #data_dir = "data/variant_calling"
    #annot_dir = os.path.join("data", "annot_vcf")
    #annot_vcf(args.cpu_number, args.conf, data_dir, annot_dir)

    # process the annot_vcf
    #variant_calling_stats = process_annot_vcf(variant_calling_stats, args.af_threshold)

    report_dir = 'report'
    if not os.path.exists(report_dir):
        os.makedirs(report_dir)
    with open(os.path.join(report_dir, 'variant_calling_stats.json'), 'w') as f:
        f.write(json.dumps(variant_calling_stats, indent=2))

    # generate tsv file
    variant_calling_stats_tsv = []
    pcawg_qc_tsv = []
    for d, v in variant_calling_stats.items():
        variant_calling_stats_tsv.append(get_dict_value(None, v, variant_calling_stats_fields))
        pcawg_qc_tsv.append(get_dict_value(None, v, pcawg_qc_fields))
    report(variant_calling_stats_tsv, os.path.join(report_dir, 'variant_calling_stats.tsv'))
    report(pcawg_qc_tsv, os.path.join(report_dir, 'pcawg_qc.tsv'))


if __name__ == "__main__":
    main()
