import json
import os
import csv
import glob
import sys
import subprocess
from collections import OrderedDict
import functools
import re
import pandas as pd
import numpy as np


def report(donor, report_name):
    keys = donor[0].keys()
    with open(report_name, 'w') as output_file:
        dict_writer = csv.DictWriter(output_file, keys, delimiter="\t")
        dict_writer.writeheader()
        dict_writer.writerows(donor)

def run_cmd(cmd):
    try:
        p = subprocess.run([cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           shell=True, check=True)

    except subprocess.CalledProcessError as e:   # this is triggered when cmd returned non-zero code
        print(e.stdout.decode("utf-8"))
        print('Execution returned non-zero code: %s. Additional error message: %s' %
              (e.returncode, e.stderr.decode("utf-8")), file=sys.stderr)
        sys.exit(e.returncode)

    except Exception as e:  # all other errors go here, like unable to find the command
        sys.exit('Execution failed: %s' % e)

    return p  # in case the caller of this funtion needs p.stdout, p.stderr etc

def udf(x, y):
    if not x:
        return
    if isinstance(x, dict):
        value = x.get(y)
    elif isinstance(x, list):
        value = [udf(a,y) for a in x]
    return value

def get_dict_value(fields, json_obj, field_map):
    tsv_obj = OrderedDict()
    if fields is None: fields = field_map.keys()
    for f in fields:
        fm = field_map.get(f)
        if not fm: continue
        value = functools.reduce(udf, fm.split('.'), json_obj)
        tsv_obj[f] = value
    return tsv_obj 

def download(song_dump, file_type, ACCESSTOKEN, METADATA_URL, STORAGE_URL, include=None, subfolder=None):

    file_type_map = { # [analysisType, dataType, data_category]
        "qc_metrics": ['qc_metrics', 'Alignment QC', 'Quality Control Metrics'],
        "timing_metrics": ['variant_calling_supplement', 'Variant Calling Supplement', None],
        "snv": ['variant_calling', 'Raw SNV Calls', 'Simple Nucleotide Variation'],
        "indel": ['variant_calling', 'Raw InDel Calls', 'Simple Nucleotide Variation'],
        "sv": ['variant_calling', 'Raw SV Calls', 'Structural Variation'],
        "cnv": ['variant_calling', 'Raw CNV Calls', 'Copy Number Variation']
    }

    if subfolder:
        data_dir = os.path.join("data", subfolder)
    else:
        data_dir = os.path.join("data", file_type_map[file_type][0])
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    downloaded = []
    for fn in glob.glob(os.path.join(data_dir, "*-*", "*.*"), recursive=True):
        downloaded.append(os.path.basename(fn))

    with open(song_dump, 'r') as fp:
        for fline in fp:
            analysis = json.loads(fline)
            if include and not analysis['analysisId'] in include: continue
            if analysis.get('analysisId') in ['fa2cc829-c646-4719-acc8-29c6465719e4']: continue
            if not analysis.get('analysisState') == 'PUBLISHED': continue
            if not analysis['analysisType']['name'] == file_type_map[file_type][0]: continue
            output_dir = os.path.join(data_dir, analysis['studyId'])
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            for fl in analysis['files']:
                if fl['fileName'] in downloaded: continue
                if not file_type in fl['fileName']: continue
                if not fl['dataType'] == file_type_map[file_type][1]: continue
                if file_type_map[file_type][2] is None and 'data_category' in fl['info']: continue
                if file_type_map[file_type][2] and not fl['info']['data_category'] == file_type_map[file_type][2]: continue
                

                cmd = 'export ACCESSTOKEN=%s && export METADATA_URL=%s \
                    && export STORAGE_URL=%s && export TRANSPORT_PARALLEL=3 \
                    && export TRANSPORT_MEMORY=8 \
                    && docker run -it --rm  -u $(id -u):$(id -g) \
                    -e ACCESSTOKEN -e METADATA_URL -e STORAGE_URL \
                    -e TRANSPORT_PARALLEL -e TRANSPORT_MEMORY \
                    -v "$PWD":"$PWD" -w "$PWD" overture/score:latest /score-client/bin/score-client download \
                    --study-id %s --object-id %s --output-dir %s' \
                    % (ACCESSTOKEN, METADATA_URL, STORAGE_URL, fl['studyId'], fl['objectId'], output_dir)

                run_cmd(cmd)


def annot_vcf(cores, conf, data_dir, annot_dir):

    annotated = []
    for fn in glob.glob(os.path.join(annot_dir, "*-*", "*.*"), recursive=True):
        annotated.append(os.path.basename(fn))

    for fp in glob.glob(os.path.join(data_dir, "*-*", "*.vcf.gz"), recursive=True):
        basename = os.path.basename(fp)
        if basename in annotated: continue

        study_id = basename.split('.')[0]
        if not os.path.exists(os.path.join(annot_dir, study_id)):
            os.makedirs(os.path.join(annot_dir, study_id))
        vcf = os.join(annot_dir, study_id, basename)

        vcfanno = f'vcfanno -p {cores} {conf} {fp}'
        bgzip = f'bgzip > {vcf}'
        tabix = f'tabix -p vcf {vcf}'
        cmd = vcfanno + ' | ' + bgzip + ' && ' + tabix
        run_cmd(cmd)
    
    for fp in glob.glob(os.path.join(annot_dir, "*-*", "*.vcf.gz"), recursive=True):
        # use bcftools to query the annotated vcf
        basename = os.path.basename(fp)
        if basename+'.query.txt' in annotated: continue
        bcftools_query(fp)

def bcftools_query(vcf):
    basename = os.path.basename(vcf)
    caller = 'sanger' if basename.split('.')[5] in ['sanger-wgs', 'sanger-wxs'] else 'mutect2'
    evtype = basename.split('.')[7]
    output_base = re.sub(r'.vcf.gz$', '', vcf)

    if caller == 'sanger':
        if evtype == 'snv':
            cmd = f"bcftools query -f '[%CHROM\\t%POS\\t%REF\\t%ALT\\t%PM\\t%gnomad_af\\t%gnomad_filter\\n]' \
                -i 'FILTER=\"PASS\"' -s 'TUMOUR' {vcf} > {output_base}.query.txt"
        elif evtype == 'indel':
            bcftools = f"bcftools query -f '[%CHROM\\t%POS\\t%REF\\t%ALT\\t%FD\\t%FC\\t%gnomad_af\\t%gnomad_filter\\n]' \
                -i 'FILTER=\"PASS\"' -s 'TUMOUR' {vcf}"
            awk = f"awk '{{printf \"%s\\t%s\\t%s\\t%s\\t%f\\t%s\\t%s\\n\",$1,$2,$3,$4,$6/$5,$7,$8}}' > {output_base}.query.txt"
            cmd = bcftools + '|' + awk
        else:
            pass
    elif caller == 'mutect2':
        cmd = f"bcftools query -f '[%CHROM\\t%POS\\t%REF\\t%ALT\\t%AF\\t%gnomad_af\\t%gnomad_filter\\n]' -i 'FILTER=\"PASS\" & GT=\"0/1\"' {vcf} > {output_base}.query.txt" 
    else:
        pass
    
    run_cmd(cmd)