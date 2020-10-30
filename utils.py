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
from datetime import date


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


def annot_vcf(cores, conf, data_dir, annot_dir, bed_dir=None):

    annotated = []
    for fn in glob.glob(os.path.join(annot_dir, "*-*", "*.*"), recursive=True):
        annotated.append(os.path.basename(fn))
    
    #annot gnomAD 
    for fp in glob.glob(os.path.join(data_dir, "*-*", "*.vcf.gz"), recursive=True):
        basename = os.path.basename(fp)
        if basename in annotated: continue

        study_id = basename.split('.')[0]
        if not os.path.exists(os.path.join(annot_dir, study_id)):
            os.makedirs(os.path.join(annot_dir, study_id))
        vcf = os.path.join(annot_dir, study_id, basename)

        vcfanno = f'vcfanno -p {cores} {conf} {fp}'
        bgzip = f'bgzip > {vcf}'
        tabix = f'tabix -p vcf {vcf}'
        cmd = vcfanno + ' | ' + bgzip + ' && ' + tabix
        run_cmd(cmd)
    
    #use bcftools to query the annotated vcf
    for fp in glob.glob(os.path.join(annot_dir, "*-*", "*.vcf.gz"), recursive=True):
        basename = os.path.basename(fp)
        if re.sub(r'.vcf.gz$', '', basename) + '.query.txt' in annotated: continue
        bcftools_query(fp, bed_dir)

    #concatenate the query results for each caller
    for fp in glob.glob(os.path.join(annot_dir, "*-*", "*.query.txt"), recursive=True):
        prefix = os.path.basename(fp).split("2020")[0].replace(".", "_")
        evtype = os.path.basename(fp).split(".")[7]
        cat = f'cat {fp}'
        awk = f'awk \'{{printf "%s_%s_%s_%s\\t%f\\t%s\\t%s\\n\",$1,$2,$3,$4,$5,$6,$7}}\''
        sed = f'sed "s/^/{prefix}/g" >> {annot_dir}.{evtype}.all'
        cmd = '|'.join([cat, awk, sed])
        run_cmd(cmd)
    

def bcftools_query(vcf, bed_dir=None):
    basename = os.path.basename(vcf)
    donorId, sampleId = basename.split('.')[1:3]
    caller = 'sanger' if basename.split('.')[5] in ['sanger-wgs', 'sanger-wxs'] else 'mutect2'
    evtype = basename.split('.')[7]
    output_base = re.sub(r'.vcf.gz$', '', vcf)
    
    if bed_dir:
        if evtype == "snv": evtype_str = evtype + "_mnv" 
        bed_filename = os.path.join(bed_dir, '.'.join([donorId, evtype_str+'_inflated', 'bed']))    
        bcftools = f"bcftools query -R {bed_filename} " if os.path.exists(bed_filename) else f"bcftools query "        
    else:
        bcftools = f"bcftools query "

    if caller == 'sanger':
        if evtype == 'snv':
            cmd = bcftools + f"-f '[%CHROM\\t%POS\\t%REF\\t%ALT\\t%PM\\t%gnomad_af\\t%gnomad_filter\\n]' \
                -i 'FILTER=\"PASS\"' -s 'TUMOUR' {vcf} > {output_base}.query.txt"
        elif evtype == 'indel':
            bcftools = bcftools + f"-f '[%CHROM\\t%POS\\t%REF\\t%ALT\\t%FD\\t%FC\\t%gnomad_af\\t%gnomad_filter\\n]' \
                -i 'FILTER=\"PASS\"' -s 'TUMOUR' {vcf}"
            awk = f"awk '{{printf \"%s\\t%s\\t%s\\t%s\\t%f\\t%s\\t%s\\n\",$1,$2,$3,$4,$6/$5,$7,$8}}' > {output_base}.query.txt"
            cmd = bcftools + '|' + awk
        else:
            pass
    elif caller == 'mutect2':
        cmd = bcftools + f"-f '[%CHROM\\t%POS\\t%REF\\t%ALT\\t%AF\\t%gnomad_af\\t%gnomad_filter\\n]' -i 'FILTER=\"PASS\" & GT=\"0/1\"' {vcf} > {output_base}.query.txt" 
    else:
        pass

    run_cmd(cmd)

def get_info(row):
    info = []
    for col in ['Caller', 'wgsTumourVAF', 'gnomadAF']:
        if pd.isna(row[col]): continue
        info.append(col+"="+str(row[col]))
    return ';'.join(info)

def union_vcf(data_dir, union_dir):
    if not os.path.exists(union_dir):
        os.makedirs(union_dir)
    donor = set()

    for fn in glob.glob(os.path.join(data_dir, "*_annot_vcf", "*-*", "*.query.txt"), recursive=True):
        donor.add(os.path.basename(fn).split("2020")[0].rstrip('.'))
   
    for evtype in ['snv', 'indel']:
        for do in donor:
            projectId, donorId, sampleId, library_strategy = do.split('.')
            df = None
            for caller in ['sanger', 'mutect2', 'mutect2-bqsr']:
                query_file = glob.glob(os.path.join(data_dir, caller+'_annot_vcf', projectId, do+'*'+evtype+'.query.txt'))
                df_caller = pd.read_table(query_file[0], sep='\t', \
                        names=["CHROM", "POS", "REF", "ALT", "AF_"+caller, "gnomad_af_"+caller, "gnomad_filter_"+caller], \
                        dtype={"CHROM": str, "POS": int, "REF": str, "ALT": str, "AF_"+caller: float, "gnomad_af_"+caller: float, "gnomad_filter_"+caller: str}, \
                        na_values=".")
                if df is None:
                    df = df_caller
                else:
                    df_all = df.join(df_caller.set_index(['CHROM', 'POS', 'REF', 'ALT']), on=['CHROM', 'POS', 'REF', 'ALT'], how='outer')
                    df = df_all

            # add columns: caller, VAF, VAF_level, gnomad_af, gnomad_af_level
            conditions = [
                df_all['AF_sanger'].notna() & df_all['AF_mutect2'].isna() & df_all['AF_mutect2-bqsr'].isna(),
                df_all['AF_sanger'].isna() & df_all['AF_mutect2'].notna() & df_all['AF_mutect2-bqsr'].isna(),
                df_all['AF_sanger'].isna() & df_all['AF_mutect2'].isna() & df_all['AF_mutect2-bqsr'].notna(),
                df_all['AF_sanger'].notna() & df_all['AF_mutect2'].notna() & df_all['AF_mutect2-bqsr'].isna(),
                df_all['AF_sanger'].notna() & df_all['AF_mutect2'].isna() & df_all['AF_mutect2-bqsr'].notna(),
                df_all['AF_sanger'].isna() & df_all['AF_mutect2'].notna() & df_all['AF_mutect2-bqsr'].notna(),
                df_all['AF_sanger'].notna() & df_all['AF_mutect2'].notna() & df_all['AF_mutect2-bqsr'].notna()
            ]
            caller_outputs = ['sanger', 'mutect2', 'mutect2-bqsr', 'sanger,mutect2', 'sanger,mutect2-bqsr', 'mutect2,mutect2-bqsr', 'sanger,mutect2,mutect2-bqsr']
            df_all['Caller'] = np.select(conditions, caller_outputs, 'other')

            VAF = df_all.loc[:, ['AF_sanger', 'AF_mutect2', 'AF_mutect2-bqsr']].mean(axis=1)
            df_all['wgsTumourVAF'] = VAF
            gnomad_af = df_all.loc[:, ['gnomad_af_sanger', 'gnomad_af_mutect2', 'gnomad_af_mutect2-bqsr']].mean(axis=1)
            df_all['gnomadAF'] = gnomad_af
            # format INFO field
            cols = ['Caller', 'wgsTumourVAF', 'gnomadAF']
            df_all['INFO'] = df[cols].apply(get_info, axis=1)

            # generate vcf from dataframe
            vcf_file = os.path.join(union_dir, '.'.join([do, evtype, 'vcf']))
            date_str = date.today().strftime("%Y%m%d")
            header = [f"##fileformat=VCFv4.3", \
                      f"##fileDate={date_str}", \
                      f"#CHROM POS ID REF ALT QUAL FILTER INFO"]
            with open(vcf_file, 'w') as vcf:
                vcf.write("\n".join(header))
                vcf.write("\n")
            cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
            df_all['ID'] = ""
            df_all['QUAL'] = ""
            df_all['FILTER'] = ""
            df_all.to_csv(vcf_file, index=False, header=False, sep="\t", mode='a', columns=cols)

            # generate bed from dataframe
            bed_file = os.path.join(union_dir, '.'.join([do, evtype, 'bed']))
            df_all['START'] = df['POS']
            df_all['END'] = df['POS'] + 1
            cols = ['CHROM', 'START', 'END']
            df_all.to_csv(bed_file, index=False, header=False, sep="\t", columns=cols)


def bam_readcount(bam_dir, union_dir, readcount_dir, ref_fa):
    if not os.path.exists(readcount_dir):
        os.makedirs(readcount_dir)
    
    for fn in glob.glob(os.path.join(bam_dir, "*.aln.bam"), recursive=True):
        bam_rc_file = os.path.join(readcount_dir, os.path.basename(fn) + ".rc")
        projectId, donorId = os.path.basename(fn).split(".")[0:2]
        target_bed = glob.glob(os.path.join(union_dir, projectId+'.'+donorId+'.*.wgs.snv.bed'))[0]

        cmd = f"bam-readcount --reference-fasta {ref_fa} \
                  --site-list {target_bed} \
                  --min-mapping-quality 1 \
                  --max-count 8000 \
                  --max-warnings 0 {fn} > {bam_rc_file}"

        run_cmd(cmd)


def snv_readcound(union_dir, readcount_dir):
    sample = {
        "DO50311": {
            "normal": "SA622503",
            "tumour": "SA622502"
        },
        "DO50342": {
            "normal": "SA622504",
            "tumour": "SA622509"
        },  
        "DO50346": {
            "normal": "SA622506",
            "tumour": "SA622505"
        }, 
        "DO50407": {
            "normal": "SA622508",
            "tumour": "SA622507"
        }
    }
    for fn in glob.glob(os.path.join(union_vcf, "*.vcf"), recursive=True):
        projectId, donorId, sampleId, library_strategy, evtype, fileformat = os.path.basename(fn).split(".")
        output_vcf = os.path.join(union_dir, '.'.join([projectId, donorId, 'validated', evtype, fileformat]))
        normal_rc = glob.glob(os.path.join(readcount_dir, '.'.join([projectId, donorId, sample[donorId]['normal'], 'targeted-seq', '*', 'aln.bam.rc'])))[0]
        tumour_rc = glob.glob(os.path.join(readcount_dir, '.'.join([projectId, donorId, sample[donorId]['tumour'], 'targeted-seq', '*', 'aln.bam.rc'])))[0]

        header = """    
##fileformat=VCFv4.3
##INFO=<ID=Callers,Number=.,Type=String,Description="Callers that made this call">
##INFO=<ID=NormalEvidenceReads,Number=2,Type=Integer,Description="Number of reads in normal sample supporting the alt allele (forward,backward)">
##INFO=<ID=TumourEvidenceReads,Number=2,Type=Integer,Description="Number of reads in tumour sample supporting the alt allele (forward,backward)">
##INFO=<ID=TumourReads,Number=1,Type=Integer,Description="Total number of reads in tumour sample">
##INFO=<ID=NormalReads,Number=1,Type=Integer,Description="Total number of reads in normal sample">
##INFO=<ID=Validation_status,Number=.,Type=String,Description="Validation status, as info field">
##FILTER=<ID=LOWDEPTH,Description="Insufficient depth to validate the call">
##FILTER=<ID=NOTSEEN,Description="Variant not seen in Tumour">
##FILTER=<ID=STRANDBIAS,Description="Too much strand bias in Tumour sample to validate the call">
##FILTER=<ID=NORMALEVIDENCE,Description="Evidence for presence in normal sample">
##FILTER=<ID=GERMLINE,Description="Germline het or hom">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""
        with open(output_vcf, 'w') as vcf:
            vcf.write(header+"\n")

        snv_rc = f"./snv_readcounts.py {fn} {normal_rc} {tumour_rc}"
        snv_indel_call = f"./scripts/snv_indel_call.py --error 0.01 --callthreshold 0.02 --strandbias 0.001 --germlineprob 0.01 --mindepth 30"
        sed = f"sed -e 's/;;/;/'"
        grep = f'grep -v "^#"'
        sort = f'sort -k1,1 -k2,2n  >> {output_vcf}'
        cmd = '|'.join([snv_rc, snv_indel_call, sed, grep, sort])
        print(cmd)

        run_cmd(cmd)
