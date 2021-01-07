#!/usr/bin/env python3

import json
import os
import requests
import json
from argparse import ArgumentParser
import time
import sys
import csv
import ast

def song_operation(endpoint, operation, token, data=None):
    
    headers = {
        "Authorization": "Bearer %s" % token,
        "Content-Type": "application/json",
        "Accept": "application/json"
    }

    # wait a bit
    time.sleep(.5)
    try:
        if data is None:
            res = requests.put(endpoint, headers=headers)
        else:
            res = requests.put(endpoint, data=data, headers=headers)
        res.raise_for_status()
        
    except requests.exceptions.HTTPError as err:
        sys.exit("SONG %s failed, HTTPError: %s" % (operation, err))
    except requests.exceptions.RequestException as err:
        sys.exit("SONG %s failed, RequestException: %s" % (operation, err))

    if res.status_code != 200:
        sys.exit("SONG %s failed HTTP status code not 200: %s" % (operation, res))


def main():
    parser = ArgumentParser()
    parser.add_argument("-r", "--report", dest="report", type=str, help="report file containing analysis info to suppress")
    parser.add_argument("-l", "--analysis_list", dest="analysis_list", type=str, help="file containing analysis info to suppress")
    parser.add_argument("-m", "--song_url", dest="song_url", type=str, default="https://song.rdpc.cancercollaboratory.org", help="SONG URL")
    parser.add_argument("-s", "--score_url", dest="score_url", type=str, default="https://score.rdpc.cancercollaboratory.org", help="SCORE URL")
    parser.add_argument("-t", "--token", dest="token", type=str, required=True)
    args = parser.parse_args()

    if args.report:
        with open(args.report, 'r') as fp:
            reader = csv.DictReader(fp, delimiter='\t')
            for run in reader:
                alist = ast.literal_eval(run.get('run_output_analysis_to_suppress')) 
                for a in alist:
                    print(a)
                    endpoint = "%s/studies/%s/analysis/suppress/%s" % (args.song_url, run.get('studyId'), a)
                    operation = 'analysis_suppress'
                    song_operation(endpoint, operation, args.token)
    
    if args.analysis_list:
        with open(args.analysis_list, 'r') as fp:
            reader = csv.DictReader(fp, delimiter='\t')
            for run in reader:
                print(run.get('analysisId'))
                endpoint = "%s/studies/%s/analysis/suppress/%s" % (args.song_url, run.get('studyId'), run.get('analysisId'))
                operation = 'analysis_suppress'
                song_operation(endpoint, operation, args.token)
    

if __name__ == "__main__":
    main()