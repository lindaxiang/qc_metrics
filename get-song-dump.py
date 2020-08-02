#!/usr/bin/env python3

import requests
import json


def download(song_url, file_dump):
    response = requests.get(song_url + '/studies/all')
    if response.status_code == 200:
        studies = response.json()
    else:
        raise Exception(response.text)

    with open(file_dump, 'w') as fp:
        for study in studies:
            if study in ['TEST-CA']: continue
            analyses = requests.get(song_url + ('/studies/%s/analysis?analysisStates=PUBLISHED' % (study))).json()
            for analysis in analyses:
                fp.write(json.dumps(analysis)+"\n")

    return

def main():
    song_url = "https://song.rdpc.cancercollaboratory.org"
    file_dump = "data/rdpc-song.jsonl"
    download(song_url, file_dump)

if __name__ == "__main__":
    main()
