# -*- coding: utf-8 -*-
#!/usr/bin/env python3

import os
import io
from EntrezDownloader import EntrezDownloader
from Bio import SeqIO
from tqdm import tqdm
# from glob import glob
# import re
from bs4 import BeautifulSoup as Bp
import requests
from multiprocessing import Pool
import argparse

from jobs import DownloadJob


def parseArgs():
    parser = argparse.ArgumentParser(description='Parallel downloading genome/nucleotide/protien sequences from NCBI using Entrez.')

    parser.add_argument('-i', '--input', type = str, required=True, default='test/' ,
                        help = 'List file with accession number seperated by \\n.')
    parser.add_argument('-o', '--output', type = str, required=True, help = 'Downloaded files saved directory.')
    parser.add_argument('-e', '--email', type = str, required=False, default = 'eacochen@163.com' , 
                        help = 'An email address. You might get blocked by the NCBI without specifying one.')
    parser.add_argument('-a', '--api', type=str, required=False, default='df1607793b0171d9f7fd39585f90a3491007', 
                        help='An API key. You can obtain one by creating an NCBI account. Speeds things up.')
    parser.add_argument('-p', '--parallel', type=int, required=False, default = 10, 
                        help='The number of parallel requests to make')
    parser.add_argument('-b', '--batch', type=int, required=False, default = 20, 
                        help='The number of IDs to fetch per request')
    parser.add_argument('-P', '--progress-bar', dest='progress_bar', type=bool, required=False, default=True,
                        help = 'Enables a progress bar, requires tqdm package')

    return parser.parse_args()


def edl_config(user_email,api,parallel,batch,progress_bar):

    edl = EntrezDownloader(
        email=user_email,
        api_key=api,
        num_threads=parallel,
        batch_size=batch,
        pbar=progress_bar
    )

    return edl


def get_uid(ids, _database, edl):
    print('Sorting the input accession number to diffenent database and changing it to uid...')
    
    ids_dict = {}
    for _db in _database:
        tqdm.write(_db)
        r_search, f_search = edl.esearch(_db, ids)

        # initialize container
        lost_ids = []
        ids_dict[_db] = ['Head']

        for result in r_search:
            if type(result) == tuple:
                _soup = Bp(result[1], 'xml')
            else:
                _soup = Bp(result, 'xml')
            lost_ids += [_.get_text() for _ in _soup.find_all('PhraseNotFound')]
            lost_ids += [_.find('Term').get_text().split('[')[0] for _ in _soup.find_all('TermSet')
                         if int(_.find('Count').get_text()) > 1]
            ids_dict[_db] += [_.get_text() for _ in _soup.find_all('Id')]

        ids = lost_ids

    ids_dict['lost'] = ['Head'] + lost_ids
    return ids_dict


def report_divide(_dict):
    print('\nReporting the download sequence number in each database...')
    for k, v in _dict.items():
        print(f'{k}: {len(v)-1}')
    print('End reporting, the file including uids are stored at output/tmp.\n')
    # return _dict[_database[0]][1:],_dict[_database[1]][1:],_dict[_database[2]][1:]


def download_assem(_assem_uids, _output, _edl, args):
    print('genome sequences downloading...')
    print('Getting the NCBI ftp address')
    r_fetch, f_fetch = _edl.efetch(db='assembly', ids=_assem_uids, 
                                    retmode='xml', retype='docsum')

    print('Falut number: %s\n%s' % (len(f_fetch), f_fetch)) if not len(f_fetch) == 0 else print('No failed')
    if len(r_fetch) == 0:
        return 1

    genbank_list = []
    for result_fetch in r_fetch:
        soup = Bp(result_fetch, 'xml')
        genbank_list += (soup.find_all('FtpPath_GenBank'))

    url_list = []
    for ftp_url in genbank_list:
        gca_id = ftp_url.get_text().split('/')[-1]
        down_url = f'{ftp_url.get_text()}/{gca_id}_genomic.fna.gz'.replace('ftp://','http://')
        url_list.append(down_url)
    
    output = f'{_output}/assembly'
    if not os.path.exists(output):
        os.makedirs(output)

    download_jobs = []
    try:
        with Pool(processes=args.parallel) as pool:
            dl_jobs = pool.imap(
                downloadjob_creator_caller, [(url, output) for url in url_list]
            )

            for index,create_dl_job in enumerate(dl_jobs):
                download_jobs.extend(create_dl_job)
            
            jobs = [pool.apply_async(worker, (_,)) for _ in download_jobs]

            try:
                if args.progress_bar:
                    tqdm.write(f"\nStart download. Parallel={args.parallel}")
                    _jobs = tqdm(jobs)
                else:
                    _jobs=jobs
                [_.get(0xFFFF) for _ in _jobs]
            except KeyboardInterrupt:
                print("Interrupted by user")
                return 1
    except requests.exceptions.ConnectionError as err:
        print('Download from NCBI failed: %r', err)
        return 75
    return 0
    

def worker(job):
    try:
        req = requests.get(job.full_url,stream=True)
        with open(f'{job.output}/{job.filename}','wb') as handle:
            for chunk in req.iter_content(4096):
                handle.write(chunk)
    except KeyboardInterrupt:
        print("Ignoring keyboard interrupt.")

    return True


def downloadjob_creator_caller(args):
    return downloadjob_creator(*args)


def downloadjob_creator(url, output):
    filename = url.split('/')[-1]
    download_jobs = [DownloadJob(url, filename, output)]
    return download_jobs


def download_edl(uids,db,output,edl):
    print(f'{db} sequences downloading...')
    try:
        r_fetch, f_fetch = edl.efetch(db=db, ids=uids, retype='fasta', retmode='text',
                                  result_func=lambda x: [SeqIO.read(io.StringIO(f), 'fasta')
                                                         for f in x.split('\n\n')[:-1]])
    except KeyboardInterrupt:
        print("Interrupted by user")
        return 1

    print('Falut number: %s\n%s' % (len(f_fetch), f_fetch)) if not len(f_fetch) == 0 else print('No failed')
    if len(r_fetch) == 0:
        return 1
    
    outfile = f'{output}/{db}.fasta'
    count = SeqIO.write(r_fetch,outfile,'fasta')
    print(f"{count} sequences have saved to {outfile}")
    return 0


def merge_dict(dict1,dict2):
    new_dict = {}
    for k, v in dict2.items():
        if k != 'lost':
            new_v = list(set(dict1[k] + v))
            new_dict[k] = new_v
    new_dict['lost'] = dict2['lost']
    return new_dict


def main():
    args = parseArgs()

    infile = args.input
    output = args.output
    parallel = args.parallel

    tmp = f'{output}/tmp'
    if not os.path.exists(tmp):
        os.makedirs(tmp)

    edl = edl_config(args.email, args.api,
                     parallel, args.batch, args.progress_bar)
    
    acc_ids = open(infile).read().split('\n')
    # delete the last empty line
    if not acc_ids[-1]:
        acc_ids = acc_ids[:-1]
    
    # assembly must be at the FIRST position, or some MAG will be more when it goes nucleotide protein database
    database = ['assembly', 'nucleotide', 'protein']

    # initialize sorted_input
    sorted_dict = get_uid(acc_ids, database, edl)

    if len(sorted_dict['lost']) > 0:
        print('\nsetting batch size to 1 for lost accession number')
        _edl = edl_config(args.email, args.api,
                         parallel, 1, args.progress_bar)
        more_sorted = get_uid(sorted_dict['lost'][1:],database,_edl)

        print('merging to a new uid dictory')
        sorted_dict = merge_dict(sorted_dict,more_sorted)
    
    for k,v in sorted_dict.items():
        with open(f'{tmp}/uid_{k}.tsv', 'w') as f:
            f.write('\n'.join(v))
    
    report_divide(sorted_dict)

    for db,uids in sorted_dict.items():
        if db != 'lost':
            if db == 'assembly':
                state = download_assem(uids,output,edl,args)
            else:
                state = download_edl(uids,db,output,edl)

            if state:
                print(f'{db} download Failed')
            else:
                print(f'{db} download Successful')


if __name__ == '__main__':
    main()
    