# -*- coding: utf-8 -*-
#!/usr/bin/env python3

import os
import io
from Bio import SeqIO
from tqdm import tqdm
import logging
import re
from bs4 import BeautifulSoup as Bp
import requests
from multiprocessing import Pool
import argparse
import sys
from ncbi_accession_download import EntrezDownloader
from ncbi_accession_download import DownloadJob
# from EntrezDownloader import EntrezDownloader
# from jobs import DownloadJob


def parseArgs():
    parser = argparse.ArgumentParser(description='Parallel downloading genome/nucleotide/protien sequences from NCBI using Entrez.')

    parser.add_argument('-i', '--input', type = str, required=True, 
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
                        help='Enables a progress bar, requires tqdm package')
    parser.add_argument('-u', '--uid', type=bool, required=False, default=False,
                        help='When input file is a uid list.')
    parser.add_argument('-d', '--database', type=str, required=False, 
                        help='When uid is given or input is a str, --databse should be \'assembly\'/\'protein\'/\'nucleotide\',')
    parser.add_argument('--retmax', type=int, required=False, default=0,
                        help='Searching mode max retrieving sequences, default=[batch*2], maxmium 100,000')
    parser.add_argument('--dry', action=argparse.BooleanOptionalAction, default= False,
                        help='Download part not running, just find the uid.')
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


def config(args):
    infile = args.input
    output = args.output

    tmp = f'{output}/tmp'
    if not os.path.exists(tmp):
        os.makedirs(tmp)

    edl = edl_config(args.email, args.api,
                     args.parallel, args.batch, args.progress_bar)
    
    return infile, tmp, edl


def get_uid(ids, _database, edl):

    logger = logging.getLogger('ncbi-accession-download')
    logger.info('Sorting the input accession number to diffenent database and changing it to uid...')
    
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
    logger = logging.getLogger('ncbi-accession-download')

    logger.info('\nReporting the download sequence number in each database...')
    for k, v in _dict.items():
        print(f'{k}: {len(v)-1}')
    logger.info('End reporting, the file including uids are stored at output/tmp.\n')
    _dict = {k:v.remove('Head') for k,v in _dict.items()}


def download_assem(_assem_uids, _output, _edl, args):
    logger = logging.getLogger('ncbi-accession-download')

    logger.info('genome sequences downloading...')
    logger.info('Getting the NCBI ftp address')
    r_fetch, f_fetch = _edl.efetch(db='assembly', ids=_assem_uids, retmode='xml', retype='docsum',
                                    result_func = lambda x:re.findall(r'<FtpPath_GenBank>(ftp:.+)</FtpPath_GenBank>',x))

    logger.error('Falut number: %s\n%s' % (len(f_fetch), f_fetch)) if not len(f_fetch) == 0 else logger.info('No failed')
    if len(r_fetch) == 0:
        return 1

    url_list = []
    for ftp_url in r_fetch:
        gca_id = ftp_url.split('/')[-1]
        down_url = f'{ftp_url}/{gca_id}_genomic.fna.gz'.replace('ftp://','http://')
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
                logger.debug("Interrupted by user")
                return 1
    except requests.exceptions.ConnectionError as err:
        logger.error('Download from NCBI failed: %r', err)
        return 75
    return 0
    

def worker(job):
    logger = logging.getLogger('ncbi-accession-download')

    try:
        req = requests.get(job.full_url,stream=True)
        with open(f'{job.output}/{job.filename}','wb') as handle:
            for chunk in req.iter_content(4096):
                handle.write(chunk)
    except KeyboardInterrupt:
        logger.debug("Ignoring keyboard interrupt.")

    return True


def downloadjob_creator_caller(args):
    return downloadjob_creator(*args)


def downloadjob_creator(url, output):
    filename = url.split('/')[-1]
    download_jobs = [DownloadJob(url, filename, output)]
    return download_jobs


def download_edl(uids,db,output,edl):
    logger = logging.getLogger('ncbi-accession-download')

    logger.info(f'{db} sequences downloading...')
    try:
        r_fetch, f_fetch = edl.efetch(db=db, ids=uids, retype='fasta', retmode='text',
                                  result_func=lambda x: [SeqIO.read(io.StringIO(f), 'fasta')
                                                         for f in x.split('\n\n')[:-1]])
    except KeyboardInterrupt:
        logger.error("Interrupted by user")
        return 1

    logger.error('Falut number: %s\n%s' % (len(f_fetch), f_fetch)) if not len(f_fetch) == 0 else logger.info('No failed')
    if len(r_fetch) == 0:
        return 1
    
    outfile = f'{output}/{db}.fasta'
    count = SeqIO.write(r_fetch,outfile,'fasta')
    logger.info(f"{count} sequences have saved to {outfile}")
    return 0


def merge_dict(dict1,dict2):
    new_dict = {}
    for k, v in dict2.items():
        if k != 'lost':
            new_v = list(set(dict1[k] + v))
            new_dict[k] = new_v
    new_dict['lost'] = dict2['lost']
    return new_dict


def download_part(sorted_dict, edl, args):
    logger = logging.getLogger('ncbi-accession-download')

    if args.dry:
        logger.info('DRY RUN. exiting...')
        sys.exit()

    for db,uids in sorted_dict.items():
        if db != 'lost':
            if db == 'assembly':
                state = download_assem(uids,args.output,edl,args)
            else:
                state = download_edl(uids,db,args.output,edl)

            if state:
                logger.error(f'{db} download Failed')
            else:
                logger.info(f'{db} download Successful')


def download_acc(args,database):
    logger = logging.getLogger('ncbi-accession-download')

    infile, tmp, edl = config(args)
    
    acc_ids = open(infile).read().split('\n')
    # delete the last empty line
    if not acc_ids[-1]:
        acc_ids = acc_ids[:-1]
    
    if args.uid:
        if args.database in database:
            sorted_dict = {args.database:acc_ids}
        else:
            sys.exit("When the uid list as an inpufile, \'--database\' should be gaven correctly.\n \"nad --help\".")
    
    else:
        # initialize sorted_input
        sorted_dict = get_uid(acc_ids, database, edl)

        if len(sorted_dict['lost']) > 0:
            logger.warning(f"There are {len(sorted_dict['lost'])} accession number not found in database.")
            logger.info('setting batch size to 1 for lost accession number')

            _edl = edl_config(args.email, args.api,
                            args.parallel, 1, args.progress_bar)
            more_sorted = get_uid(sorted_dict['lost'][1:],database,_edl)

            logger.warning('merging to a new uid dictory')
            sorted_dict = merge_dict(sorted_dict,more_sorted)
        
        for k,v in sorted_dict.items():
            with open(f'{tmp}/uid_{k}.tsv', 'w') as f:
                f.write('\n'.join(v))
    
    report_divide(sorted_dict)

    download_part(sorted_dict, edl, args)


def download_spe(args,database):
    logger = logging.getLogger('ncbi-accession-download')

    infile, tmp, edl = config(args)

    if not args.database in database:
        sys.exit("When the uid list as an inpufile, \'--database\' should be gaven correctly.\n Check \"nad --help\".")

    if not args.retmax:
        retmax = args.batch * 2
    else:
        retmax = args.retmax
    if retmax > 100000:
        retmax = 100000

    logger.info(f'Searching the term \"{infile}\"')
    r_search, f_search = edl.esearch(args.database,infile,retmax)

    result_count = re.search(r'<Count>([0-9]+)', r_search[0][1])[1]

    if result_count:
        logger.info(f'Find {result_count} records')
    else:
        logger.warning(f'Find {result_count} records')
        sys.exit()

    if int(result_count) > retmax:
        logger.warning(f'The retrieving record {retmax} is samller than records, the result maybe incomplete')

    uids = re.findall(r'<Id>([0-9]+)', r_search[0][1])
    
    sorted_dict = {args.database:uids}
    with open(f'{tmp}/uid_{args.database}_{infile}.tsv', 'w') as f:
        f.write('\n'.join(uids))

    download_part(sorted_dict, edl, args)


def main():
    args = parseArgs()

    infile = args.input

    # assembly must be at the FIRST position, or some MAG will be more when it goes nucleotide protein database
    database = ['assembly', 'nucleotide', 'protein']

    logger = logging.getLogger("ncbi-genome-download")
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    if os.path.isfile(infile):
        download_acc(args,database)
    elif '/' not in infile and not os.path.isdir(infile):
        download_spe(args,database)


if __name__ == '__main__':
    main()
    