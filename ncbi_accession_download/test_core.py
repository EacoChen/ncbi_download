# -*- coding: utf-8 -*-
#!/usr/bin/env python3

import os
import io
from xml.etree.ElementInclude import default_loader
from Bio import SeqIO,Entrez
from tqdm import tqdm
import logging
import re
from bs4 import BeautifulSoup as Bp
import requests
from multiprocessing import Pool
import argparse
import sys
from collections import OrderedDict
from EntrezDownloader import EntrezDownloader
from jobs import DownloadJob


_FORMATS = OrderedDict([
        ('genbank', '_genomic.gbff.gz'),
        ('fasta', '_genomic.fna.gz'),
        ('rm', '_rm.out.gz'),
        ('features', '_feature_table.txt.gz'),
        ('gff', '_genomic.gff.gz'),
        ('protein-fasta', '_protein.faa.gz'),
        ('genpept', '_protein.gpff.gz'),
        ('wgs', '_wgsmaster.gbff.gz'),
        ('cds-fasta', '_cds_from_genomic.fna.gz'),
        ('rna-fna', '_rna.fna.gz'),
        ('rna-fasta', '_rna_from_genomic.fna.gz'),
        ('assembly-report', '_assembly_report.txt'),
        ('assembly-stats', '_assembly_stats.txt'),
    ])


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
    parser.add_argument('-u', '--uid', action='store_true',
                        help='When input file is a uid list.')
    parser.add_argument('-d', '--database', type=str, required=False, 
                        help='When uid is given or input is a str, --databse should ncbi entrez database,')
    parser.add_argument('--retmax', type=int, required=False, default=0,
                        help='Searching mode max retrieving sequences, default=[batch*2], maxmium 100,000')
    parser.add_argument('--dry', action='store_true',
                        help='Download part not running, just find the uid.')
    parser.add_argument('-f', '--format', type=str, required=False, default='all', 
                        help='Genome download format, multi format sepearte with comma')
    parser.add_argument('--retries', type = int, required = False, default=3,
                        help='Retry download time, if download failed. Default: 3')

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


def download_genome(_uids, _output, _edl, args):
    logger = logging.getLogger('ncbi-accession-download')

    logger.info('Convert genome ids to assembly ids.')
    r_link, f_link = _edl.elink(dbfrom='genome', db='assembly', ids=_uids,
                                result_func= lambda x:re.findall(r'\t\t\t\t<Id>([0-9]+)</Id>',x))
    
    logger.info(f'There are {len(r_link)} assembly records in this genome.')

    return download_assem(r_link, _output, _edl, args)


def download_assem(_assem_uids, _output, _edl, args):
    logger = logging.getLogger('ncbi-accession-download')

    logger.info('genome sequences downloading...')
    logger.info('Getting the NCBI ftp address')
    r_fetch, f_fetch = _edl.efetch(db='assembly', ids=_assem_uids, retmode='xml', retype='docsum',
                                    result_func = lambda x:re.findall(r'<FtpPath_GenBank>(ftp:.+)</FtpPath_GenBank>',x))

    logger.error('Falut number: %s\n%s' % (len(f_fetch), f_fetch)) if not len(f_fetch) == 0 else logger.info('No failed')
    if len(r_fetch) == 0:
        return 1

    output = f'{_output}/{args.database}'
    if not os.path.exists(output):
        os.makedirs(output)

    url_dict = {}
    for ftp_url in r_fetch:
        gca_id = ftp_url.split('/')[-1]

        if not os.path.exists(f"{output}/{gca_id}"):
            os.makedirs(f"{output}/{gca_id}")

        if args.format == 'all':
            for fmt in _FORMATS:
                down_url = f'{ftp_url}/{gca_id}{_FORMATS[fmt]}'.replace('ftp://','http://')
                url_dict[down_url] = gca_id
        else:
            for fmt in args.format.split(','):
                if fmt not in _FORMATS:
                    sys.exit('The format you given is not in database')
                else:
                    if os.path.exists(f"{output}/{gca_id}/{gca_id}{_FORMATS[fmt]}"):
                        continue
                    down_url = f'{ftp_url}/{gca_id}{_FORMATS[fmt]}'.replace('ftp://','http://')
                    url_dict[down_url] = gca_id

    download_jobs = []
    try:
        with Pool(processes=args.parallel) as pool:
            dl_jobs = pool.imap(
                downloadjob_creator_caller, [(url, f"{output}/{gca_id}") for url,gca_id in url_dict.items()]
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
        if not req.status_code == 404:
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
            elif db == 'genome':
                state = download_genome(uids,args.output,edl,args)
            else:
                state = download_edl(uids,db,args.output,edl)

            if state:
                logger.error(f'{db} download Failed')
            else:
                logger.info(f'{db} download Successful')
    
    return state


def download_acc(args,database):
    logger = logging.getLogger('ncbi-accession-download')

    infile, tmp, edl = config(args)
    
    acc_ids = open(infile).read().split('\n')
    # delete the last empty line
    if not acc_ids[-1]:
        acc_ids = acc_ids[:-1]
    
    if args.uid:
        if args.database in database:
            sorted_dict = {args.database:['Head']+acc_ids}
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

    return download_part(sorted_dict, edl, args)


def download_spe(args,database):
    logger = logging.getLogger('ncbi-accession-download')

    infile, tmp, edl = config(args)

    if not args.retmax:
        retmax = args.batch * 2
    else:
        retmax = args.retmax
    if retmax > 100000:
        retmax = 100000

    logger.info(f'Searching the term \"{infile}\"')
    r_search, f_search = edl.esearch(database, infile)

    result_count = re.search(r'<Count>([0-9]+)', r_search[0][1])[1]

    if result_count:
        logger.info(f'Find {result_count} records')
    else:
        logger.warning(f'Find {result_count} records')
        sys.exit()

    if int(result_count) > retmax:
        logger.warning(f'The retrieving record {retmax} is samller than records, the result maybe incomplete')

    uids = re.findall(r'<Id>([0-9]+)', r_search[0][1])
    
    sorted_dict = {database:uids}
    name = infile.replace(" ","_").replace(":","-").replace("[","(").replace("]",")")
    with open(f'{tmp}/uid_{database}_{name}.tsv', 'w') as f:
        f.write('\n'.join(uids))

    return download_part(sorted_dict, edl, args)


def main():
    args = parseArgs()

    infile = args.input

    # assembly must be at the FIRST position, or some MAG will be more when it goes nucleotide protein database
    # up to 2022.01.17 the database are as following:
    # all_database = ['pubmed', 'protein', 'nuccore', 'ipg', 'nucleotide', 'structure', 'genome', 
    #             'annotinfo', 'assembly', 'bioproject', 'biosample', 'blastdbinfo', 'books', 
    #             'cdd', 'clinvar', 'gap', 'gapplus', 'grasp', 'dbvar', 'gene', 'gds', 'geoprofiles', 
    #             'homologene', 'medgen', 'mesh', 'ncbisearch', 'nlmcatalog', 'omim', 'orgtrack', 
    #             'pmc', 'popset', 'proteinclusters', 'pcassay', 'protfam', 'biosystems', 
    #             'pccompound', 'pcsubstance', 'seqannot', 'snp', 'sra', 'taxonomy', 'biocollections', 'gtr']
    Entrez.email = args.email
    all_database = Entrez.read(Entrez.einfo())["DbList"]

    database = ['assembly', 'nucleotide', 'protein']
    
    args.database = args.database.lower()
    if args.database in all_database and ',' not in args.database:
        database.append(args.database)
    elif args.database:
        sys.exit(f'The -d database you given is not in ncbi entrez databases\n {",".join(all_database)}')
    else:
        pass

    logger = logging.getLogger("ncbi-genome-download")
    logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

    attempt = 0

    if os.path.isfile(infile):
        state = download_acc(args,database)
    elif '/' not in infile and not os.path.isdir(infile):
        state = download_spe(args,args.database)

    while state == 75 and attempt < args.retries:
        attempt += 1
        logger.error(f"Download failed due to connection error. Retries so far {attempt}")
        if os.path.isfile(infile):
            state = download_acc(args,database)
        elif '/' not in infile and not os.path.isdir(infile):
            state = download_spe(args,args.database)

    return state


if __name__ == '__main__':
    main()
    