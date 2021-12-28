from Bio import Entrez
from ..ncbi_accession_download.EntrezDownloader import EntrezDownloader
from bs4 import BeautifulSoup as Bp

Entrez.email = 'eacochen@163.com'
Entrez.api_key = 'df1607793b0171d9f7fd39585f90a3491007'
handle = Entrez.esearch('nucleotide',term='AHKK01012539 OR AHKK01003915 OR AY324363')
print(Entrez.read(handle))

handle = Entrez.esearch('nucleotide',term='AHKK01012539')
print(Entrez.read(handle))

infile = 'ncbi_download/test/s2_a_genome_or_and_protein_id.txt'
acc_ids = open(infile).read().split('\n')
# delete the last empty line
if not acc_ids[-1]:
    acc_ids = acc_ids[:-1]

edl = EntrezDownloader(
        email='eacochen@163.com',
        api_key='df1607793b0171d9f7fd39585f90a3491007',
        num_threads=10,
        batch_size=20,
        pbar=True
    )

database = ['assembly', 'nucleotide', 'protein']

ids_dict = {}

_db = database[0]
r_search, f_search = edl.esearch(_db, ['AHKK01012539','AHKK01003915'])

lost_ids = []
ids_dict[_db] = ['Head']

print(r_search[0])
# for result in r_search:
#     _soup = Bp(result, 'xml')
#     lost_ids += [_.get_text() for _ in _soup.find_all('PhraseNotFound')]
#     lost_ids += [_.find('Term').get_text().split('[')[0] for _ in _soup.find_all('TermSet')
#                     if int(_.find('Count').get_text()) > 1]
#     ids_dict[_db] += [_.get_text() for _ in _soup.find_all('Id')]
