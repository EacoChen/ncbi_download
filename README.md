# NCBI Mixed Accession Download

Downloading mixed up NCBI assembly/protein/nucleotide accession number,
or search for relative term, and finish the download part

The assembly download script frame is from [Dr. Blin, Kai ncbi-genome-download](https://github.com/kblin/ncbi-genome-download).
The protein and nucleotide part is mainly from [Leon Kuchenbecke EntrezDownloader](https://github.com/lkuchenb/EntrezDownloader) efetch part and mainly frame. And modified by [Tianhua Liao](https://github.com/444thLiao) who write the esearch esummary elink part.

This script is focus on a MIXED Accession Number, which may come from difference literature.
This script is also focus on download some genes or bioproject sequences.

This repository contains:

1. [sequences_tag.txt](sequences_tag.txt) for ncbi sequence searching tag. The Usage can be found below.
2. A test directory [test](test/) which contain a test accession ids, and step by step [test.py](test/test.py)
3. A program directory [ncbi_accession_download](ncbi_accession_download/)


## Table of Contents

- [Install](#install)
- [Usage](#usage)
	- [Accession List Download](#AccessionListDownload)
    - [Searching term and Download ALL](#SearchAndDownloadALL)
    - [BioProject Download](#BioProjectDownload)
- [Function coming soon](#FunctionCommingSoon)
- [License](#license)


## Install

This script is written by python 3.7. And test in >= python 3.6 successfully.
Make sure you have python >= 3.6

Those package will be installed.

'requests >= 2.24.0',
'tqdm >= 4.46.1',
'biopython >= 1.77',
'beautifulsoup4 >= 4.9.1'

Others are Standard Library in python

If you meet some question when in higher version, please let me know.

You can clone this repository from github and go to directory
```
git clone https://github.com/kblin/ncbi-genome-download.git
cd ncbi_download
```

And you can install with pip
```
pip install .
```

or install through python setup.py
```
python setup.py install
```
if you failed with old python version, try updating pip first\
```
pip install --upgrade pip
```

## Usage

This project have to different functions.
1. accession list download
2. search for term and download all accession about it
3. BioProject or BioSample download from literature


```
  -h, --help                    Show this help message and exit
  -i INPUT, --input INPUT
                                List file with accession number seperated by \n.
  -o OUTPUT, --output OUTPUT
                                Downloaded files saved directory.
  -e EMAIL, --email EMAIL
                                An email address. You might get blocked by the NCBI without specifying one.
  -a API, --api API             An API key. You can obtain one by creating an NCBI account. Speeds things up.
  -p PARALLEL, --parallel PARALLEL
                                The number of parallel requests to make
  -b BATCH, --batch BATCH
                                The number of IDs to fetch per request
  -P PROGRESS_BAR, --progress-bar PROGRESS_BAR
                                Enables a progress bar, requires tqdm package
  -u UID, --uid UID             When input file is a uid list.
  -d DATABASE, --database DATABASE
                                When uid is given or input is a str, --databse should be 'assembly'/'protein'/'nucleotide',
  --retmax RETMAX               Searching mode max retrieving sequences, default=[batch*2], maxmium 100,000
  --dry                         Download part not running, just find the uid.
```


### Accession_list_download
You can use test/ data to find if the script worked.

```
ncbi-acc-download -i ./test/s2_a_genome_or_and_protein_id.txt -o ./test_out
```

This [test list](test/s2_a_genome_or_and_protein_id.txt) have 175 accession number in total.
Firstly, the script will find the assembly database to avoid too much nucleotide sequences' uids.
And then will search nucleotide and protien database.

Some WGS sequences in nueclotide(nuccore) database, can not get uid in a parallel way, like 'AHKK01003915', 'AHKK01012539'
So the script will set the batch size to 1, search again.

Accessions like 'AHKK01003915' can use efetch in [E-utilities on Unix Command Line](https://www.ncbi.nlm.nih.gov/books/NBK179288/)

```
assembly: 20
nucleotide: 38
protein: 98
lost: 6
```

In the report, the total sequences is 162. There are two possible situation of reduce,
First is redundancy: The searching process will deduplicate the uid.
Second is Record removed: In the batch searching, this record can not be found. Meanwhile in single searching, this record will not return normally, can not be put in the 'lost' key.


### Search_And_Download_ALL

If I want to download all McrA genes' protein sequences.

```
ncbi-acc-download -i "mcra" -d nucleotide -o test_term
```

The log shows:
```
INFO: Searching the term "mcra"
100%|██████████████████████████████████████████████████████████████████| 1/1 [00:00<00:00,  2.62records/s]
INFO: Find 52035 records
WARNING: The retriving record 40 is samller than records, the result maybe incomplete
INFO: nucleotide sequences downloading...
 50%|███████████████████████████████████                               | 20/40 [00:01<00:01, 12.89records/s][UNKNOWN ERROR] No records found in handle
100%|██████████████████████████████████████████████████████████████████| 40/40 [00:08<00:00,  6.46records/s][UNKNOWN ERROR] No records found in handle
100%|██████████████████████████████████████████████████████████████████| 40/40 [00:08<00:00,  4.85records/s]
ERROR: Falut number: 40
['2169261366', '2169261351', '2169261075', '2169261069', '2169261067', '2169258118', '2169214139', '2169212952', '2169070516', '2169070475', '2169070474', '2169070473', '2169070472', '2169070471', '2169070470', '2169070468', '2169070467', '2169070466', '2169070465', '2169070464', '1371864510', '1371864491', '330506319', '116753325', '2169890971', '224021198', '2169831867', '2169831858', '2169831853', '2169831831', '2169802213', '2169792465', '2169652950', '2169467948', '2169395137', '2169278492', '2169278485', '2169269759', '2169269756', '2169261367']
INFO: nucleotide download Failed
```

The results shows '52035' records. The default batch size is 20, so the retrieving size is '40'.
And all this uid have no nucleotide sequences. In fact, they are WGS sequences head.

So We add some tag behind the mcra, like \[GENE\] or \[PROT\]

>The tag is used in [squence search](sequence_tag.txt). 
There also are [pubmed tag, info tag](https://www.ncbi.nlm.nih.gov/books/NBK179288/)

```
ncbi-acc-download -i "mcra[GENE]" -d nucleotide -o test_term --dry

INFO: Searching the term "mcra[GENE]"
100%|██████████████████████████████████████████████████████████████████| 1/1 [00:00<00:00,  6.22records/s]
INFO: Find 43184 records
WARNING: The retrieving record 40 is samller than records, the result maybe incomplete
INFO: DRY RUN. exiting...
```
'--dry' used for not download. But the uid file will save at test_term/tmp/{term}.txt.
Check the uid file, we can find the most result is the same.

```
ncbi-acc-download -i "mcra[PROT]" -d nucleotide -o test_term --dry

INFO: Searching the term "mcra[PROT]"
100%|██████████████████████████████████████████████████████████████████| 1/1 [00:00<00:00,  6.21records/s]
INFO: Find 682 records
WARNING: The retrieving record 40 is samller than records, the result maybe incomplete
INFO: DRY RUN. exiting...
```

The number of records decrease dramatically, but 682 is a reasonable number.

We can also constrain the length of sequence `300:800[SLEN]`, 
combined with Boolean Operators AND, OR and NOT

```
ncbi-acc-download -i "mcra[PROT] AND 300:800[SLEN]" -d nucleotide -o test_term --dry

INFO: Find 663 records
```

We always get `WARNING: The retrieving record 40 is samller than records, the result maybe incomplete`
`--retmax` show be given correctly like

```
ncbi-acc-download -i "mcra[PROT] AND 300:800[SLEN]" -d nucleotide -o test_term --retmax 663

INFO: Searching the term "mcra[PROT] AND 300:800[SLEN]"
100%|██████████████████████████████████████████████████████████████████| 1/1 [00:00<00:00,  1.78records/s]
INFO: Find 663 records
INFO: nucleotide sequences downloading...
100%|██████████████████████████████████████████████████████████████████| 663/663 [00:07<00:00, 94.56records/s]
INFO: No failed
INFO: 663 sequences have saved to test_term/nucleotide.fasta
INFO: nucleotide download Successful
```

The sequences download successful.
Check the download file:

```
$ grep -c '^>' test_term/nucleotide.fasta
663
```

### BioProject_Download

If you want to download BioProject assembly

```
ncbi-acc-download -i "prjna475886[BioProject]" -d assembly -o test_proj --dry

100%|██████████████████████████████████████████████████████████████████| 1/1 [00:00<00:00,  6.39records/s]
INFO: Find 11 records
INFO: DRY RUN. exiting...
```

No warnings, so delete the `--dry`
```
ncbi-acc-download -i "prjna475886[BioProject]" -d assembly -o test_proj

INFO: Searching the term "prjna475886[BioProject]"
100%|██████████████████████████████████████████████████████████████████| 1/1 [00:00<00:00,  4.76records/s]
INFO: Find 11 records
INFO: genome sequences downloading...
INFO: Getting the NCBI ftp address
100%|██████████████████████████████████████████████████████████████████| 11/11 [00:00<00:00, 48.79records/s]
INFO: No failed

Start download. Parallel=10
100%|██████████████████████████████████████████████████████████████████| 11/11 [00:00<00:00, 26.31it/s]
INFO: assembly download Successful
```

## Function_coming_soon
Searching in more database. Like SRA, GENE, PUBMED and so on.

Other functions. Create an issue, to let me know


## License

[MIT](LICENSE) © Honeyu Chen