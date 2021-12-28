import os
from setuptools import setup, find_packages


def read_version():
    for line in open(os.path.join('ncbi_accession_download', '__init__.py'), 'r'):
        if line.startswith('__version__'):
            return line.split('=')[-1].strip().strip("'")


setup(
    name="ncbi_accession_download",
    version=read_version(),
    author="Honeyu Chen",
    author_email="eacochen@163.com",
    description="Fast, Download mixed ncbi accession number in a batch way.",

    url="https://github.com/EacoChen/ncbi_download",
    packages=find_packages(),
    classifiers = [
                  'Development Status :: 3 - Alpha',
                  'Intended Audience :: Science/Research',
                  'Topic :: Scientific/Engineering :: Bio-Informatics',
                  'License :: OSI Approved :: MIT License',
                  'Programming Language :: Python :: 3',
                  'Natural Language :: Chinese (Simplified)',
                  'Operating System :: Microsoft :: Windows :: Windows 10'
                  ],
    
    install_requires = ['requests >= 2.24.0',
                        'tqdm >= 4.46.1',
                        'biopython >= 1.77',
                        'beautifulsoup4 >= 4.9.1'],
    
    python_requires='>=3.7',
    
    entry_points={
        'console_scripts': [
            'ncbi-acc-download=ncbi_accession_download.core:main',
            'nad=ncbi_accession_download.core:main'
        ],
    },

)