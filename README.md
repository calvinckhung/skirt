# SKIRT

## Description
This is a project for a pipeline of high-resolution KIR allele annotations on human genome assemblies.

It uses several libraries including argparse, pandas, re, gzip, and Biopython.
The script provides options for outputting in different formats, including BED, VCF, and HAP(CSV).

## Installation

This project requires Python 3.6+ and the following Python libraries installed:

- BioPython
- pandas
- argparse

You can install these packages using pip:
```bash
pip install biopython pandas argparse
```

To install this project, clone the SKIRT repository to your local machine:
```bash
git clone https://github.com/calvinckhung/skirt.git
```

Also, clone the IPDKIR repository to the SKIRT directory:
```bash
cd skirt
git clone https://github.com/ANHIG/IPDKIR.git
```

Don't forget to install tblastn, part of the BLAST+ suite:
```bash
cd ~/bin
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.14.0/ncbi-blast-2.14.0+-x64-linux.tar.gz
tar -xzf ncbi-blast-2.14.0+-x64-linux.tar.gz
echo "export PATH=$PATH:~/bin" >> ~/.bashrc
```
For detail installation description of BLAST+, please refer to https://www.ncbi.nlm.nih.gov/books/NBK569861/

