# SKIRT

## Description
This is a project for a pipeline of high-resolution KIR allele annotations on human genome assemblies.

It uses several libraries including argparse, pandas, re, gzip, and Biopython.
The script provides options for outputting in different formats, including BED, VCF, and HAP(CSV).

## Library Installation

### This project requires Python 3.6+ and the following Python libraries installed:

- BioPython
- pandas
- argparse

You can install these packages using pip:
```bash
pip install biopython pandas argparse
```

### To install this project, clone the SKIRT repository to your local machine:
```bash
git clone https://github.com/calvinckhung/skirt.git
```
Export your SKIRT directory
```bash
echo "export SKIRT_WD=<Your SKIRT directory>"
```

## Update IPDKIR
We already include the 2.12.0 IPDKIR in the project. 
In case you want to update from scratch, clone the IPDKIR repository to the SKIRT directory:
```bash
cd skirt
git clone https://github.com/ANHIG/IPDKIR.git
```
This will overwrite the included IPDKIR files. Be careful.

## Install Minimap2
For installation of minimap2, refer to its GitHub page https://github.com/lh3/minimap2.
Ensure your $PATH variable contains the path to the minimap2 executable binary. 
Otherwise, export this variable:
```bash
echo "export MM2_BIN_PATH=<Your minimap2 binary path>"
echo "export PATH=$PATH:$MM2_BIN_PATH" >> ~/.bashrc
```

## Install BLAST+
Don't forget to install tblastn, part of the BLAST+ suite:
```bash
cd ~/bin
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.14.0/ncbi-blast-2.14.0+-x64-linux.tar.gz
tar -xzf ncbi-blast-2.14.0+-x64-linux.tar.gz
echo "export BLAST_BIN_PATH=<Path to tblastn executable binary>"
echo "export PATH=$PATH:$BLAST_BIN_PATH" >> ~/.bashrc
```
For a detailed installation description of BLAST+, please refer to https://www.ncbi.nlm.nih.gov/books/NBK569861/
If you have BLAST+ installed already, don't forget to specify your path as above.

## Usage of SKIRT
To run SKIRT KIR allele annotation for your own assembly sequence, run below command under SKIRT directory:
```bash
./scripts/miniskirt.sh <path to the assembly sequence (.fa or .fa.gz)> <output path>
```

You will get *.BED, *.VCF and *.Hap.CSV files. Usually takes between 5~20 minutes.
