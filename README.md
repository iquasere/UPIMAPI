# UniProt Id Mapping through API

A tool for retrieving huge ammounts of information from UniProt! 
UPIMAPI is a command line interface for using UniProt's API, which allows to access [UniProt's ID mapping](https://www.uniprot.org/uploadlists/) programmatically!
UPIMAPI was developed as part of the [MOSCA](https://github.com/iquasere/MOSCA) pipeline. It is best used when having a big number of UniProt IDs (like, millions) for which information is required.

## Setting up UPIMAPI

UPIMAPI files can be retrieved from git.
```
git clone https://github.com/iquasere/UPIMAPI.git
```

UPIMAPI requires several packages already installed. They can all be installed with pip.
```
cd UPIMAPI
pip install requirements.txt
```

## Base arguments for running UPIMAPI

UPIMAPI takes either a list of UniProt IDs (one per line) or a BLAST file as input.
```
python upimapi.py -i ids.txt -o uniprotinfo.tsv
python upimapi.py -i aligned.blast --blast -o uniprotinfo.tsv
```

## Additional parameters

```
usage: upimapi.py [-h] -i INPUT [-o OUTPUT] [--excel]
                  [-anncols ANNOTATION_COLUMNS] [-anndbs ANNOTATION_DATABASES]
                  [--blast] [--entry_name] [--fasta]

UniProt Id Mapping through API

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input filename - can be a list of IDs (one per line)
                        or a BLAST TSV file - if so, specify with the --blast
                        parameter
  -o OUTPUT, --output OUTPUT
                        filename of output
  --excel               Will produce output in EXCEL format (default is TSV)
  -anncols ANNOTATION_COLUMNS, --annotation-columns ANNOTATION_COLUMNS
                        List of UniProt columns to obtain information from
  -anndbs ANNOTATION_DATABASES, --annotation-databases ANNOTATION_DATABASES
                        List of databases to cross-check with UniProt
                        information
  --blast               If input file is in BLAST TSV format
  --entry_name          If IDs are in 'Entry name' format: tr|XXX|XXX
  --fasta               Output will be generated in FASTA format

A tool for retrieving information from UniProt.
```