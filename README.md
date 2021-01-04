# UniProt Id Mapping through API

A tool for retrieving huge ammounts of information from UniProt! 
UPIMAPI is a command line interface for using UniProt's API, which allows to access [UniProt's ID mapping](https://www.uniprot.org/uploadlists/) programmatically!
UPIMAPI can handle big numbers of UniProt IDs (like, millions) for which information can be obtained in a single command.
UPIMAPI also allows to first perform annotation with DIAMOND, connecting its powerfull annotation with the convenience of directly obtaining information from UniProt.

## Installing UPIMAPI

To install UPIMAPI through Bioconda, run
```
conda install -c conda-forge -c bioconda upimapi
```
To check if it was installed correctly, run
```
upimapi.py --version
```

## Running UPIMAPI

UPIMAPI can be used to first perform annotation with DIAMOND, or directly inputing UniProt IDs to it.

### Running DIAMOND first

To run DIAMOND, the argument ```--use-diamond``` must be specified, and the input most come in FASTA protein format (.faa).
The reference database for aligning the query sequences must be specified with the ```--database``` argument. 
Optionally, parameters for ```--threads```, ```--block-size``` and ```--index-chunks``` can also be specified to speed DIAMOND annotation. 
If not specified, UPIMAPI will automatically determine best values for each of them.

### Outputs

Information obtained with UPIMAPI can come in two forms:
* the base (default) workflow obtains information for the list of columns and databases inputted
* the "fasta" workflow, specified with the ```--fasta``` argument, results in a FASTA file with the protein sequences correspondent to the inputted IDs

## Base arguments for running UPIMAPI

UPIMAPI takes either a list of UniProt IDs (one per line) or a BLAST file as input.
```
python upimapi.py -i ids.txt -o uniprotinfo.tsv
python upimapi.py -i aligned.blast --blast -o uniprotinfo.tsv
```

## Additional parameters

```
usage: upimapi.py [-h] [-i INPUT] [-o OUTPUT] [--excel]
                  [-anncols ANNOTATION_COLUMNS] [-anndbs ANNOTATION_DATABASES]
                  [--blast] [--full-id] [--fasta] [--step STEP] [-v]

UniProt Id Mapping through API

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input filename - can be a list of IDs (one per line)
                        or a BLAST TSV file - if so, specify with the --blast
                        parameter. If no file is given as input, will read
                        from command line input
  -o OUTPUT, --output OUTPUT
                        filename of output
  --excel               Will produce output in EXCEL format (default is TSV)
  -anncols ANNOTATION_COLUMNS, --annotation-columns ANNOTATION_COLUMNS
                        List of UniProt columns to obtain information from
  -anndbs ANNOTATION_DATABASES, --annotation-databases ANNOTATION_DATABASES
                        List of databases to cross-check with UniProt
                        information
  --blast               If input file is in BLAST TSV format
  --full-id             If IDs are in 'full' format: tr|XXX|XXX
  --fasta               Output will be generated in FASTA format
  --step STEP           How many IDs to submit per request to the API (default
                        is 1000)
  -v, --version         show program's version number and exit

A tool for retrieving information from UniProt.
```

#