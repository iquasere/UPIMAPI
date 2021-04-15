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

### Example commands

Get (default) columns and databases information from a list of ids (one per line): will produce uniprotinfo.tsv in ```output_directory```
```
upimapi.py -i ids.txt -o output_directory
```
Get same information for an annotation result, where IDs are in full form (tr|XXX|XXX)
```
upimapi.py -i aligned.blast -o output_directory --blast --full-id
```
Get FASTA sequences from a list of ids: will produce ```output_directory/uniprotinfo.fasta```
```
upimapi.py -i ids.txt --fasta -o output_directory
```
Annotate FASTA protein sequences (.faa) and get information in EXCEL format. Will produce, in ```output_directory```: 
1. uniprotinfo.xlsx - UniProt information in EXCEL format  
2. aligned.blast - BLAST result, output format 6
3. unaligned.blast - non annotated sequences, in FASTA format
```
upimapi.py -i sequences.fasta -o output_directory --use-diamond -db path/to/database.fasta
```

## Additional parameters

```
optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input filename - can be: 1. a file containing a list
                        of IDs (one per line) 2. a BLAST TSV result file
                        (requires to be specified with the --blast parameter
                        3. a protein FASTA file to be annotated (requires the
                        --use-diamond and -db parameters) 4. nothing! If so,
                        will read input from command line, and parse as CSV
                        (id1,id2,...)
  -o OUTPUT, --output OUTPUT
                        Folder to store outputs
  -anncols ANNOTATION_COLUMNS, --annotation-columns ANNOTATION_COLUMNS
                        List of UniProt columns to obtain information from
                        (separated by &)
  -anndbs ANNOTATION_DATABASES, --annotation-databases ANNOTATION_DATABASES
                        List of databases to cross-check with UniProt
                        information (separated by &)
  --blast               If input file is in BLAST TSV format (will consider
                        one ID per line if not set)
  --full-id             If IDs in database are in 'full' format: tr|XXX|XXX
  --fasta               Output will be generated in FASTA format
  --step STEP           How many IDs to submit per request to the API (default
                        is 1000)
  --max-tries MAX_TRIES
                        How many times to try obtaining information from
                        UniProt before giving up
  -v, --version         show program's version number and exit

DIAMOND arguments:
  --use-diamond         Use DIAMOND to annotate sequences before mapping IDs.
                        Requires protein FASTA files as input for "-db" and
                        "-i" parameters
  -db DATABASE, --database DATABASE
                        Reference database for annotation with DIAMOND.
                        NOTICE: if database's IDs are in 'full' format
                        (tr|XXX|XXX), specify with ""--full-id" parameter.
  -t THREADS, --threads THREADS
                        Number of threads to use in annotation steps
  -mts MAX_TARGET_SEQS, --max-target-seqs MAX_TARGET_SEQS
                        Number of annotations to output per sequence inputed
  -b BLOCK_SIZE, --block-size BLOCK_SIZE
                        Billions of sequence letters to be processed at a time
                        (UPIMAPI determines best value for this parameter if
                        not set
  -c INDEX_CHUNKS, --index-chunks INDEX_CHUNKS
                        Number of chunks for processing the seed index
                        (UPIMAPI determines best value for this parameter if
                        not set
```