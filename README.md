# UniProt Id Mapping through API

A tool for retrieving huge ammounts of information from UniProt! 
UPIMAPI is a command line interface for using UniProt's API, which allows to access [UniProt's ID mapping](https://www.uniprot.org/uploadlists/) programmatically!
UPIMAPI can handle big numbers of UniProt IDs (like, millions) for which information can be obtained in a single command.
UPIMAPI also allows to first perform annotation with DIAMOND, connecting its powerfull annotation with the convenience of directly obtaining information from UniProt.

### Index

1. [Installing UPIMAPI](https://github.com/iquasere/UPIMAPI#installing-upimapi)
2. [Annotation with UPIMAPI](https://github.com/iquasere/UPIMAPI#annotation-with-upimapi)
3. [Information retrieval from UniProt](https://github.com/iquasere/UPIMAPI#information-retrieval-from-uniprot)
4. [Output](https://github.com/iquasere/UPIMAPI#output)
5. [Additional parameters](https://github.com/iquasere/UPIMAPI#additional-parameters)
6. [Referencing UPIMAPI](https://github.com/iquasere/UPIMAPI#referencing-upimapi)

## Installing UPIMAPI

To install UPIMAPI through Bioconda, run
```
conda install -c conda-forge -c bioconda upimapi
```
To check if UPIMAPI was installed correctly, run
```
upimapi --version
```

## Annotation with UPIMAPI

UPIMAPI can be used to perform homology-based annotation with DIAMOND. Main advantages of using UPIMAPI are that it determines optimal values for the most important search parameters, and directly links annotation to UniProt ID mapping.
To annotate protein sequences and get information from UniProt, UPIMAPI can be run as
```
upimapi -i path/to/sequences.fasta -o path/to/output_directory -db database -t threads
```
where:
* ```sequences.fasta``` is a FASTA file with aminoacid sequences of query proteins
* ```output_directory``` can be any folder, existent or not
* ```database``` can be either "uniprot" (default), "swissprot", "taxids" or the filename of a FASTA file with the reference sequences (see below).

### Reference database

Several points to take notice about the reference database:
* It must be either UniProt or a subsection of it (e.g. SwissProt, or all proteins of a specific taxon). UPIMAPI performs ID mapping with UniProt IDs, so the database must have those;
* It can be supplied in either FASTA (.fasta) or DIAMOND (.dmnd) format. If in FASTA, UPIMAPI will create a new database in DIAMOND format for annotation;
* There are four different ways to input reference databases to UPIMAPI:

#### Use the entire UniProt (or just SwissProt)

Using the UniProt database is a valid choice if the case study is a metagenome with a mostly unknown community composition.

To use the entire UniProt database as reference for UPIMAPI, specify the database as ```--database uniprot```.

If alternatively you only want to use SwissProt (the manually curated part of UniProt), specify the database as ```--database swissprot```.

#### Input tax IDs to build a more specific database

If, for both pure and mixed cultures, the taxonomic composition is known, UPIMAPI can build a database with the reference proteomes of the known taxa. 

To build a reference for specific taxa, specify the database as ```--database taxids```, and the tax IDs as ```--tax-ids taxid1 taxid2 taxid3 ...```.

#### Input a custom database

A custom database can be inputted if, for example, there is only interest in annotating proteins of a specific family (e.g. hydrogenases). Such a database must be manually built from UniProt.

To input a custom database into UPIMAPI, specify it as ```--database path/to/database.fasta```.

## Information retrieval from UniProt

### Columns of information from UniProt

UniProt provides information for many different fields of information and cross-references. For the user's convenience, a default selection is provided: ```Entry```, ```Entry name```, ```Gene names```, ```Protein names```, ```EC number```, ```Function[CC]```, ```Pathway```, ```Keywords```, ```Protein existence```, ```Gene ontology (GO)```, ```Protein families```, ```Taxonomic lineage```, ```Organism```, ```Organism ID```, ```BioCyc```, ```BRENDA```, ```CDD```, ```eggNOG```, ```Ensembl```, ```InterPro```, ```KEGG```, ```Pfam```, ```Reactome```, ```RefSeq``` and ```UniPathway```

If another selection of columns/databases is desired, it can be specified, for example, as 
```
--columns "Coiled coil&Compositional bias"
```
where ```--columns``` takes as input the names of the fields of information required. The complete list of fields available can be consulted at [UniProtKB return fields](https://www.uniprot.org/help/return_fields).

#### UPIMAPI offers a few additional columns for taxonomic information

Previous to the Summer 2022 UniProt release, the API provided fields for taxonomic information, but these have been condensed into the ```Taxonomic lineage``` and ```Taxonomic lineage (IDs)``` columns. Since ```1.8.6```, UPIMAPI provides this information again, properly organized. Additional available columns for taxonomy are as follows:

* ```Taxonomic lineage (LEVEL OF TAXONOMY)```: the taxonomic lineage of the organism, with the specified level of taxonomy. For example, ```--columns "Taxonomic lineage (SPECIES)"``` will return the species of the organism. Other possible values are ```SUPERKINGDOM```, ```PHYLUM```, ```CLASS```, ```ORDER```, ```FAMILY```, ```GENUS```, ```SPECIES```, [among others](https://en.wikipedia.org/wiki/Taxonomic_rank).

* ```Taxonomic lineage IDs (LEVEL OF TAXONOMY)```: the TaxIDs of the organism, with the specified level of taxonomy. For example, ```--columns "Taxonomic lineage IDs (SPECIES)"``` will return the TaxID of the species of the organism. Other possible values are as above.

## ID mapping without annotation

If only retrieval of information from UniProt is required (no annotation step), IDs can be inputted to UPIMAPI directly through several different inputs.

### Annotation BLAST file

The result of an annotation with some database with UniProt IDs can be directy inputted for ID mapping with the command
```
upimapi -i aligned.blast -o output_directory --blast
```

### CSV file

A CSV file with UniProt IDs (separated by commas) can be inputted to UPIMAPI with the command
```
upimapi -i ids.txt -o output_directory
```
This repo provides an [example](https://github.com/iquasere/UPIMAPI/blob/master/ids.txt) of this file.

### Directly from the command line

IDs can also be directly inputted through the command line by not specifying an input. They must be inputted as a comma separated value:
```
>>> upimapi -o output_directory

IDs to perform mapping on (comma separated values):
```

## Output

Information obtained with UPIMAPI can come in two forms:
1. The **Base** (default) workflow obtains information for the list of columns and databases inputted. It produces the following outputs, in the output folder:
    * ```uniprotinfo.tsv```, contains information of the columns and databases specified
    * if annotation was performed, ```aligned.blast``` and ```unaligned.fasta``` contain the annotated and unannotated proteins, respectively.

2. The **Fasta** workflow, specified with the ```--fasta``` argument, results in a FASTA file with the protein sequences correspondent to the inputted IDs

## From/To ID mapping

The ID mapping available at https://www.uniprot.org/id-mapping triggered when "From database" and "To database" are different to the default values - "UniProtKB AC/ID" and "UniProtKB" - is also implemented since UPIMAPI `1.12`.

As an example, this command would convert IDs from UniProtKB to EMBL/Genbank/DDBJ CDS: 
```
upimapi -i ids.txt -o output_directory --from-db 'UniProtKB AC/ID' --to-db 'EMBL/GenBank/DDBJ CDS'
```

Possible values for parameters `--from-db` and `--to-db` can be consulted through the browser (https://www.uniprot.org/id-mapping), at https://rest.uniprot.org/configure/idmapping/fields, or by inputting a wrong value to one of those parameters. Possible options will show up.

This new ID mapping can't be combined with the ID mapping that obtains columns of information from UniProt. UPIMAPI will exit after ID mapping.

## Additional parameters

```
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input filename - can be: 1. a file containing a list of IDs (comma-separated values, no spaces) 2. a BLAST TSV result file (requires to be specified with the
                        --blast parameter 3. a protein FASTA file to be annotated (requires the -db parameter) 4. nothing! If so, will read input from command line, and parse as CSV
                        (id1,id2,...)
  -o OUTPUT, --output OUTPUT
                        Folder to store outputs
  -ot OUTPUT_TABLE, --output-table OUTPUT_TABLE
                        Filename of table output, where UniProt info is stored. If set, will override 'output' parameter just for that specific file
  -rd RESOURCES_DIRECTORY, --resources-directory RESOURCES_DIRECTORY
                        Directory to store resources of UPIMAPI [~/upimapi_resources]
  -cols COLUMNS, --columns COLUMNS
                        List of UniProt columns to obtain information from (separated by &)
  --blast               If input file is in BLAST TSV format (will consider one ID per line if not set) [false]
  --full-id FULL_ID     If IDs in database are in 'full' format: tr|XXX|XXX [auto]
  --fasta               Output will be generated in FASTA format [false]
  --step STEP           How many IDs to submit per request to the API [1000]
  --max-tries MAX_TRIES
                        How many times to try obtaining information from UniProt before giving up [3]
  --sleep SLEEP         Time between requests (in seconds) [3]
  --no-annotation       Do not perform annotation - input must be in one of BLAST result or TXT IDs file or STDIN [false]
  --local-id-mapping    Perform local ID mapping of SwissProt IDs. Advisable if many IDs of SwissProt are present [false]
  --skip-id-mapping     If true, UPIMAPI will not perform ID mapping [false]
  --skip-id-checking    If true, UPIMAPI will not check if IDs are valid before mapping [false]
  --skip-db-check       So UPIMAPI doesn't check for (FASTA) database existence [false]
  --mirror {expasy,uniprot,ebi}
                        From where to download UniProt database [expasy]
  -v, --version         show program's version number and exit

DIAMOND arguments:
  -db DATABASE, --database DATABASE
                        How the reference database is inputted to UPIMAPI. 1. uniprot - UPIMAPI will download the entire UniProt and use it as reference 2. swissprot - UPIMAPI will
                        download SwissProt and use it as reference 3. taxids - Reference proteomes will be downloaded for the taxa specified with the --taxids, and those will be used as
                        reference 4. a custom database - Input will be considered as the database, and will be used as reference
  -t THREADS, --threads THREADS
                        Number of threads to use in annotation steps [all available]
  --evalue EVALUE       Maximum e-value to report annotations for [1e-3]
  --pident PIDENT       Minimum pident to report annotations for.
  --bitscore BITSCORE   Minimum bit score to report annotations for (overrides e-value).
  -mts MAX_TARGET_SEQS, --max-target-seqs MAX_TARGET_SEQS
                        Number of annotations to output per sequence inputed [1]
  -b BLOCK_SIZE, --block-size BLOCK_SIZE
                        Billions of sequence letters to be processed at a time [memory / 20]
  -c INDEX_CHUNKS, --index-chunks INDEX_CHUNKS
                        Number of chunks for processing the seed index [dependant on block size]
  --max-memory MAX_MEMORY
                        Maximum memory to use (in Gb) [all available]
  --taxids TAXIDS       Tax IDs to obtain protein sequences of for building a reference database.
  --diamond-mode {fast,mid_sensitive,sensitive,more_sensitive,very_sensitive,ultra_sensitive}
                        Mode to run DIAMOND with [fast]
```

## Referencing UPIMAPI

If you use UPIMAPI, please cite its [publication](https://www.sciencedirect.com/science/article/pii/S2001037022001179).