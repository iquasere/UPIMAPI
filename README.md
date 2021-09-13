# UniProt Id Mapping through API

A tool for retrieving huge ammounts of information from UniProt! 
UPIMAPI is a command line interface for using UniProt's API, which allows to access [UniProt's ID mapping](https://www.uniprot.org/uploadlists/) programmatically!
UPIMAPI can handle big numbers of UniProt IDs (like, millions) for which information can be obtained in a single command.
UPIMAPI also allows to first perform annotation with DIAMOND, connecting its powerfull annotation with the convenience of directly obtaining information from UniProt.

* [Installing UPIMAPI](https://github.com/iquasere/UPIMAPI#installing-upimapi)
* [Annotation with UPIMAPI](https://github.com/iquasere/UPIMAPI#annotation-with-upimapi)
* [Columns and databases of information from UniProt](https://github.com/iquasere/UPIMAPI#columns-and-databases-of-information-from-uniprot)
* [ID mapping without annotation](https://github.com/iquasere/UPIMAPI#id-mapping-without-annotation)
* [Additional parameters](https://github.com/iquasere/UPIMAPI#additional-parameters)
* [Output](https://github.com/iquasere/UPIMAPI#output)

## Installing UPIMAPI

To install UPIMAPI through Bioconda, run
```
conda install -c conda-forge -c bioconda upimapi
```
To check if it was installed correctly, run
```
upimapi.py --version
```

## Annotation with UPIMAPI

UPIMAPI can be used to perform homology-based annotation with DIAMOND. Main advantages of using UPIMAPI are that it determines optimal values for the most important search parameters, and directly links annotation to UniProt ID mapping.
To annotate protein sequences and get information from UniProt, UPIMAPI can be run as
```
upimapi.py -i sequences.fasta -o output_directory --use-diamond -db database.fasta -t threads
```
where:
* ```sequences.fasta``` is a FASTA file with aminoacid sequences of query proteins
* ```output_directory``` can be any folder, existent or not, and
* ```database.fasta``` is the filename of the reference database

#### Reference database

Several points to take notice about the reference database:
* It must be either UniProt or a subsection of it (e.g. SwissProt, or all proteins of a specific taxon). UPIMAPI performs ID mapping with UniProt IDs, so the database must have those;
* It can be supplied in either FASTA (.fasta) or DIAMOND (.dmnd) format. If in FASTA, UPIMAPI will create a new database in DIAMOND format for annotation;
* There are four different ways to input reference databases to UPIMAPI:
    ##### 1. Use the entire UniProt (or just SwissProt)
    Using the UniProt database is a valid choice if the case study is a metagenome with a mostly unknown community composition.
    
    To use the entire UniProt database as reference for UPIMAPI, specify the database as ```--database uniprot```.
    
    If alternatively you only want to use SwissProt (the manually curated part of UniProt), specify the database as ```--database swissprot```.
    ##### 2. Input tax IDs to build a more specific database
    If, for both pure and mixed cultures, the taxonomic composition is known, UPIMAPI can build a database with the reference proteomes of the known taxa. 
    
    To build a reference for specific taxa, specify the database as ```--database taxids```, and the tax IDs as ```--tax-ids taxid1 taxid2 taxid3 ...```.
    ##### 3. Input a custom database
    A custom database can be inputted if, for example, there is only interest in annotating proteins of a specific family (e.g. hydrogenases). Such a database must be manually built from UniProt.
    
    To input a custom database into UPIMAPI, specify it as ```--database path/to/database.fasta```.

## Columns and databases of information from UniProt

UniProt provides information for many different fields of information and cross-references. For the user's convenience, a default selection is provided:
* **default columns:** ```Entry```, ```Entry name```, ```Gene names```, ```Protein names```, ```EC number```, ```Function[CC]```, ```Pathway```, ```Keywords```, ```Protein existence```, ```Gene ontology (GO)```, ```Protein families```, ```Taxonomic lineage (SUPERKINGDOM)```, ```Taxonomic lineage (PHYLUM)```, ```Taxonomic lineage (CLASS)```, ```Taxonomic lineage (ORDER)```, ```Taxonomic lineage (FAMILY)```, ```Taxonomic lineage (GENUS)```, ```Taxonomic lineage (SPECIES)```

* **default databases:** ```BioCyc Collection of Pathway/Genome Databases```, ```BRENDA Comprehensive Enzyme Information System```, ```Conserved Domains Database```, ```evolutionary genealogy of genes: Non-supervised Orthologous Groups```, ```Ensembl eukaryotic genome annotation project```, ```Integrated resource of protein families, domains and functional sites```, ```KEGG: Kyoto Encyclopedia of Genes and Genomes```, ```KEGG Orthology (KO)```, ```Pfam protein domain database```, ```Reactome - a knowledgebase of biological pathways and processes```, ```NCBI Reference Sequences```, ```UniPathway: a resource for the exploration and annotation of metabolic pathways```

If another selection of columns/databases is desired, it can be specified, for example, as 
```
--columns Coiled coil&Compositional bias --databases UCSC genome browser&Xenopus laevis and tropicalis biology and genomics resource
```
where ```--columns``` and ```--databases``` take as input the names of the [columns](https://www.uniprot.org/help/uniprotkb_column_names) and [databases](https://www.uniprot.org/docs/dbxref) required, respectively. The links provided take to the pages of UniProt where the possible columns and databases values are listed.

## ID mapping without annotation

If only retrieval of information from UniProt is required (no annotation step), IDs can be inputted to UPIMAPI directly through several different inputs.

#### Annotation BLAST file

The result of an annotation with some database with UniProt IDs can be directy inputted for ID mapping with the command
```
upimapi.py -i aligned.blast -o output_directory --blast
```

#### CSV file

A CSV file with UniProt IDs (separated by commas) can be inputted to UPIMAPI with the command
```
upimapi.py -i ids.txt -o output_directory
```
This repo provides an [example](https://github.com/iquasere/UPIMAPI/blob/master/ids.txt) of this file.

#### Directly from the command line

IDs can also be directly inputted through the command line by not specifying an input. They must be inputted as a comma separated value:
```
>>> upimapi.py -o output_directory

IDs to perform mapping on (comma separated values):
```

## Output

Information obtained with UPIMAPI can come in two forms:
1. The **Base** (default) workflow obtains information for the list of columns and databases inputted. It produces the following outputs, in the output folder:
    * ```uniprotinfo.tsv```, contains information of the columns and databases specified
    * if annotation was performed, ```aligned.blast``` and ```unaligned.fasta``` contain the annotated and unannotated proteins, respectively.

2. The **Fasta** workflow, specified with the ```--fasta``` argument, results in a FASTA file with the protein sequences correspondent to the inputted IDs

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
  -ot OUTPUT_TABLE, --output-table OUTPUT_TABLE
                        Filename of table output, where UniProt info is
                        stored. If set, will override 'output' parameter just
                        for that specific file
  -cols COLUMNS, --columns COLUMNS
                        List of UniProt columns to obtain information from
                        (separated by &)
  -dbs ANNOTATION_DATABASES, --annotation-databases ANNOTATION_DATABASES
                        List of databases to cross-check with UniProt
                        information (separated by &)
  --blast               If input file is in BLAST TSV format (will consider
                        one ID per line if not set)
  --full-id FULL_ID     If IDs in database are in 'full' format: tr|XXX|XXX
  --fasta               Output will be generated in FASTA format
  --step STEP           How many IDs to submit per request to the API [1000]
  --max-tries MAX_TRIES
                        How many times to try obtaining information from
                        UniProt before giving up
  --sleep SLEEP         Time between requests (in seconds) [10]
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
                        Number of threads to use in annotation steps [total -
                        2]
  --evalue EVALUE       Maximum e-value to report annotations for [1e-3]
  --pident PIDENT       Minimum pident to report annotations for.
  --bitscore BITSCORE   Minimum bit score to report annotations for (overrides
                        e-value).
  -mts MAX_TARGET_SEQS, --max-target-seqs MAX_TARGET_SEQS
                        Number of annotations to output per sequence inputed
                        [1]
  -b BLOCK_SIZE, --block-size BLOCK_SIZE
                        Billions of sequence letters to be processed at a time
                        (default: auto determine best value)
  -c INDEX_CHUNKS, --index-chunks INDEX_CHUNKS
                        Number of chunks for processing the seed index
                        (default: auto determine best value)
```