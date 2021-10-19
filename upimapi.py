#!/usr/bin/env python
"""
UPIMAPI - UniProt Id Mapping through API

By João Sequeira

Mar 2020
"""

from argparse import ArgumentParser, ArgumentTypeError
import os
from time import strftime, gmtime, time, sleep
import urllib.error
import urllib.parse
import urllib.request
from subprocess import run, Popen, PIPE, check_output

import requests
from psutil import virtual_memory
from pathlib import Path
from multiprocessing import cpu_count
from io import StringIO

import pandas as pd
import xml.etree.ElementTree as ET
from tqdm import tqdm
from datetime import datetime
from Bio import SwissProt as SP

from uniprot_support import UniprotSupport

__version__ = '1.6.0'

upmap = UniprotSupport()


def get_arguments():
    parser = ArgumentParser(description="UniProt Id Mapping through API",
                            epilog="A tool for retrieving information from UniProt.")
    parser.add_argument("-i", "--input", help="""Input filename - can be:
                            1. a file containing a list of IDs (one per line)
                            2. a BLAST TSV result file (requires to be specified with the --blast parameter
                            3. a protein FASTA file to be annotated (requires the --use-diamond and -db parameters)
                            4. nothing! If so, will read input from command line, and parse as CSV (id1,id2,...)""")
    parser.add_argument("-o", "--output", help="Folder to store outputs", default="UPIMAPI_output")
    parser.add_argument(
        "-ot", "--output-table",
        help="Filename of table output, where UniProt info is stored. If set, will override 'output' parameter "
             "just for that specific file")
    parser.add_argument(
        "-rd", "--resources-directory", default=os.path.expanduser("~/upimapi_resources"),
        help="Directory to store resources of UPIMAPI [~/upimapi_resources]")
    parser.add_argument(
        "-cols", "--columns", default=None, help="List of UniProt columns to obtain information from (separated by &)")
    parser.add_argument(
        "-dbs", "--databases", default=None,
        help="List of databases to cross-check with UniProt information (separated by &)")
    parser.add_argument(
        "--blast", action="store_true", default=False,
        help="If input file is in BLAST TSV format (will consider one ID per line if not set)")
    parser.add_argument(
        "--full-id", type=str2bool, default="auto", help="If IDs in database are in 'full' format: tr|XXX|XXX")
    parser.add_argument(
        "--fasta", help="Output will be generated in FASTA format", action="store_true", default=False)
    parser.add_argument(
        "--step", type=int, default=1000, help="How many IDs to submit per request to the API [1000]")
    parser.add_argument(
        "--max-tries", default=3, type=int,
        help="How many times to try obtaining information from UniProt before giving up")
    parser.add_argument("--sleep", default=10, type=int, help="Time between requests (in seconds) [10]")
    parser.add_argument(
        "--no-annotation", action="store_true", default=False,
        help="Do not perform annotation - input must be in one of BLAST result or TXT IDs file or STDIN")
    parser.add_argument(
        "--no-local-mapping", action="store_true", default=False,
        help="Do not perform local ID mapping of SwissProt IDs. Advisable if none or few IDs of SwissProt are present.")
    parser.add_argument('-v', '--version', action='version', version=f'UPIMAPI {__version__}')

    diamond_args = parser.add_argument_group('DIAMOND arguments')
    diamond_args.add_argument(
        "-db", "--database", default='uniprot',
        help="How the reference database is inputted to UPIMAPI."
             "1. uniprot - UPIMAPI will download the entire UniProt and use it as reference"
             "2. swissprot - UPIMAPI will download SwissProt and use it as reference"
             "3. taxids - Reference proteomes will be downloaded for the taxa specified with the --taxids, and those "
             "will be used as reference"
             "4. a custom database - Input will be considered as the database, and will be used as reference")
    diamond_args.add_argument(
        "-t", "--threads", type=int, default=cpu_count() - 2,
        help="Number of threads to use in annotation steps [total available - 2]")
    diamond_args.add_argument(
        "--evalue", type=float, default=1e-3, help="Maximum e-value to report annotations for [1e-3]")
    diamond_args.add_argument(
        "--pident", type=float, default=None, help="Minimum pident to report annotations for.")
    diamond_args.add_argument(
        "--bitscore", type=float, default=None, help="Minimum bit score to report annotations for (overrides e-value).")
    diamond_args.add_argument(
        "-mts", "--max-target-seqs", default=1, help="Number of annotations to output per sequence inputed [1]")
    diamond_args.add_argument(
        "-b", "--block-size", type=int,
        help="Billions of sequence letters to be processed at a time (default: auto determine best value)")
    diamond_args.add_argument(
        "-c", "--index-chunks", type=int,
        help="Number of chunks for processing the seed index (default: auto determine best value)")
    diamond_args.add_argument(
        "--taxids", nargs="+", help="Tax IDs to obtain protein sequences of for building a reference database.")

    args = parser.parse_args()
    args.output = args.output.rstrip('/')
    args.resources_directory = args.resources_directory.rstrip('/')

    return args


def timed_message(message):
    print(f'{strftime("%Y-%m-%d %H:%M:%S", gmtime())}: {message}')


def str2bool(v):
    if v.lower() == 'auto':
        return 'auto'
    elif v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ArgumentTypeError('Boolean value expected.')


def get_fasta_ids(filename):
    return [line[1:-1] for line in open(filename) if line.startswith('>')]


def parse_blast(blast):
    result = pd.read_csv(blast, sep='\t', header=None)
    result.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart',
                      'send', 'evalue', 'bitscore']
    return result


def uniprot_request(ids, original_database='ACC+ID', database_destination='',
                    output_format='tab', columns=None, databases=None):
    """
    Input:
        ids: list of UniProt IDs to query
        original_database: database from where the IDs are
        database_destination: database to where to map (so far, only works with 'ACC'
        output_format: format of response to get
        columns: names of UniProt columns to get info on
        databases: names of databases to cross-reference with
    Output:
        Returns the content of the response from UniProt
    """

    base_url = 'https://www.uniprot.org/uploadlists/'

    params = {
        'from': original_database,
        'format': output_format,
        'query': ' '.join(ids),
        'columns': upmap.string4mapping(columns=columns, databases=databases)
    }

    if database_destination == '' or original_database == 'ACC+ID':
        params['to'] = 'ACC'

    data = urllib.parse.urlencode(params).encode("utf-8")
    request = urllib.request.Request(base_url, data)
    response = urllib.request.urlopen(request)
    return response.read().decode("utf-8")


def get_uniprot_information(
        ids, original_database='ACC+ID', database_destination='', step=1000, sleep_time=30, columns=None,
        databases=None, max_tries=3):
    """
    Input:
        ids: list of UniProt IDs to query
        original_database: database from where the IDs are
        database_destination: database to where to map (so far, only works with 'ACC'
        chunk: INT, number of IDs to send per request
        sleep_time: INT, number of seconds to wait between requests
        columns: list - names of UniProt columns to get info on
        databases: list - names of databases to cross-reference with
    Output:
        pd.DataFrame will be returned with the information about the IDs queried.
    """
    result = pd.DataFrame()
    for i in tqdm(range(0, len(ids), step), desc=f'Retrieving UniProt information from {str(len(ids))} IDs'):
        tries = 0
        done = False
        j = min(i + step, len(ids))
        while not done and tries < max_tries:
            try:
                data = uniprot_request(
                    ids[i:j], original_database, database_destination, columns=columns, databases=databases)
                if len(data) > 0:
                    uniprotinfo = pd.read_csv(StringIO(data), sep='\t')
                    k = (len(uniprotinfo.columns) - (30 if (
                            len(columns + databases) == 0 or (columns == [''] and databases == [''])
                    ) else (len(columns) + len(databases))))  # Removes the "yourlist:" and "isomap:" columns
                    uniprotinfo = uniprotinfo[uniprotinfo.columns.tolist()[:-k]]
                    result = pd.concat([result, uniprotinfo[uniprotinfo.columns.tolist()]])
                sleep(sleep_time)
                done = True
            except ConnectionError:
                print(f'ID mapping failed. Remaining tries: {max_tries - tries}')
                tries += 1
                sleep(10)
    return result


def get_uniprot_fasta(ids, step=1000, sleep_time=30):
    """
    Input:
        ids: list of UniProt IDs to query
        chunk: INT, number of IDs to send per request
        sleep_time: INT, number of seconds to wait between requests
    Output:
        str object containing the fasta sequences and headers
        of the proteis belonging to the IDs queried will be returned
    """

    print(f'Building FASTA from {len(ids)} IDs.')
    result = str()
    for i in tqdm(range(0, len(ids), step), desc="UniProt ID mapping"):
        j = min(i + step, len(ids))
        data = uniprot_request(ids[i:j], original_database='ACC+ID', database_destination='', output_format='fasta')
        if len(data) > 0:
            result += data
        sleep(sleep_time)
    return result


def uniprot_fasta_workflow(all_ids, output, max_iter=5, step=1000, sleep_time=10):
    if os.path.isfile(output):
        print(f'{output} was found. Will perform mapping for the remaining IDs.')
        ids_done = get_fasta_ids(output)
    else:
        print(f'{output} not found. Will perform mapping for all IDs.')
        ids_done = list()

    ids_missing = list(set(all_ids) - set(ids_done))

    tries = 0
    ids_done = ([ide.split('|')[1] for ide in get_fasta_ids(output)] if os.path.isfile(output) else list())
    while len(ids_done) < len(all_ids) and tries < max_iter:
        print('Checking which IDs are missing information.')
        ids_missing = list(set([ide for ide in tqdm(all_ids, desc='Checking which IDs are missing information.')
                                if ide not in ids_done]))
        print(f'Information already gathered for {len(ids_done)} ids. Still missing for {len(ids_missing)}.')
        uniprotinfo = get_uniprot_fasta(ids_missing, step=step, sleep_time=sleep_time)
        with open(output, 'a') as file:
            file.write(uniprotinfo)
        ids_done = [ide.split('|')[1] for ide in get_fasta_ids(output)]
        tries += 1
    if len(ids_done) == len(all_ids):
        print(f'Results for all IDs are available at {output}')
    else:
        ids_unmapped_output = f"{'/'.join(output.split('/')[:-1])}/ids_unmapped.txt"
        handler = open(ids_unmapped_output, 'a')
        handler.write('\n'.join(ids_missing))
        print(f'Maximum iterations were made. Results related to {str(len(ids_missing))} IDs were not obtained. '
              f'IDs with missing information are available at {ids_unmapped_output} and information obtained is '
              f'available at {output}')


def uniprot_information_workflow(ids, output, max_iter=5, columns=None, databases=None, step=1000,
                                 sleep_time=10):
    if os.path.isfile(output) and os.stat(output).st_size > 1:
        try:
            result = pd.read_csv(output, sep='\t', low_memory=False).drop_duplicates()
            print(f'{output} was found. Will perform mapping for the remaining IDs.')
            ids_done = list(set(result['Entry'].tolist() + result['Entry name'].tolist()))
        except OSError:
            print(f'{output} was found. However, it could not be parsed. Will restart mapping.')
            result = pd.DataFrame()
            ids_done = list()
    else:
        print(f'{output} not found or empty. Will perform mapping for all IDs.')
        result = pd.DataFrame()
        ids_done = list()
    tries = 0
    ids_unmapped_output = f"{'/'.join(output.split('/')[:-1])}/ids_unmapped.txt"
    ids_missing = list(set(ids) - set(ids_done))
    last_ids_missing = None

    print(f'IDs present in uniprotinfo file: {len(ids_done)}')
    print(f'IDs missing: {len(ids_missing)}')

    while len(ids_missing) > 0 and tries < max_iter and ids_missing != last_ids_missing:
        print(f'Information already gathered for {len(ids_done)} ids. '
              f'Still missing for {len(ids_missing)}.')
        last_ids_missing = ids_missing
        uniprotinfo = get_uniprot_information(
            ids_missing, step=step, columns=columns, databases=databases, max_tries=max_iter, sleep_time=sleep_time)
        if len(uniprotinfo) > 0:
            ids_done += list(set(uniprotinfo['Entry'].tolist() + uniprotinfo['Entry name'].tolist()))
            result = pd.concat([result, uniprotinfo], ignore_index=True)
        ids_missing = list(set(ids) - set(ids_done))

        if len(ids_missing) > 0:
            if last_ids_missing == ids_missing:
                print("Could not map additional IDs for this mapping. There were probably some outdated IDs. "
                      "For more questions, please contact through https://github.com/iquasere/UPIMAPI/issues")
            else:
                print('Failed to retrieve information for some IDs. Retrying request.')
                tries += 1

    result.to_csv(output, sep='\t', index=False)

    if len(ids_missing) == 0:
        print(f'Results for all IDs are available at {output}')
    else:
        open(ids_unmapped_output, 'a').write('\n'.join(ids_missing))
        print(f"Maximum iterations were made. Results related to {str(len(ids_missing))} IDs were not obtained. "
              f"IDs with missing information are available at {ids_unmapped_output} and information obtained is "
              f"available at {output}")


def determine_full_id(ids):
    for ide in ids:
        if '|' in ide:
            return True
    return False


def get_ids(inpute, input_type='blast', full_id='auto'):
    if input_type == 'blast':
        ids = parse_blast(inpute)['sseqid']
    elif input_type == 'txt':
        ids = open(inpute).read().split('\n')
    else:
        ids = inpute.split(',')
    if full_id == 'auto':
        full_id = determine_full_id(ids)
        print(f'Auto determined "full id" as: {full_id}')
    if full_id:
        return [ide.split('|')[1] for ide in ids if ide != '*'], full_id
    return [ide for ide in ids if ide != '*'], full_id


def run_command(bash_command, print_message=True):
    if print_message:
        print(bash_command)
    run(bash_command.split(), check=True)


def run_pipe_command(bash_command, output='', mode='w', print_message=True):
    if print_message:
        print(bash_command)
    if output == '':
        Popen(bash_command, stdin=PIPE, shell=True).communicate()
    elif output == 'PIPE':
        return Popen(bash_command, stdin=PIPE, shell=True,
                     stdout=PIPE).communicate()[0].decode('utf8')
    else:
        with open(output, mode) as output_file:
            Popen(bash_command, stdin=PIPE, shell=True, stdout=output_file).communicate()


def make_diamond_database(fasta, dmnd):
    run_command(f'diamond makedb --in {fasta} -d {dmnd}')


def block_size_and_index_chunks(argsb, argsc):
    if argsb:
        b = argsb
    else:
        b = virtual_memory().available / (1024.0 ** 3) / 20  # b = memory in Gb / 20
    if argsc:
        return b, argsc
    if b > 3:
        return b, 1
    if b > 2:
        return b, 2
    if b > 1:
        return b, 3
    return b, 4


def run_diamond(query, aligned, unaligned, database, threads=12, max_target_seqs=50, b=1, c=4, e_value=0.01,
                bit_score=None, pident=None):
    command = (
        f"diamond blastp --query {query} --out {aligned} --un {unaligned} --db {database} --outfmt 6 --unal 1 "
        f"--threads {threads} --max-target-seqs {max_target_seqs} -b {b} -c {c} --evalue {e_value} --very-sensitive"
    )
    if bit_score:
        command += f' --min-score {bit_score}'
    if pident:
        command += f' --id {pident}'
    run_command(command)


def get_proteome_for_taxid(taxid, max_tries=3):
    tries = 0
    res = None
    done = False
    while tries < max_tries and not done:
        try:
            res = requests.get(f'https://www.uniprot.org/uniprot/?query=taxonomy:{taxid}&format=fasta')
            done = True
        except ConnectionError:
            tries += 1
            sleep(10)
    return res.content.decode('utf8')


def local_uniprot_is_outdated(local_reldate_file):
    local = open(local_reldate_file).readlines()
    [sp_date, tr_date] = [datetime.strptime(local[i][:-1].split()[-1], '%d-%b-%Y') for i in [1, 2]]
    current = requests.get("https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/reldate.txt"
                           ).content.decode('utf8').split('\n')
    [c_sp_date, c_tr_date] = [datetime.strptime(current[i][:-1].split()[-1], '%d-%b-%Y') for i in [1, 2]]
    return c_sp_date < sp_date or c_tr_date < tr_date


def download_with_progress_bar(url, output_folder):
    # Streaming, so we can iterate over the response.
    response = requests.get(url, stream=True)
    total_size_in_bytes = int(response.headers.get('content-length', 0))
    block_size = 1024  # 1 Kibibyte
    progress_bar = tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True, desc=f'Downloading {url.split("/")[-1]}')
    with open(f'{output_folder}/{url.split("/")[-1]}', 'wb') as file:
        for data in response.iter_content(block_size):
            progress_bar.update(len(data))
            file.write(data)
    progress_bar.close()


def download_uniprot(output_folder):
    for url in [
        "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz",
        "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_trembl.fasta.gz",
        "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/reldate.txt"
    ]:
        print(f'Downloading and writing: {url}')
        download_with_progress_bar(url, output_folder)
    run_pipe_command(f'zcat {output_folder}/uniprot_trembl.fasta.gz {output_folder}/uniprot_sprot.fasta.gz > '
                     f'{output_folder}/uniprot.fasta')
    for file in [f'{output_folder}/uniprot_trembl.fasta.gz', f'{output_folder}/uniprot_sprot.fasta.gz']:
        os.remove(file)


def build_reference_database(database, output_folder, taxids=None, max_tries=3):
    if database == 'uniprot':
        download_uniprot(output_folder)
    elif database == 'swissprot':
        download_with_progress_bar(
            "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz",
            output_folder)
        run_command(f'gunzip {output_folder}/uniprot_sprot.fasta.gz')
    elif database == 'taxids':
        for taxid in tqdm(taxids, desc=f'Retrieving reference proteomes for {len(taxids)} taxa from UniProt.'):
            with open(f'{output_folder}/upimapi_database.fasta', 'a') as f:
                f.write(get_proteome_for_taxid(taxid, max_tries=max_tries))


def must_build_database(database, resources_folder):
    db2suffix = {'uniprot': 'uniprot.fasta', 'swissprot': 'uniprot_sprot.fasta', 'taxids': 'upimapi_database.fasta'}
    if database in db2suffix.keys():
        if os.path.isfile(f'{resources_folder}/{db2suffix[database]}'):
            return str2bool(input(f'{resources_folder}/{db2suffix[database]} exists. Overwrite? [Y/N] '))
    return True


def get_tabular_taxonomy(output):
    res = requests.get('https://ftp.uniprot.org/pub/databases/uniprot/current_release/rdf/taxonomy.rdf.xz')
    with open('taxonomy.rdf.xz', 'wb') as f:
        f.write(res.content)
    run_command(f'unxz taxonomy.rdf.xz')
    print('Reading RDF taxonomy')
    root = ET.parse('taxonomy.rdf').getroot()
    elems = root.findall('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}Description')
    with open(output, 'w') as f:
        written = f.write('\t'.join(['taxid', 'name', 'rank', 'parent_taxid']) + '\n')
        for elem in tqdm(elems, desc='Converting XML taxonomy.rdf to TSV format'):
            info = [elem.get('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about').split('/')[-1]]
            scientific_name = elem.find('{http://purl.uniprot.org/core/}scientificName')
            info.append(scientific_name.text if scientific_name is not None else '')
            rank_elem = elem.find('{http://purl.uniprot.org/core/}rank')
            info.append(rank_elem.get('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource').split('/')[-1]
                        if rank_elem is not None else '')
            upper_taxon = elem.find('{http://www.w3.org/2000/01/rdf-schema#}subClassOf')
            info.append(upper_taxon.get('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource').split('/')[-1]
                        if upper_taxon is not None else '')
            written = f.write('\t'.join(info) + '\n')


def get_match_id(record, ids):
    if record.entry_name in ids:
        return record.entry_name
    if record.accessions[0] in ids:
        return record.accessions[0]
    return None


def count_on_file(expression, file, compressed=False):
    return int(check_output(f"{'zgrep' if compressed else 'grep'} -c '{expression}' {file}", shell=True))


def get_local_swissprot_data(sp_dat_filename, ids):
    sp_dat = SP.parse(open(sp_dat_filename))
    result = list()
    i = 1
    record = next(sp_dat)
    number_of_entries = count_on_file('Reviewed;', sp_dat_filename)
    while record is not None and len(ids) > 0:
        match_id = get_match_id(record, ids)
        if match_id is not None:
            result.append(record.__dict__)
            ids.remove(match_id)
        if i % 100000 == 0:
            print(f'[{i}/{number_of_entries}] SwissProt entries queried')
        record = next(sp_dat, None)
        i += 1
    print(f'[{i}/{number_of_entries}] SwissProt entries queried')
    return pd.DataFrame(result), ids


# TODO - someday, so this becomes like "protein names" column of UniProt's API
def description2names(description):
    data = [key_value.split('=') for key_value in description.split('; ')]
    data = {v[0]: v[1] for v in data if len(v) == 2}


def lineage_to_columns(lineage, tax_tsv):
    result = dict()
    for taxon in lineage:
        level = tax_tsv.loc[taxon]["rank"]
        if type(level) == str:
            result[f'Taxonomic lineage ({level.upper()})'] = taxon
        elif type(level) == float:
            continue
        else:   # some taxIDs have multiple levels (e.g. "Craniata")
            for lvl in level:
                if type(lvl) == str:
                    result[f'Taxonomic lineage ({lvl.upper()})'] = taxon
    return pd.Series(result)


def get_taxonomy_in_columns(data, tax_tsv):
    tax_tsv = pd.read_csv(tax_tsv, sep='\t')
    tax_tsv = tax_tsv[tax_tsv.name.notnull()]
    tax_tsv.set_index('name', inplace=True)
    tax_df = pd.DataFrame()
    for lineage in data['organism_classification']:
        tax_df = tax_df.append(lineage_to_columns(lineage, tax_tsv), ignore_index=True)
    return tax_df


def append_value(result, comment):
    for key in result.keys():
        if comment.startswith(key):
            result[key].append(comment)
            return result


def get_comments(comments_data):
    result = {key: [] for key in [
        'FUNCTION:', 'SUBUNIT:', 'INTERACTION:', 'SUBCELLULAR LOCATION:', 'ALTERNATIVE PRODUCTS:',
        'TISSUE SPECIFICITY:', 'PTM:', 'POLYMORPHISM:', 'DISEASE:', 'MISCELLANEOUS:', 'SIMILARITY:', 'CAUTION:',
        'SEQUENCE CAUTION:', 'WEB RESOURCE:']}
    for comments in comments_data:
        for comment in comments:
            append_value(result, comment)
    result = pd.DataFrame.from_dict(result, orient='index').reset_index()
    for col in result.columns.tolist():
        result[col] = result[col].apply('; '.join)
    result.columns = [
        'Function [CC]', 'Subunit structure [CC]', 'Interacts with', 'Subcellular location [CC]',
        'Alternative products (isoforms)', 'Tissue specificity', 'Post-translational modification', 'Polymorphism',
        'Involvement in disease', 'Miscellaneous [CC]', 'Sequence similarities', 'Caution', 'Sequence caution',
        'Web resources']
    return result


def parse_sp_data(sp_data, tax_tsv):
    result = pd.DataFrame()
    result['Entry'] = sp_data['accessions'].apply(lambda x: x[0])
    for k, v in upmap.local2api.items():
        if v not in [None, False]:
            result[v] = sp_data[k]
    result['Taxonomic identifier (SPECIES)'] = sp_data['taxonomy_id'].apply(lambda x: x[0])
    result['Virus hosts'] = sp_data['host_organism'].apply(lambda x: x[0])
    result['Gene names (primary )'] = sp_data['gene_name'].apply(lambda x: x.split('=')[1].split('; ')[0])
    result['Organism'] = sp_data['organism'].str.rstrip('.')
    result = pd.merge(result, get_taxonomy_in_columns(sp_data, tax_tsv), left_index=True, right_index=True, how='left')
    result = pd.merge(result, get_comments(sp_data['comments']), left_index=True, right_index=True, how='left')
    result['Keywords'] = sp_data['keywords'].apply(';'.join)
    result['Date of creation'] = sp_data['created'].apply(
        lambda x: datetime.strptime(x[0], '%d-%b-%Y').strftime('%Y-%m-%d'))
    result['Date of last modification'] = sp_data['annotation_update'].apply(
        lambda x: datetime.strptime(x[0], '%d-%b-%Y').strftime('%Y-%m-%d'))
    result['Version (entry)'] = sp_data['annotation_update'].apply(lambda x: x[1])
    result['Date of last sequence modification'] = sp_data['sequence_update'].apply(
        lambda x: datetime.strptime(x[0], '%d-%b-%Y').strftime('%Y-%m-%d'))
    result['Version (sequence)'] = sp_data['sequence_update'].apply(lambda x: x[1])
    # TODO - work gene names for columns "Gene names" and "Gene names (ORF )"
    return result


def get_sprot_dat(sp_dat):
    r = requests.get(
        'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz',
        allow_redirects=True)
    with open(sp_dat, 'wb') as f:
        f.write(r.content)


def local_id_mapping(ids, sp_dat, tax_tsv, output):
    if not os.path.isfile(sp_dat):
        get_sprot_dat(sp_dat)
    if not os.path.isfile(tax_tsv):
        get_tabular_taxonomy(tax_tsv)
    sp_data, ids_not_found = get_local_swissprot_data(sp_dat, ids)
    parse_sp_data(sp_data, tax_tsv).to_csv(output, sep='\t', index=False)
    return ids_not_found


def upimapi():
    args = get_arguments()
    Path(args.output).mkdir(parents=True, exist_ok=True)
    Path(args.resources_directory).mkdir(parents=True, exist_ok=True)

    # Using annotation with DIAMOND
    if not args.no_annotation:
        db2file = {'uniprot': f'{args.resources_directory}/uniprot.fasta',
                   'swissprot': f'{args.resources_directory}/uniprot_sprot.fasta',
                   'taxids': f'{args.resources_directory}/upimapi_database.fasta'}
        if args.database in db2file.keys():
            database = db2file[args.database]
        else:
            database = args.database

        if must_build_database(args.database, args.resources_directory):
            build_reference_database(
                args.database, args.resources_directory, taxids=args.taxids, max_tries=args.max_tries)
        if not database.endswith(".dmnd"):
            if not os.path.isfile(f"{'.'.join(database.split('.')[:-1])}.dmnd"):
                make_diamond_database(database, f"{'.'.join(database.split('.')[:-1])}.dmnd")
            database = f"{'.'.join(database.split('.')[:-1])}.dmnd"
        (b, c) = block_size_and_index_chunks(argsb=args.block_size, argsc=args.index_chunks)

        run_diamond(
            args.input, f'{args.output}/aligned.blast', f'{args.output}/unaligned.blast',
            database, threads=args.threads, max_target_seqs=args.max_target_seqs, b=b, c=c,
            e_value=args.evalue, bit_score=args.bitscore, pident=args.pident)
        args.input = f'{args.output}/aligned.blast'
        args.blast = True

    # Get IDs from STDIN and set input type
    if args.input is None:
        args.input = input('IDs to perform mapping on (comma separated values):')
        input_type = 'stdin'
    elif args.blast:
        input_type = 'blast'
    else:
        input_type = 'txt'

    # Get the IDs
    ids, full_id = get_ids(args.input, input_type=input_type, full_id=args.full_id)

    # Get UniProt information
    if not args.fasta:
        if args.columns is not None:
            args.columns = args.columns.split('&')
        else:
            args.columns = upmap.get_default_columns()
        if args.databases is not None:
            args.databases = args.databases.split('&')
        else:
            args.databases = upmap.get_default_databases()

        if args.output_table:
            table_output = args.output_table
            print(f'Overrided table output to {table_output}')
            Path('/'.join(args.output_table.split('/')[:-1])).mkdir(parents=True, exist_ok=True)
        else:
            table_output = f'{args.output}/uniprotinfo.tsv'

        # ID mapping through local information
        if not args.no_local_mapping:
            ids = local_id_mapping(
                ids, f'{args.resources_directory}/uniprot_sprot.dat', f'{args.resources_directory}/taxonomy.tsv',
                table_output)

        # ID mapping through API
        uniprot_information_workflow(
            ids, table_output, columns=args.columns, databases=args.databases, step=args.step, max_iter=args.max_tries,
            sleep_time=args.sleep)

        if not args.no_annotation:
            blast = parse_blast(f'{args.output}/aligned.blast')
            if full_id:
                blast.sseqid = [ide.split('|')[1] if ide != '*' else ide for ide in blast.sseqid]
            result = pd.merge(blast, pd.read_csv(table_output, sep='\t'), left_on='sseqid', right_on='Entry')
            result.to_csv(f'{args.output}/UPIMAPI_results.tsv', index=False, sep='\t')
    else:
        uniprot_fasta_workflow(ids, f'{args.output}/uniprotinfo.fasta', step=args.step, sleep_time=args.sleep)


if __name__ == '__main__':
    start_time = time()
    upimapi()
    timed_message(f'UPIMAPI analysis finished in {strftime("%Hh%Mm%Ss", gmtime(time() - start_time))}')
