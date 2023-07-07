#!/usr/bin/env python
"""
UPIMAPI - UniProt Id Mapping through API

By João Sequeira

Mar 2020
"""

import json
from argparse import ArgumentParser, ArgumentTypeError
import os
import sys
from time import strftime, gmtime, time, sleep
from subprocess import run, Popen, PIPE, check_output
import requests
import yaml
from psutil import virtual_memory
from pathlib import Path
from multiprocessing import cpu_count, Pool, Manager
from io import StringIO
import pandas as pd
import xml.etree.ElementTree as ET
from tqdm import tqdm
from datetime import datetime
from Bio import SwissProt as SP
import numpy as np
from functools import partial
import re

__version__ = '1.12.0'


def load_api_info():
    return yaml.safe_load(requests.get('https://rest.uniprot.org/docs/uniprot-openapi3.yaml').text)


def get_url(url, **kwargs):
    response = requests.get(url, **kwargs)
    if not response.ok:
        response.raise_for_status()
        sys.exit(response.txt)
    return response


def get_uniprot_columns():
    res = get_url('https://rest.uniprot.org/configure/uniprotkb/result-fields')
    obj = json.loads(res.text)
    result = {}
    for i in range(len(obj)):
        for col in obj[i]['fields']:
            result[col['label']] = col['name']
    return result


def get_id_mapping_fields():
    res = get_url('https://rest.uniprot.org/configure/idmapping/fields')
    obj = json.loads(res.text)
    froms, tos = {}, {}
    for group in range(len(obj['groups'])):
        for item in obj['groups'][group]['items']:
            if item['from']:
                froms[item['displayName']] = item['name']
            if item['to']:
                tos[item['displayName']] = item['name']
    return froms, tos


api_info = load_api_info()
columns_dict = get_uniprot_columns()
from_fields, to_fields = get_id_mapping_fields()


def get_arguments():
    parser = ArgumentParser(description="UniProt Id Mapping through API",
                            epilog="A tool for retrieving information from UniProt.")
    parser.add_argument(
        "-i", "--input", help="""Input filename - can be:\n
        \t1. a file containing a list of IDs (comma-separated values, no spaces)\n
        \t2. a BLAST TSV result file (requires to be specified with the --blast parameter\n
        \t3. a protein FASTA file to be annotated (requires the -db parameter)\n
        \t4. nothing! If so, will read input from command line, and parse as CSV (id1,id2,...)""")
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
        "--from", default="UniProtKB AC/ID", choices=from_fields.keys(),
        help="Which database are the IDs from. If from UniProt, default is fine [UniProtKB AC/ID]")
    parser.add_argument(
        "--to", default="UniProtKB", choices=to_fields.keys(),
        help="To which database the IDs should be mapped. If only interested in columns information "
             "(which include cross-references), default is fine [UniProtKB]")
    parser.add_argument(
        "--blast", action="store_true", default=False,
        help="If input file is in BLAST TSV format (will consider one ID per line if not set) [false]")
    parser.add_argument(
        "--full-id", type=str2bool, default="auto", help="If IDs in database are in 'full' format: tr|XXX|XXX [auto]")
    parser.add_argument(
        "--fasta", help="Output will be generated in FASTA format [false]", action="store_true", default=False)
    parser.add_argument(
        "--step", type=int, default=1000, help="How many IDs to submit per request to the API [1000]")
    parser.add_argument(
        "--max-tries", default=3, type=int,
        help="How many times to try obtaining information from UniProt before giving up [3]")
    parser.add_argument("--sleep", default=3, type=int, help="Time between requests (in seconds) [3]")
    parser.add_argument(
        "--no-annotation", action="store_true", default=False,
        help="Do not perform annotation - input must be in one of BLAST result or TXT IDs file or STDIN [false]")
    parser.add_argument(
        "--local-id-mapping", action="store_true", default=False,
        help="Perform local ID mapping of SwissProt IDs. Advisable if many IDs of SwissProt are present [false]")
    parser.add_argument(
        "--skip-id-mapping", action="store_true", default=False,
        help="If true, UPIMAPI will not perform ID mapping [false]")
    parser.add_argument(
        "--skip-id-checking", action="store_true", default=False,
        help="If true, UPIMAPI will not check if IDs are valid before mapping [false]")
    parser.add_argument(
        "--skip-db-check", action="store_true", default=False,
        help="So UPIMAPI doesn't check for (FASTA) database existence [false]")
    parser.add_argument(
        "--mirror", choices=['expasy', 'uniprot', 'ebi'], default='expasy',
        help="From where to download UniProt database [expasy]")
    parser.add_argument('-v', '--version', action='version', version=f'UPIMAPI {__version__}')

    diamond_args = parser.add_argument_group('DIAMOND arguments')
    diamond_args.add_argument(
        "-db", "--database", default='uniprot',
        help="How the reference database is inputted to UPIMAPI.\n"
             "\t1. uniprot - UPIMAPI will download the entire UniProt and use it as reference\n"
             "\t2. swissprot - UPIMAPI will download SwissProt and use it as reference\n"
             "\t3. taxids - Reference proteomes will be downloaded for the taxa specified with the --taxids, and those "
             "will be used as reference\n"
             "\t4. a custom database - Input will be considered as the database, and will be used as reference")
    diamond_args.add_argument(
        "-t", "--threads", type=int, default=cpu_count(),
        help="Number of threads to use in annotation steps [all available]")
    diamond_args.add_argument(
        "--evalue", type=float, default=1e-3, help="Maximum e-value to report annotations for [1e-3]")
    diamond_args.add_argument(
        "--pident", type=float, default=None, help="Minimum pident to report annotations for.")
    diamond_args.add_argument(
        "--bitscore", type=float, default=None, help="Minimum bit score to report annotations for (overrides e-value).")
    diamond_args.add_argument(
        "-mts", "--max-target-seqs", type=int, default=1,
        help="Number of annotations to output per sequence inputed [1]")
    diamond_args.add_argument(
        "-b", "--block-size", type=float,
        help="Billions of sequence letters to be processed at a time [memory / 20]")
    diamond_args.add_argument(
        "-c", "--index-chunks", type=int,
        help="Number of chunks for processing the seed index [dependant on block size]")
    diamond_args.add_argument(
        "--max-memory", type=float, default=virtual_memory().available / (1024.0 ** 3),
        help="Maximum memory to use (in Gb) [all available]")
    diamond_args.add_argument(
        "--taxids", default=None, help="Tax IDs to obtain protein sequences of for building a reference database.")
    diamond_args.add_argument(
        '--diamond-mode', help="Mode to run DIAMOND with [fast]", default='fast',
        choices=['fast', 'mid_sensitive', 'sensitive', 'more_sensitive', 'very_sensitive', 'ultra_sensitive'])
    args = parser.parse_args()
    args.output = args.output.rstrip('/')
    args.resources_directory = args.resources_directory.rstrip('/')
    args.columns = args.columns.split('&') if args.columns else None
    if args.taxids is not None:
        args.taxids = args.taxids.split(',')
    return args


def timed_message(message):
    print(f'{strftime("%Y-%m-%d %H:%M:%S", gmtime())}: {message}')


def human_time(seconds):
    days = round(seconds // 86400)
    if days > 0:
        return strftime(f"{days}d%Hh%Mm%Ss", gmtime(seconds))
    return strftime("%Hh%Mm%Ss", gmtime(seconds))


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
    result.columns = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue',
        'bitscore']
    return result


def string4mapping(columns=None):
    if columns is None or columns == []:    # if no columns are inputted, UPIMAPI uses defaults
        valid_columns = [
            'Entry', 'Entry Name', 'Organism', 'Organism (ID)', 'Taxonomic lineage', 'Taxonomic lineage (Ids)',
            'Gene Names', 'Protein names', 'EC number', 'Function [CC]', 'Pathway', 'Keywords',
            'Protein existence', 'Gene Ontology (GO)', 'Protein families', 'BRENDA', 'BioCyc', 'CDD', 'eggNOG',
            'Ensembl', 'InterPro', 'KEGG', 'Pfam', 'Reactome', 'RefSeq', 'UniPathway']
    else:                                   # check what columns are valid
        valid_columns = [column for column in columns if column in columns_dict.keys()]
        invalid_columns = [column for column in columns if column not in columns_dict.keys()]
        for col in invalid_columns:
            print(f'WARNING: "{col}" is not a valid column name. '
                  f'Check https://www.uniprot.org/help/return_fields (Label* column) for valid column names '
                  f'or raise an issue at https://github.com/iquasere/UPIMAPI/issues')
    for col in ['Entry', 'Entry Name']:     # UPIMAPI requires these two columns to be present
        if col not in valid_columns:
            valid_columns.insert(0, col)
    return ','.join([columns_dict[column] for column in valid_columns])


def parallelize(data, func, num_of_processes=8):
    data_split = np.array_split(data, num_of_processes)
    pool = Pool(num_of_processes)
    data = pd.concat(pool.map(func, data_split))
    pool.close()
    pool.join()
    return data


def run_on_subset(func, data_subset, **kwargs):
    return data_subset.apply(func, kwargs)


def parallelize_on_rows(data, func, num_of_processes=8, **kwargs):
    return parallelize(data, partial(run_on_subset, func, kwargs), num_of_processes)


def uniprot_request(ids, columns=None, output_format='tsv'):
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
    fields = f'&fields={string4mapping(columns=columns)}'
    WEBSITE_API = api_info['servers'][0]['url']
    resp = get_url(f"{WEBSITE_API}/uniprotkb/accessions?accessions={','.join(ids)}{fields}&format={output_format}")
    return resp.text


def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))


def submit_id_mapping(from_db, to_db, ids):
    """
    Get info from one database to the other
    :param from_db: Available options at https://rest.uniprot.org/configure/idmapping/fields
    :param to_db: Available options at https://rest.uniprot.org/configure/idmapping/fields
    :param ids:
    :return:
    """
    data = {"from": from_fields[from_db], "to": to_fields[to_db], "ids": ids}
    r = requests.post(f"{api_info['servers'][0]['url']}/idmapping/run", data=data)
    r.raise_for_status()
    return r.json()["jobId"]


def get_id_mapping_results(job_id):
    while True:
        r = get_url(f"{api_info['servers'][0]['url']}/idmapping/status/{job_id}")
        job = r.json()
        if "jobStatus" in job:
            if job["jobStatus"] == "RUNNING":
                sleep(3)
        else:
            return r


def basic_idmapping(ids, from_db, to_db):
    """
    Get info from one database to the other
    :param ids:
    :param from_db: Available options at https://rest.uniprot.org/configure/idmapping/fields
    :param to_db: Available options at https://rest.uniprot.org/configure/idmapping/fields
    :return:
    """
    job_id = submit_id_mapping(from_db, to_db, ids)
    r = get_id_mapping_results(job_id)
    result = pd.DataFrame().from_dict(r.json()["results"])
    while r.links.get("next", {}).get("url"):
        r = get_url(r.links["next"]["url"])
        result = pd.concat([result, pd.DataFrame().from_dict(r.json()["results"])])
    return result


def basic_idmapping_batch(ids, from_db, to_db, step=1000):
    """
    Allows to retrieve millions of IDs at once, there seems to be some limit causing UniProt's API to fail with
    "Request Entity Too Large for url".
    :param step:
    :param ids:
    :return:
    """
    result = pd.DataFrame()
    for i in tqdm(range(0, len(ids), step), desc='Getting valid UniProt IDs', ascii=' >='):
        done = False
        while not done:
            j = min(i + step, len(ids))
            try:
                result = pd.concat([result, basic_idmapping(ids[i:j], from_db, to_db)])
                done = True
            except:
                sleep(3)
    return result


def basic_idmapping_multiprocess(ids, output, from_db, to_db, step=1000, threads=15):
    result = pd.DataFrame()
    ids_groups = split(ids, threads)
    with Manager() as m:
        with m.Pool() as p:
            mapping_results = p.starmap(basic_idmapping_batch, [(
                ids_group, from_db, to_db, step) for ids_group in ids_groups])
    for res in mapping_results:
        result = pd.concat([result, res])
    timed_message(f'{result["from"].unique().size} IDs were successfully mapped.')
    result.to_csv(output, sep='\t', index=False)
    timed_message(f'Results saved at {output}.')


def get_valid_entries(ids):
    job_id = submit_id_mapping("UniProtKB_AC-ID", "UniProtKB", ids)
    r = get_id_mapping_results(job_id)
    valid_entries = [res["from"] for res in r.json()["results"] if '_' not in res["from"]]
    while r.links.get("next", {}).get("url"):
        r = get_url(r.links["next"]["url"])
        valid_entries += [res["from"] for res in r.json()["results"] if '_' not in res["from"]]
    return valid_entries


def get_valid_entries_batch(ids, step=1000):
    """
    Allows to retrieve millions of IDs at once, there seems to be some limit causing UniProt's API to fail with
    "Request Entity Too Large for url".
    :param step:
    :param ids:
    :return:
    """
    valid_entries = []
    for i in tqdm(range(0, len(ids), step), desc='Getting valid UniProt IDs', ascii=' >='):
        done = False
        while not done:
            j = min(i + step, len(ids))
            try:
                valid_entries += get_valid_entries(ids[i:j])
                done = True
            except:
                sleep(3)
    return valid_entries


def get_valid_entries_multiprocess(ids, step=1000, threads=15):
    valid_entries = []
    ids_groups = split(ids, threads)
    with Manager() as m:
        with m.Pool() as p:
            result = p.starmap(get_valid_entries_batch, [(ids_group, step) for ids_group in ids_groups])
    for res in result:
        valid_entries += res
    not_valid = [ide for ide in ids if ide not in valid_entries]
    # take, from the valid IDs, the part after the dot, as this invalidates them
    valid_entries = [entry.split('.')[0] for entry in valid_entries]
    timed_message(f'{len(valid_entries)} UniProt IDs identified as valid.')
    return valid_entries, not_valid


def get_uniprot_information(ids, step=1000, sleep_time=30, columns=None, max_tries=3):
    """
    Input:
        ids: list of UniProt IDs to query
        step: INT, number of IDs to send per request
        sleep_time: INT, number of seconds to wait between requests
        columns: list - names of UniProt columns to get info on
    Output:
        pd.DataFrame will be returned with the information about the IDs queried.
    """
    result = pd.DataFrame()
    for i in tqdm(range(0, len(ids), step), desc=f'Retrieving UniProt information from {len(ids)} IDs'):
        tries = 0
        done = False
        j = min(i + step, len(ids))
        while not done and tries < max_tries:
            try:
                data = uniprot_request(ids[i:j], columns=columns)
                if len(data) > 0:
                    uniprotinfo = pd.read_csv(StringIO(data), sep='\t')
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
        step: INT, number of IDs to send per request
        sleep_time: INT, number of seconds to wait between requests
    Output:
        str object containing the fasta sequences and headers
        of the proteis belonging to the IDs queried will be returned
    """
    result = ''
    for i in tqdm(range(0, len(ids), step), desc=f"Building FASTA from {len(ids)} IDs."):
        j = min(i + step, len(ids))
        data = uniprot_request(ids[i:j], output_format='fasta')
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
        ids_done = []
    ids_missing = list(set(all_ids) - set(ids_done))

    tries = 0
    ids_done = ([ide.split('|')[1] for ide in get_fasta_ids(output)] if os.path.isfile(output) else [])
    while len(ids_done) < len(all_ids) and tries < max_iter:
        ids_missing = list(set([ide for ide in tqdm(all_ids, desc='Checking which IDs are missing information.')
                                if ide not in ids_done]))
        print(f'Information already gathered for {int(len(ids_done) / 2)} ids. Still missing for {len(ids_missing)}.')
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


def check_ids_already_done(output, ids):
    if os.path.isfile(output) and os.stat(output).st_size > 1:
        try:
            result = pd.read_csv(output, sep='\t', low_memory=False).drop_duplicates()
            print(f'{output} was found. Will perform mapping for the remaining IDs.')
            ids_done = list(set(result['Entry'].tolist() + result['Entry Name'].tolist()))
        except OSError:
            print(f'{output} was found. However, it could not be parsed. Will restart mapping.')
            result = pd.DataFrame()
            ids_done = []
    else:
        print(f'{output} not found or empty. Will perform mapping for all IDs.')
        result = pd.DataFrame()
        ids_done = []
    ids_missing = list(set(ids) - set(ids_done))
    print(f'IDs present in uniprotinfo file: {int(len(ids_done) / 2)}')      # entry and entry name count by 2
    print(f'IDs missing: {len(ids_missing)}')
    return ids_done, ids_missing, result


def extract_lineage_columns(columns):
    tax_cols, taxids_cols = [], []
    if columns is not None:
        tax_cols = [col for col in columns if ('Taxonomic lineage (' in col and col != 'Taxonomic lineage (Ids)')]
        taxids_cols = [col for col in columns if ('Taxonomic lineage IDs (' in col)]
        columns = [col for col in columns if col not in tax_cols + taxids_cols]
        if len(tax_cols) > 0:
            columns.append('Taxonomic lineage')
        if len(taxids_cols) > 0:
            columns.append('Taxonomic lineage (Ids)')
    return columns, tax_cols, taxids_cols, tax_cols + taxids_cols


def make_taxonomic_lineage_df(tax_lineage_col, prefix='Taxonomic lineage IDs'):
    """
    Parses the taxonomic lineage column of the uniprotinfo dataframe and returns a dataframe with the taxonomic lineage
    separated in columns.
    :param tax_lineage_col: pd.Series with the taxonomic lineage column of the uniprotinfo dataframe
    :param prefix: str, prefix to use for the columns of the new dataframe
    :return: pd.DataFrame with the taxonomic lineage separated in columns
    """
    # First, split records by ', '
    split_regex = r"(?<=\)),"
    result = pd.DataFrame.from_records(tax_lineage_col.apply(lambda x: re.split(split_regex, x)).apply(
        # Then, split each record by ' ('
        lambda x: [part[:-1].split(' (') for part in x]).apply(
        # Finally, build dictionary with the taxonomic level as key and the taxonomy as value. ' ('.join avoids cases
        # where the taxonomy has a '(' in it (e.g. 'Clostridium scindens (strain JCM 10418 / VPI 12708) (species)')
        lambda x: {part[-1]: ' ('.join(part[:-1]) for part in x if part[-1] != 'no rank'}))
    # Rename columns in old UniProt fashion
    result = result.rename(columns={col: f'{prefix} ({col.upper()})' for col in result.columns})
    return result


def uniprot_information_workflow(ids, output, max_iter=5, columns=None, step=1000, sleep_time=10):
    ids_done, ids_missing, result = check_ids_already_done(output, ids)
    add_tax_cols = columns is None
    tries, last_ids_missing = 0, None
    ids_unmapped_output = f"{'/'.join(output.split('/')[:-1])}/ids_unmapped.txt"
    columns, tax_cols, taxids_cols, all_tax_cols = extract_lineage_columns(columns)
    default_tax_cols = [
        'Taxonomic lineage (SUPERKINGDOM)', 'Taxonomic lineage (PHYLUM)', 'Taxonomic lineage (CLASS)',
        'Taxonomic lineage (ORDER)', 'Taxonomic lineage (FAMILY)', 'Taxonomic lineage (GENUS)']
    if add_tax_cols:
        tax_cols += default_tax_cols
        all_tax_cols += default_tax_cols
    uniprotinfo = pd.DataFrame()
    while len(ids_missing) > 0 and tries < max_iter and ids_missing != last_ids_missing:
        print(f'Information already gathered for {int(len(ids_done) / 2)} ids. Still missing for {len(ids_missing)}.')
        last_ids_missing = ids_missing
        info = get_uniprot_information(
            ids_missing, step=step, columns=columns, max_tries=max_iter, sleep_time=sleep_time)
        info.reset_index(inplace=True)
        del info['index']
        if len(info) > 0:
            ids_done += list(set(info['Entry'].tolist() + info['Entry Name'].tolist()))
            uniprotinfo = pd.concat([uniprotinfo, info], ignore_index=True)
        ids_missing = list(set(last_ids_missing) - set(ids_done))
        if len(ids_missing) > 0:
            if last_ids_missing == ids_missing:
                print("Could not map additional IDs for this mapping. There were probably some outdated IDs. "
                      "For more questions, please contact through https://github.com/iquasere/UPIMAPI/issues")
            else:
                print('Failed to retrieve information for some IDs. Retrying request.')
                tries += 1
    if len(uniprotinfo) > 0:
        tax_df = pd.DataFrame()
        if len(tax_cols) > 0:
            tax_df = make_taxonomic_lineage_df(uniprotinfo['Taxonomic lineage'], prefix='Taxonomic lineage')
        if len(taxids_cols) > 0:
            tax_df = pd.concat([tax_df, make_taxonomic_lineage_df(
                uniprotinfo['Taxonomic lineage (Ids)'], prefix='Taxonomic lineage IDs')], axis=1)
        # remove taxonomic lineage and taxonomic lineage (ids) columns if they were not supposed to be added
        for col in all_tax_cols:
            if col not in tax_df.columns:
                del tax_df[col]
        uniprotinfo = pd.concat([uniprotinfo, tax_df[all_tax_cols]], axis=1).rename(columns={
            'Organism': 'Taxonomic lineage (SPECIES)'})
        cols = uniprotinfo.columns.tolist() + [col for col in result.columns.tolist() if col not in uniprotinfo.columns]
        if 'Taxonomic lineage (SPECIES)' in cols:
            cols.remove('Taxonomic lineage (SPECIES)')
            cols.append('Taxonomic lineage (SPECIES)')
        if 'index' in result.columns:
            del result['index']
    else:
        cols = result.columns.tolist()
    result = pd.concat([result, uniprotinfo], ignore_index=True)
    result[cols].to_csv(output, sep='\t', index=False)
    if len(ids_missing) == 0:
        print(f'Results for all IDs are available at {output}')
    else:
        open(ids_unmapped_output, 'a').write('\n'.join(ids_missing))
        print(f"Maximum iterations were made. Results related to {str(len(ids_missing))} IDs were not obtained. "
              f"IDs with missing information are available at {ids_unmapped_output} and information obtained is "
              f"available at {output}")
    return result


def determine_full_id(ids):
    for ide in ids:
        if '|' in ide:
            return True
    return False


def parse_fasta(file):
    with open(file) as f:
        lines = [line.rstrip('\n') for line in f]
    i = 0
    sequences = {}
    while i < len(lines):
        if lines[i].startswith('>'):
            name = lines[i][1:]
            sequences[name] = ''
            i += 1
            while i < len(lines) and not lines[i].startswith('>'):
                sequences[name] += lines[i]
                i += 1
    return sequences


def get_ids(args_input, input_type, full_id='auto'):
    if input_type == 'blast':
        ids = parse_blast(args_input)['sseqid'].tolist()
    elif input_type == 'txt':
        ids = []
        with open(args_input) as f:
            preids = f.read().split('\n')
        for preid in preids:
            ids += preid.split(',')
    elif input_type == 'fasta':
        ids = parse_fasta(args_input).keys()
    else:       # if PIPE
        ids = args_input.split(',')
    if full_id == 'auto':
        full_id = determine_full_id(ids)
        print(f'Auto determined "full id" as: {full_id}')
    if full_id:
        return_ids = [ide.split('|')[1] for ide in ids if ide not in ['*', '']]
        sp_ids = [ide.split('|')[1] for ide in ids if ide.startswith('sp')]
        return return_ids, full_id, sp_ids
    return_ids = [ide for ide in ids if ide not in ['*', '']]
    return return_ids, full_id, return_ids
    # second return_ids is just mock to return the same type of output as for sp_ids


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


def block_size_and_index_chunks(argsb, argsc, memory):
    if argsb:
        b = argsb
    else:
        b = memory / 20  # b = memory in Gb / 20
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
                bit_score=None, pident=None, mode='fast'):
    command = (
        f"diamond blastp --query {query} --out {aligned} --un {unaligned} --db {database} --outfmt 6 --unal 1 "
        f"--threads {threads} --max-target-seqs {max_target_seqs} -b {b} -c {c} --evalue {e_value} "
        f"--{mode.replace('_', '-')}")
    if bit_score:
        command += f' --min-score {bit_score}'
    if pident:
        command += f' --id {pident}'
    run_command(command)


def get_proteome_for_taxid_slow(taxid, max_tries=3):
    """
    Get proteome for taxid the "proper" way. It is very slow, though, so not used.
    :param taxid:
    :param max_tries:
    :return:
    """
    tries = 0
    res = requests.get(f'https://rest.uniprot.org/uniprotkb/search?format=fasta&query=%28taxonomy_id%3A{taxid}%29')
    result = res.content.decode('utf8')
    pages = 0
    while tries < max_tries and res.links != {}:
        try:
            res = requests.get(res.links['next']['url'])
            result += res.content.decode('utf8')
            pages += 1
        except:
            tries += 1
            sleep(10)
    return result


def get_proteome_for_taxid(taxid, max_tries=3):
    tries = 0
    url = f'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=taxonomy_id:{taxid}'
    while tries < max_tries:
        try:
            return requests.get(url).content.decode('utf8')
        except:
            print(f'Failed! {max_tries - tries} tries remaining.')
            tries += 1
            sleep(10)


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
    block_size = 102400  # 100 Kibibytes
    progress_bar = tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True, desc=f'Downloading {url.split("/")[-1]}')
    with open(f'{output_folder}/{url.split("/")[-1]}', 'wb') as file:
        for data in response.iter_content(block_size):
            progress_bar.update(len(data))
            file.write(data)
    progress_bar.close()


def download_uniprot(output_folder, mirror='expasy'):
    base_urls = {
        'expasy': 'https://ftp.expasy.org',
        'uniprot': 'https://ftp.uniprot.org/pub',
        'ebi': 'https://ftp.ebi.ac.uk/pub'
    }
    for file in ["uniprot_sprot.fasta.gz", "uniprot_trembl.fasta.gz", "reldate.txt"]:
        print(f'Downloading and writing: {file}')
        download_with_progress_bar(
            f'{base_urls[mirror]}/databases/uniprot/current_release/knowledgebase/complete/{file}', output_folder)
    run_pipe_command(f'zcat {output_folder}/uniprot_trembl.fasta.gz {output_folder}/uniprot_sprot.fasta.gz > '
                     f'{output_folder}/uniprot.fasta')
    for file in [f'{output_folder}/uniprot_trembl.fasta.gz', f'{output_folder}/uniprot_sprot.fasta.gz']:
        os.remove(file)


def build_reference_database(database, output_folder, taxids=None, max_tries=3, mirror='expasy'):
    if database == 'uniprot':
        download_uniprot(output_folder, mirror=mirror)
    elif database == 'swissprot':
        download_with_progress_bar(
            "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz",
            output_folder)
        run_command(f'gunzip {output_folder}/uniprot_sprot.fasta.gz')
    elif database == 'taxids':
        for taxid in tqdm(taxids, desc=f'Retrieving reference proteomes for {len(taxids)} taxa from UniProt.'):
            with open(f'{output_folder}/taxids_database.fasta', 'a') as f:
                f.write(get_proteome_for_taxid(taxid, max_tries=max_tries))


def must_build_database(database, resources_folder):
    db2suffix = {'uniprot': 'uniprot.fasta', 'swissprot': 'uniprot_sprot.fasta', 'taxids': 'taxids_database.fasta'}
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
    result, ids_found = [], []
    i = 1
    record = next(sp_dat)
    number_of_entries = count_on_file('Reviewed;', sp_dat_filename)
    while record is not None and len(ids) > 0:
        match_id = get_match_id(record, ids)
        if match_id is not None:
            result.append(record.__dict__)
            ids_found.append(match_id)
        if i % 100000 == 0:
            print(f'[{i}/{number_of_entries}] SwissProt entries queried')
        record = next(sp_dat, None)
        i += 1
    print(f'[{i}/{number_of_entries}] SwissProt entries queried')
    return pd.DataFrame(result), ids_found


def lineage_to_columns(lineage, tax_tsv):
    l2c_result, l2c_taxids = {}, {}
    for taxon in lineage:
        match = tax_tsv.loc[taxon, ["rank", "taxid"]]
        if type(match) == pd.core.series.Series:
            rank, taxid = match[["rank", "taxid"]]
            if type(rank) == str:
                l2c_result[f'Taxonomic lineage ({rank.upper()})'] = taxon
                l2c_taxids[f'Taxonomic identifier ({rank.upper()})'] = taxid
        else:   # some taxIDs have multiple levels (e.g. "Craniata")
            for i in range(len(match)):
                rank, taxid = match.iloc[i][["rank", "taxid"]]
                if type(rank) == str:
                    l2c_result[f'Taxonomic lineage ({rank.upper()})'] = taxon
                    l2c_taxids[f'Taxonomic identifier ({rank.upper()})'] = taxid
    l2c_result['Taxonomic lineage (ALL)'] = ', '.join(set(l2c_result.values()))
    l2c_taxids['Taxonomic identifier (ALL)'] = ', '.join(set(l2c_taxids.values()))
    l2c_result = {**l2c_result, **l2c_taxids, 'index': lineage}
    return l2c_result


def lineages_to_columns(lineages, tax_tsv):
    """
    Does the same as lineage_to_columns, but to all lineages, instead of a single one
    :param lineages:
    :param tax_tsv:
    :return:
    """
    return [lineage_to_columns(lineage, tax_tsv) for lineage in lineages]


def get_upper_taxids(taxid, tax_df):
    """
    :param taxid: str - taxID to get upper taxIDs from
    :param tax_df: pd.DataFrame - of read taxonomy.tsv (from taxonomy.rdf)
    :returns list - of upper taxIDs
    """
    if taxid == '0':
        return []
    taxids = []
    while taxid != '1' and taxid != 'Taxon':
        taxids.append(taxid)
        taxid = tax_df.loc[taxid]['parent_taxid']
    return taxids


def parse_taxonomy(data, tax_tsv_df, threads=15):
    tax_tsv_df.set_index('name', inplace=True)
    tax_tsv_df['taxid'] = tax_tsv_df['taxid'].astype(str)
    all_classifications = split(data['organism_classification'].drop_duplicates().tolist(), threads)
    with Manager() as m:
        with m.Pool() as p:
            result = p.starmap(lineages_to_columns, [(classifications, tax_tsv_df)
                                                     for classifications in all_classifications])
    decompacted = []
    for res in result:
        decompacted += res
    return pd.DataFrame(decompacted).set_index('index')


def parse_comments(sp_data):
    result = []
    bpc_list = []
    for comments in sp_data['comments']:
        partial = {key: '' for key in [
            'FUNCTION', 'SUBUNIT', 'INTERACTION', 'SUBCELLULAR LOCATION', 'ALTERNATIVE PRODUCTS', 'TISSUE SPECIFICITY',
            'PTM', 'POLYMORPHISM', 'DISEASE', 'MISCELLANEOUS', 'SIMILARITY', 'CAUTION', 'SEQUENCE CAUTION',
            'WEB RESOURCE', 'MASS SPECTROMETRY', 'RNA EDITING', 'CATALYTIC ACTIVITY', 'COFACTOR', 'ACTIVITY REGULATION',
            'PATHWAY', 'DEVELOPMENTAL STAGE', 'INDUCTION', 'ALLERGEN', 'BIOTECHNOLOGY', 'DISRUPTION PHENOTYPE',
            'PHARMACEUTICAL', 'TOXIC DOSE', 'DOMAIN']}
        bpc_dict = {'Kinetic parameters': []}
        for comment in comments:
            comment = comment.split(': ')
            if comment[0] in partial.keys():
                partial[comment[0]] += f'{": ".join(comment)} '
            else:
                if comment[0] in ['BIOPHYSICOCHEMICAL PROPERTIES']:
                    if comment[1] not in bpc_dict.keys():
                        bpc_dict[comment[1]] = [f'{comment[0]}: {comment[1]}: {comment[2]}']
                    else:
                        bpc_dict[comment[1]].append(f'{comment[0]}: {comment[1]}: {comment[2]}')
                else:
                    print(f'Comment still not implemented: [{comment[0]}]')
        result.append(partial)
        bpc_dict['Kinetics'] = bpc_dict.pop('Kinetic parameters')
        bpc_list.append(bpc_dict)
    result = pd.DataFrame(result, columns=[
        'Function [CC]', 'Subunit structure [CC]', 'Interacts with', 'Subcellular location [CC]',
        'Alternative products (isoforms)', 'Tissue specificity', 'Post-translational modification', 'Polymorphism',
        'Involvement in disease', 'Miscellaneous [CC]', 'Sequence similarities', 'Caution', 'Sequence caution',
        'Web resources', 'Mass spectrometry', 'RNA editing', 'Catalytic activity', 'Cofactor', 'Activity regulation',
        'Pathway', 'Developmental stage', 'Induction', 'Allergenic properties', 'Biotechnological use',
        'Disruption phenotype', 'Pharmaceutical use', 'Toxic dose', 'Domain [CC]'])
    result['Erroneous gene model prediction'] = result['Sequence caution']
    bpc_df = pd.DataFrame(bpc_list)
    for col in bpc_df:
        bpc_df[col] = bpc_df[col].apply(lambda x: '; '.join(x) if type(x) == list else x)
    return pd.concat([result, bpc_df], axis=1)


def add_to_dict(dictionary, key, value):
    if key in dictionary.keys():
        dictionary[key] += value
    else:
        dictionary[key] = value


def cross_references_to_columns(cross_refs):
    result = {}
    go_dict = {}
    go_rel = {'C': 'cellular component', 'F': 'molecular function', 'P': 'biological process'}
    for ref in cross_refs:
        if ref[0] == 'GO':
            refie = ref[2].split(':')
            add_to_dict(go_dict, f'Gene ontology ({go_rel[refie[0]]})', f'{refie[1]} [{ref[1]}]; ')
            add_to_dict(go_dict, 'Gene ontology (GO)', f'{refie[1]} [{ref[1]}]; ')
            add_to_dict(go_dict, 'Gene ontology IDs', f'{ref[1]}; ')
        else:
            if ref[0] == 'Proteomes':
                value = f'{ref[1]}: {ref[2]}'
            else:
                value = f'{ref[1]};'
            add_to_dict(result, ref[0], value)
    return result, go_dict


def parse_cross_references(sp_data):
    ref_result = [cross_references_to_columns(cross_refs) for cross_refs in sp_data['cross_references']]
    ref_dict, go_dict = zip(*ref_result)
    ref_df = pd.DataFrame(ref_dict)
    ref_df.columns = map(lambda x: f'Cross-reference ({x})' if x != 'Proteomes' else x, ref_df.columns)
    go_df = pd.DataFrame(go_dict)
    ref_df = pd.concat([ref_df, go_df], axis=1)
    return ref_df


def gene_name_to_columns(genes):
    if genes == '':
        info = {}
    else:
        info = [pair.split('=') for pair in genes.rstrip(';').split('; ')]
        info = {pair[0]: pair[1] for pair in info}
    return {'Gene names': ' '.join(info.values()) if info != {} else '',
            'Gene names  (ordered locus )': info['OrderedLocusNames'] if 'OrderedLocusNames' in info else '',
            'Gene names  (ORF )': info['ORFNames'] if 'ORFNames' in info else '',
            'Gene names  (primary )': info['Name'] if 'Name' in info else '',
            'Gene names  (synonym )': info['Synonyms'] if 'Synonyms' in info else ''}


def parse_gene_names(sp_data):
    return pd.DataFrame([gene_name_to_columns(genes) for genes in sp_data['gene_name']])


def parse_description_text(description):
    result = {}
    parts = description[:-1].split('; ')
    i = 0
    while i < len(parts):
        parted = parts[i].split('=')
        if parted[0].startswith('RecName: Full'):
            result['RecName'] = {}
            result['RecName']['Full'] = parted[1]
            i += 1
            while i < len(parts) and ':' not in parts[i].split()[0]:
                parted = parts[i].split('=')
                result['RecName'][parted[0]] = parted[1]
                i += 1
        elif parted[0].startswith('AltName: Full'):
            if 'AltName' not in result.keys():
                result['AltName'] = []
            altname = {'Full': parted[1]}
            i += 1
            while i < len(parts) and ':' not in parts[i].split()[0]:
                parted = parts[i].split('=')
                altname[parted[0]] = parted[1]
                i += 1
            result['AltName'].append(altname)
        elif parted[0].startswith('Contains: RecName'):
            if 'Contains' not in result.keys():
                result['Contains'] = []
            contains = {'RecName': {'Full': parted[1]}}
            i += 1
            while i < len(parts) and ':' not in parts[i].split()[0]:
                parted = parts[i].split('=')
                contains['RecName'][parted[0]] = parted[1]
                i += 1
            result['Contains'].append(contains)
        elif parts[i].startswith('Flags'):
            parted = parts[i].split(': ')
            if 'Flags' in result.keys():
                result['Flags'] = parted[1]
            else:
                result['Flags'] = [parted[1]]
            i += 1
        else:
            #print('A description UPIMAPI cannot yet handle!')
            #print(parts[i])
            i += 1
    return result


def fix_term(term):
    return term if '{ECO:' not in term else ' '.join(term.split()[:-1])


def parse_descriptions(sp_data):
    desc_data_df = sp_data['description'].apply(parse_description_text)
    description_df = pd.DataFrame()
    description_df['Protein names'] = desc_data_df.apply(
        lambda x: '{}{}{}{}{}{}'.format(
            fix_term(x['RecName']['Full']), f" ({fix_term(x['RecName']['Short'])})" if 'Short' in x['RecName'].keys()
            else "", f" (EC {fix_term(x['RecName']['EC'])})" if 'EC' in x['RecName'].keys() else "",
            ' ' + ' '.join(' '.join([f"({fix_term(value)})" for value in altname.values()]) for altname in x['AltName'])
            if 'AltName' in x.keys() else "",
            f" [Cleaved into: {'; '.join([fix_term(v['RecName']['Full']) for v in x['Contains']])}]"
            if 'Contains' in x.keys() else "",
            ' '.join([f" ({flag})" for flag in x['Flags']]) if 'Flags' in x.keys() else ""))
    description_df['EC number'] = desc_data_df.apply(
        lambda x: x['RecName']['EC'].split()[0] if type(x) != float and 'RecName' in x.keys()
        and 'EC' in x['RecName'].keys() else np.nan)
    return description_df


def parse_feature(feature, position, qualifiers=True, ide=True):
    """
    :param feature: str - the feature itself
    :param position: str - position information
    :param qualifiers: bool - add qualifiers information?
    :param ide: bool - add id information?
    :return: str - the term to add
    """
    result = f'{feature.type} {position}'
    if qualifiers:
        result += '  ' + "  ".join([f'/{key}="{value}";' for key, value in feature.qualifiers.items()])
    if ide:
        result += '  ' + f'  /id="{feature.id}";'
    return result


def parse_features(sp_data):
    feats_list = []
    pos_funcs = {
        'all': lambda x: f'{"?" if x.location.start.position is None else x.location.start.position + 1}..'
                         f'{x.location.end.position};  ',
        'end': lambda x: f'{feature.location.end.position};  '}
    prefix2info = {
        'VAR_SEQ': ('Alternative sequence', 'all', True, True),
        'VARIANT': ('Natural variant', 'end', True, True),
        'NON_CONS': ('Non-adjacent residues', 'all', True, False),
        'NON_STD': ('Non-standard residue', 'end', True, False),
        'NON_TER': ('Non-terminal residue', 'end', False, False),
        'CONFLICT': ('Sequence conflict', 'all', True, False),
        'UNSURE': ('Sequence uncertainty', 'end', True, False),
        'ACT_SITE': ('Active site', 'end', True, False),
        'BINDING': ('Binding site', 'end', True, False),
        'DNA_BIND': ('DNA binding', 'all', True, False),
        'METAL': ('Metal binding', 'end', True, False),
        'NP_BIND': ('Nucleotide binding', 'all', True, False),
        'SITE': ('Site', 'end', True, False),
        'INTRAMEM': ('Intramembrane', 'all', True, False),
        'TOPO_DOM': ('Topological domain', 'all', True, False),
        'TRANSMEM': ('Transmembrane', 'all', True, False),
        'CHAIN': ('Chain', 'all', True, True),
        'CROSSLNK': ('Cross-link', 'all', True, False),
        'DISULFID': ('Disulfide bond', 'all', True, False),
        'CARBOHYD': ('Glycosylation', 'end', True, False),
        'INIT_MET': ('Initiator methionine', 'end', True, False),
        'LIPID': ('Lipidation', 'end', True, False),
        'MOD_RES': ('Modified residue', 'end', True, False),
        'PEPTIDE': ('Peptide', 'all', True, True),
        'PROPEP': ('Propeptide', 'all', False, True),
        'SIGNAL': ('Signal peptide', 'all', True, False),
        'TRANSIT': ('Transit peptide', 'all', True, False),
        'STRAND': ('Beta strand', 'all', True, False),
        'HELIX': ('Helix', 'all', True, False),
        'TURN': ('Turn', 'all', True, False),
        'COILED': ('Coiled coil', 'all', True, False),
        'COMPBIAS': ('Compositional bias', 'all', True, False),
        'DOMAIN': ('Domain [FT]', 'all', True, False),
        'MOTIF': ('Motif', 'all', True, False),
        'REGION': ('Region', 'all', True, False),
        'REPEAT': ('Repeat', 'all', True, False),
        'ZN_FING': ('Zinc finger', 'all', True, False),
        'MUTAGEN': ('Mutagenesis', 'end', True, False),
        'CA_BIND': ('Calcium binding', 'all', True, False)}
    count_features = {}
    for features in sp_data['features']:
        feats_dict = {}
        for feature in features:
            if feature.type in prefix2info.keys():
                parameters = prefix2info[feature.type]
                if parameters[0] not in feats_dict.keys():
                    feats_dict[parameters[0]] = parse_feature(
                    feature, pos_funcs[parameters[1]](feature), qualifiers=parameters[2], ide=parameters[3])
                    count_features[parameters[0]] = 1
                else:
                    feats_dict[parameters[0]] += ' ' + parse_feature(
                        feature, pos_funcs[parameters[1]](feature), qualifiers=parameters[2], ide=parameters[3])
                    count_features[parameters[0]] += 1
            else:
                print(f'A feature UPIMAPI can yet not handle! [{feature.type}]')
        feats_dict['Features'] = '; '.join([f'{feat_type} ({count})' for feat_type, count in count_features.items()])
        feats_list.append(feats_dict)
    return pd.DataFrame(feats_list)


def parse_host_taxonomy_id(sp_data, tax_tsv):
    tax_tsv = tax_tsv.reset_index().set_index('taxid')
    return sp_data['host_taxonomy_id'].apply(
        lambda x: '; '.join([f'{tax_tsv.loc[tid, "name"]} [TaxID: {tid}]' for tid in x]) if len(x) > 0 else np.nan)


def parse_sp_data(sp_data, tax_tsv, threads=15):
    """
    Parses data from local ID mapping through DAT file
    :param threads:
    :param sp_data: pandas.DataFrame
    :param tax_tsv: str - filename of taxonomy in TSV format
    :return: pandas.DataFrame - organized in same columns as data from UniProt's API
    """
    if len(sp_data) == 0:
        return pd.DataFrame()
    tax_tsv_df = pd.read_csv(tax_tsv, sep='\t', dtype={'taxid': str, 'name': str, 'rank': str, 'parent_taxid': str})
    tax_tsv_df = tax_tsv_df[tax_tsv_df.name.notnull()]
    result = pd.DataFrame()
    result['Entry'] = sp_data['accessions'].apply(lambda x: x[0])
    local2api = {
        'entry_name': 'Entry name',
        'data_class': 'Status',
        'sequence_length': 'Length',
        'sequence': 'Sequence'
    }
    for k, v in local2api.items():
        if v not in [None, False]:
            result[v] = sp_data[k]
    result['Organism ID'] = result['Taxonomic identifier (SPECIES)'] = \
        sp_data['taxonomy_id'].apply(lambda x: x[0] if len(x) > 0 else x)
    result['Virus hosts'] = sp_data['host_organism'].apply(lambda x: x[0] if len(x) > 0 else x)
    result['Keywords'] = sp_data['keywords'].apply(';'.join)
    result['Organism'] = sp_data['organism'].str.rstrip('.')
    result['Taxonomic lineage (SPECIES)'] = result['Organism'].apply(lambda x: ' '.join(x.split()[:2]))
    timed_message('Parsing taxonomy (this may take a while)')
    tax_df = parse_taxonomy(sp_data, tax_tsv_df, threads=threads).reset_index()
    rel_df = sp_data['organism_classification'].apply(','.join)
    tax_df['index'] = tax_df['index'].apply(','.join)
    rel_df = pd.merge(rel_df, tax_df, left_on='organism_classification', right_on='index', how='left')
    del rel_df['organism_classification']
    del rel_df['index']
    result['Virus hosts'] = parse_host_taxonomy_id(sp_data, tax_tsv_df)
    result = pd.concat([result, rel_df], axis=1)
    timed_message('Parsing genes')
    result = pd.concat([result, parse_gene_names(sp_data)], axis=1)
    timed_message('Parsing cross-references')
    result = pd.concat([result, parse_cross_references(sp_data)], axis=1)
    timed_message('Parsing comments')
    result = pd.concat([result, parse_comments(sp_data)], axis=1)
    timed_message('Parsing features')
    result = pd.concat([result, parse_features(sp_data)], axis=1)
    result = pd.concat([result, parse_descriptions(sp_data)], axis=1)
    result['Gene encoded by'] = sp_data['organelle'].str.rstrip('.')
    result['Mass'] = sp_data['seqinfo'].apply(lambda x: x[1])
    result['Date of creation'] = sp_data['created'].apply(
        lambda x: datetime.strptime(x[0], '%d-%b-%Y').strftime('%Y-%m-%d'))
    result['Date of last modification'] = sp_data['annotation_update'].apply(
        lambda x: datetime.strptime(x[0], '%d-%b-%Y').strftime('%Y-%m-%d'))
    result['Version (entry)'] = sp_data['annotation_update'].apply(lambda x: x[1])
    result['Date of last sequence modification'] = sp_data['sequence_update'].apply(
        lambda x: datetime.strptime(x[0], '%d-%b-%Y').strftime('%Y-%m-%d'))
    result['Version (sequence)'] = sp_data['sequence_update'].apply(lambda x: x[1])
    result['PubMed ID'] = sp_data['references'].apply(
        lambda x: '; '.join([ref.references[0][1] for ref in x if len(ref.references) > 0]))
    return result


def get_sprot_dat(sp_dat):
    run_command(
        f'wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/'
        f'uniprot_sprot.dat.gz -O {sp_dat}.gz')
    run_command(f'gunzip {sp_dat}.gz')


def local_id_mapping(ids, sp_dat, tax_tsv, output, columns=None, databases=None, threads=15):
    ids_done, ids_missing, result = check_ids_already_done(output, ids)
    if len(ids_missing) == 0:
        return set()
    if not os.path.isfile(sp_dat):
        timed_message(f'Creating {sp_dat}')
        get_sprot_dat(sp_dat)
    if not os.path.isfile(tax_tsv):
        timed_message(f'Creating {tax_tsv}')
        get_tabular_taxonomy(tax_tsv)
    timed_message('Searching for IDs in SwissProt DAT')
    sp_data, ids_found = get_local_swissprot_data(sp_dat, ids_missing)
    timed_message('Parsing SwissProt resuls')
    sp_parsed = parse_sp_data(sp_data, tax_tsv, threads=threads)
    result = pd.concat([result, sp_parsed])
    columns = [col for col in columns if col in result.columns.tolist()]
    databases = [db for db in databases if db in result.columns.tolist()]
    result[columns + databases].to_csv(output, sep='\t', index=False)
    return ids_found


def get_input_type(args_input, blast=True):
    if args_input is None:
        return input('IDs to perform mapping on (comma separated values):'), 'stdin'
    if blast:
        return args_input, 'blast'
    if check_output(f"head -c 1 {args_input}", shell=True).decode('utf8') == '>':
        return args_input, 'fasta'
    return args_input, 'txt'


def check_no_annotation(args_input, no_annotation):
    if args_input is None:
        is_fasta = False
    else:
        is_fasta = check_output(f"head -c 1 {args_input}", shell=True).decode('utf8') == '>'
    if not is_fasta:
        no_annotation = True
    return no_annotation


def blast_consensus(alignment_file):
    blast = parse_blast(alignment_file)
    query_to_ref, ref_to_query, res = {}, {}, {}
    with open(alignment_file) as file:
        line = file.readline()
        while line:
            line = line.strip('\n').split('\t')
            query_seq, ref_seq, evalue = line[0], line[1], float(line[-2])
            if query_seq not in query_to_ref:
                if ref_seq not in ref_to_query:
                    query_to_ref[query_seq] = {'ref_seq': ref_seq, 'evalue': evalue}
                    ref_to_query[ref_seq] = {'query_seq': query_seq, 'evalue': evalue}
                else:
                    if ref_to_query[ref_seq]['evalue'] > evalue:
                        ref_to_query[ref_seq] = {'query_seq': query_seq, 'evalue': evalue}
                        query_to_ref[query_seq] = {'ref_seq': ref_seq, 'evalue': evalue}
            else:
                if ref_seq not in ref_to_query:
                    if query_to_ref[query_seq]['evalue'] > evalue:
                        ref_to_query[ref_seq] = {'query_seq': query_seq, 'evalue': evalue}
                        query_to_ref[query_seq] = {'ref_seq': ref_seq, 'evalue': evalue}
                else:
                    if ref_to_query[ref_seq]['evalue'] > evalue:
                        if query_to_ref[query_seq]['evalue'] > evalue:
                            ref_to_query[ref_seq] = {'query_seq': query_seq, 'evalue': evalue}
                            query_to_ref[query_seq] = {'ref_seq': ref_seq, 'evalue': evalue}
            line = file.readline()
    for query_seq in query_to_ref:
        ref_seq = query_to_ref[query_seq]['ref_seq']
        if query_seq == ref_to_query[ref_seq]['query_seq']:
            res[query_seq] = query_to_ref[query_seq]['ref_seq']
    res = pd.DataFrame.from_dict(res, orient='index').reset_index()
    res.columns = ['qseqid', 'sseqid']
    return blast.set_index(['qseqid', 'sseqid']).loc[res.set_index(['qseqid', 'sseqid']).index].reset_index()


def upimapi():
    args = get_arguments()
    Path(args.output).mkdir(parents=True, exist_ok=True)
    Path(args.resources_directory).mkdir(parents=True, exist_ok=True)
    args.no_annotation = check_no_annotation(args.input, args.no_annotation)

    # Annotation with DIAMOND
    if not args.no_annotation:
        db2file = {'uniprot': f'{args.resources_directory}/uniprot.fasta',
                   'swissprot': f'{args.resources_directory}/uniprot_sprot.fasta',
                   'taxids': f'{args.resources_directory}/taxids_database.fasta'}
        if args.database in db2file.keys():
            database = db2file[args.database]
        else:
            database = args.database

        if not args.skip_db_check:
            if must_build_database(args.database, args.resources_directory):
                build_reference_database(
                    args.database, args.resources_directory, taxids=args.taxids, max_tries=args.max_tries,
                    mirror=args.mirror)
        if not database.endswith(".dmnd"):
            diamond_formatted = f"{'.'.join(database.split('.')[:-1])}.dmnd"
            if not os.path.isfile(diamond_formatted):
                make_diamond_database(database, diamond_formatted)
            database = diamond_formatted
        (b, c) = block_size_and_index_chunks(
            argsb=args.block_size, argsc=args.index_chunks, memory=args.max_memory)
        run_diamond(
            args.input, f'{args.output}/aligned.blast', f'{args.output}/unaligned.blast', database,
            threads=args.threads, max_target_seqs=args.max_target_seqs, b=b, c=c, e_value=args.evalue,
            bit_score=args.bitscore, pident=args.pident, mode=args.diamond_mode)
        if args.max_target_seqs > 1:
            blast_consensus(f'{args.output}/aligned.blast').to_csv(
                f'{args.output}/consensus.blast', sep='\t', index=False)
        args.input = f'{args.output}/aligned.blast'
        args.blast = True

    if args.skip_id_mapping:
        exit('Not performing ID mapping as specified.')

    timed_message('ID mapping has begun.')
    args_input, input_type = get_input_type(args.input, blast=args.blast)

    # Get the IDs
    ids, full_id, sp_ids = get_ids(args_input, input_type=input_type, full_id=args.full_id)

    if args.output_table:
        table_output = args.output_table
        print(f'Overrided table output to {table_output}')
        Path('/'.join(args.output_table.split('/')[:-1])).mkdir(parents=True, exist_ok=True)
    else:
        table_output = f'{args.output}/uniprotinfo.tsv'

    if args.from_db != 'UniProtKB AC/ID' or args.to_db != 'UniProtKB':
        basic_idmapping_multiprocess(ids, table_output, args.from_db, args.to_db, threads=args.threads)
        return

    if not args.skip_id_checking:
        ids, not_valid = get_valid_entries_multiprocess(ids, threads=args.threads)
        # UniProt's API now fails if outdated IDs or entry names are submitted. This function removes those
        with open(f'{args.output}/valid_ids.txt', 'w') as f:
            f.write('\n'.join(ids))
        with open(f'{args.output}/not_valid_ids.txt', 'w') as f:
            f.write('\n'.join(not_valid))

    # Get UniProt information
    if not args.fasta:
        # ID mapping through local information
        if args.local_id_mapping:
            ids = set(ids) - set(local_id_mapping(
                sp_ids, f'{args.resources_directory}/uniprot_sprot.dat', f'{args.resources_directory}/taxonomy.tsv',
                table_output, columns=args.columns, databases=args.databases, threads=15))

        # ID mapping through API
        uniprot_information_workflow(
            ids, table_output, columns=args.columns, step=args.step, max_iter=args.max_tries,
            sleep_time=args.sleep)
        result = pd.read_csv(table_output, sep='\t', low_memory=False)

        if not args.no_annotation:
            blast = parse_blast(f'{args.output}/aligned.blast')
            if full_id:
                blast.sseqid = [ide.split('|')[1] if ide not in ['*', ''] else ide for ide in blast.sseqid]
            result = pd.merge(blast, result, left_on='sseqid', right_on='Entry')
        sort_columns = ['Entry'] if args.no_annotation else ['qseqid', 'evalue']
        result.sort_values(by=sort_columns, ascending=False).to_csv(
            f'{args.output}/UPIMAPI_results.tsv', index=False, sep='\t')
    else:
        uniprot_fasta_workflow(
            ids, f'{args.output}/uniprotinfo.fasta', step=args.step, sleep_time=args.sleep)


if __name__ == '__main__':
    start_time = time()
    upimapi()
    timed_message(f'UPIMAPI analysis finished in {human_time(time() - start_time)}')
