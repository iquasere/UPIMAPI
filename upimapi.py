#!/usr/bin/env python
"""
UPIMAPI - UniProt Id Mapping through API

By João Sequeira

Mar 2020
"""

import argparse
import os
from time import strftime, gmtime, time, sleep
import urllib.error
import urllib.parse
import urllib.request
import subprocess
import psutil
import pathlib
import multiprocessing
from io import StringIO

import pandas as pd
from tqdm import tqdm

from uniprot_support import UniprotSupport

__version__ = '1.4.0'

upmap = UniprotSupport()


class UPIMAPI:

    def __init__(self, **kwargs):
        self.__dict__ = kwargs

    def get_arguments(self):
        parser = argparse.ArgumentParser(description="UniProt Id Mapping through API",
                                         epilog="A tool for retrieving information from UniProt.")
        parser.add_argument("-i", "--input", help="""Input filename - can be:
                                1. a file containing a list of IDs (one per line)
                                2. a BLAST TSV result file (requires to be specified with the --blast parameter
                                3. a protein FASTA file to be annotated (requires the --use-diamond and -db parameters)
                                4. nothing! If so, will read input from command line, and parse as CSV (id1,id2,...)""")
        parser.add_argument("-o", "--output", help="Folder to store outputs", default="UPIMAPI_output")
        parser.add_argument("-ot", "--output-table",
                            help="Filename of table output, where UniProt info is stored. If set, will override "
                                 "'output' parameter just for that specific file")
        parser.add_argument("-cols", "--columns", default='',
                            help="List of UniProt columns to obtain information from (separated by &)")
        parser.add_argument("-dbs", "--databases", default='',
                            help="List of databases to cross-check with UniProt information (separated by &)")
        parser.add_argument("--blast", action="store_true", default=False,
                            help="If input file is in BLAST TSV format (will consider one ID per line if not set)")
        parser.add_argument("--full-id", type=self.str2bool, default="auto",
                            help="If IDs in database are in 'full' format: tr|XXX|XXX")
        parser.add_argument("--fasta", help="Output will be generated in FASTA format", action="store_true",
                            default=False)
        parser.add_argument("--step", type=int, default=1000,
                            help="How many IDs to submit per request to the API [1000]")
        parser.add_argument("--max-tries", default=3, type=int,
                            help="How many times to try obtaining information from UniProt before giving up")
        parser.add_argument("--sleep", default=10, type=int, help="Time between requests (in seconds) [10]")
        parser.add_argument('-v', '--version', action='version', version=f'UPIMAPI {__version__}')

        diamond_args = parser.add_argument_group('DIAMOND arguments')
        diamond_args.add_argument("--use-diamond", action="store_true", default=False,
                                  help='''Use DIAMOND to annotate sequences before mapping IDs. Requires protein FASTA 
                                  files as input for "-db" and "-i" parameters''')
        diamond_args.add_argument("-db", "--database", default=None,
                                  help="""Reference database for annotation with DIAMOND. NOTICE: if database's IDs are 
                                  in 'full' format (tr|XXX|XXX), specify with ""--full-id" parameter.""")
        diamond_args.add_argument("-t", "--threads", type=int, default=multiprocessing.cpu_count() - 2,
                                  help="Number of threads to use in annotation steps [total - 2]")
        diamond_args.add_argument("--evalue", type=float, default=1e-3,
                                  help="Maximum e-value to report annotations for [1e-3]")
        diamond_args.add_argument("--pident", type=float, default=None,
                                  help="Minimum pident to report annotations for.")
        diamond_args.add_argument("--bitscore", type=float, default=None,
                                  help="Minimum bit score to report annotations for (overrides e-value).")
        diamond_args.add_argument("-mts", "--max-target-seqs", default=1,
                                  help="Number of annotations to output per sequence inputed [1]")
        diamond_args.add_argument("-b", "--block-size",
                                  help="Billions of sequence letters to be processed at a time "
                                       "(default: auto determine best value)")
        diamond_args.add_argument("-c", "--index-chunks",
                                  help="Number of chunks for processing the seed index "
                                       "(default: auto determine best value)")

        args = parser.parse_args()
        args.output = args.output.rstrip('/')

        return args

    def str2bool(self, v):
        if v.lower() == 'auto':
            return 'auto'
        elif v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    def get_fasta_ids(self, filename):
        return [line[1:-1] for line in open(filename) if line.startswith('>')]

    def parse_blast(self, blast):
        result = pd.read_csv(blast, sep='\t', header=None)
        result.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                          'gapopen', 'qstart', 'qend', 'sstart', 'send',
                          'evalue', 'bitscore']
        return result

    def uniprot_request(self, ids, original_database='ACC+ID', database_destination='',
                        output_format='tab', columns=list(), databases=list()):
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

    def get_uniprot_information(self, ids, original_database='ACC+ID',
                                database_destination='', step=1000, sleep_time=30,
                                columns=list(), databases=list(), max_tries=3):
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

        print('Retrieving UniProt information from ' + str(len(ids)) + ' IDs.')
        result = pd.DataFrame()
        for i in tqdm(range(0, len(ids), step), desc="UniProt ID mapping"):
            tries = 0
            done = False
            j = min(i + step, len(ids))
            while not done and tries < max_tries:
                try:
                    data = self.uniprot_request(ids[i:j], original_database,
                                                database_destination, columns=columns,
                                                databases=databases)
                    if len(data) > 0:
                        uniprotinfo = pd.read_csv(StringIO(data), sep='\t')
                        k = (len(uniprotinfo.columns) - (30 if (
                                len(columns + databases) == 0 or (columns == [''] and databases == [''])
                        ) else (len(columns) + len(databases))))  # Removes the "yourlist:" and "isomap:" columns
                        uniprotinfo = uniprotinfo[uniprotinfo.columns.tolist()[:-k]]
                        result = pd.concat([result, uniprotinfo[uniprotinfo.columns.tolist()]])
                    sleep(sleep_time)
                    done = True
                except:
                    print(f'ID mapping failed. Remaining tries: {max_tries - tries}')
                    tries += 1
                    sleep(10)
        return result

    def get_uniprot_fasta(self, ids, step=1000, sleep_time=30):
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
            data = self.uniprot_request(ids[i:j], original_database='ACC+ID',
                                        database_destination='', output_format='fasta')
            if len(data) > 0:
                result += data
            sleep(sleep_time)
        return result

    def uniprot_fasta_workflow(self, all_ids, output, max_iter=5, step=1000, sleep_time=10):
        if os.path.isfile(output):
            print(f'{output} was found. Will perform mapping for the remaining IDs.')
            ids_done = self.get_fasta_ids(output)
        else:
            print(f'{output} not found. Will perform mapping for all IDs.')
            ids_done = list()

        ids_missing = list(set(all_ids) - set(ids_done))

        tries = 0
        ids_done = ([ide.split('|')[1] for ide in self.get_fasta_ids(output)] if os.path.isfile(output) else list())
        while len(ids_done) < len(all_ids) and tries < max_iter:
            print('Checking which IDs are missing information.')
            ids_missing = list(set([ide for ide in tqdm(all_ids, desc='Checking which IDs are missing information.')
                                    if ide not in ids_done]))
            print(f'Information already gathered for {len(ids_done)} ids. Still missing for {len(ids_missing)}.')
            uniprotinfo = self.get_uniprot_fasta(ids_missing, step=step, sleep_time=sleep_time)
            with open(output, 'a') as file:
                file.write(uniprotinfo)
            ids_done = [ide.split('|')[1] for ide in self.get_fasta_ids(output)]
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

    def uniprot_information_workflow(self, ids, output, max_iter=5, columns=list(), databases=list(), step=1000,
                                      sleep_time=10):
        if os.path.isfile(output) and os.stat(output).st_size > 1:
            try:
                result = pd.read_csv(output, sep='\t', low_memory=False).drop_duplicates()
                print(f'{output} was found. Will perform mapping for the remaining IDs.')
                ids_done = list(set(result['Entry'].tolist() + result['Entry name'].tolist()))
            except:
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

        print(f'IDs present in uniprotinfo file: {str(len(ids_done))}')
        print(f'IDs missing: {str(len(ids_missing))}')

        while len(ids_missing) > 0 and tries < max_iter and ids_missing != last_ids_missing:
            print(f'Information already gathered for {str(len(ids_done))} ids. '
                  f'Still missing for {str(len(ids_missing))}.')
            last_ids_missing = ids_missing
            uniprotinfo = self.get_uniprot_information(
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

    def determine_full_id(self, ids):
        for ide in ids:
            if '|' in ide:
                return True
        return False

    def get_ids(self, inpute, input_type='blast', full_id='auto'):
        if input_type == 'blast':
            ids = self.parse_blast(inpute)['sseqid']
        elif input_type == 'txt':
            ids = open(inpute).read().split('\n')
        else:
            ids = inpute.split(',')
        if full_id == 'auto':
            full_id = self.determine_full_id(ids)
            print(f'Auto determined "full id" as: {full_id}')
        if full_id:
            return [ide.split('|')[1] for ide in ids if ide != '*'], full_id
        return [ide for ide in ids if ide != '*'], full_id

    def run_command(self, bash_command, print_message=True):
        if print_message:
            print(bash_command)
        subprocess.run(bash_command.split(), check=True)

    def generate_diamond_database(self, fasta, dmnd):
        self.run_command(f'diamond makedb --in {fasta} -d {dmnd}')

    def b_n_c(self, argsb, argsc):
        if argsb:
            b = argsb
        else:
            b = psutil.virtual_memory().available / (1024.0 ** 3) / 20  # b = memory in Gb / 20
        if argsc:
            return b, argsc
        if b > 3:
            return b, 1
        if b > 2:
            return b, 2
        if b > 1:
            return b, 3
        return b, 4

    def run_diamond(self, query, aligned, unaligned, database, threads=12, max_target_seqs=50, b=1, c=4, e_value=0.01,
                    bit_score=None, pident=None):
        command = (
            f"diamond blastp --query {query} --out {aligned} --un {unaligned} --db {database} --outfmt 6 --unal 1 "
            f"--threads {threads} --max-target-seqs {max_target_seqs} -b {b} -c {c} --evalue {e_value} --very-sensitive"
        )
        if bit_score:
            command += f' --min-score {bit_score}'
        if pident:
            command += f' --id {pident}'
        self.run_command(command)

    def upimapi(self):
        args = self.get_arguments()
        pathlib.Path(args.output).mkdir(parents=True, exist_ok=True)

        # Using annotation with DIAMOND
        if args.use_diamond:
            if not args.database.endswith(".dmnd"):
                self.generate_diamond_database(args.database, f"{'.'.join(args.database.split('.')[:-1])}.dmnd")
                args.database = f"{'.'.join(args.database.split('.')[:-1])}.dmnd"
            (b, c) = self.b_n_c(argsb=args.block_size, argsc=args.index_chunks)
            self.run_diamond(
                args.input, f'{args.output}/aligned.blast', f'{args.output}/unaligned.blast',
                args.database, threads=args.threads, max_target_seqs=args.max_target_seqs, b=b, c=c,
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
        ids, full_id = self.get_ids(args.input, input_type=input_type, full_id=args.full_id)

        # Get UniProt information
        if not args.fasta:
            columns = args.columns.split('&') if args.columns != '' else list()
            databases = args.databases.split('&') if args.databases != '' else list()

            if args.output_table:
                table_output = args.output_table
                pathlib.Path('/'.join(args.output_table.split('/')[:-1])).mkdir(parents=True, exist_ok=True)
            else:
                table_output = f'{args.output}/uniprotinfo.tsv'
            print(f'Overrided table output to {table_output}')
            self.uniprot_information_workflow(
                ids, table_output, columns=columns, databases=databases, step=args.step, max_iter=args.max_tries,
                sleep_time=args.sleep)

            if args.use_diamond:
                blast = self.parse_blast(f'{args.output}/aligned.blast')
                if full_id:
                    blast.sseqid = [ide.split('|')[1] if ide != '*' else ide for ide in blast.sseqid]
                result = pd.merge(blast, pd.read_csv(table_output, sep='\t'), left_on='sseqid', right_on='Entry')
                result.to_excel(f'{args.output}/UPIMAPI_results.xlsx', index=False)
                result.to_csv(f'{args.output}/UPIMAPI_results.tsv', index=False, sep='\t')
        else:
            self.uniprot_fasta_workflow(ids, f'{args.output}/uniprotinfo.fasta', step=args.step, sleep_time=args.sleep)


if __name__ == '__main__':
    start_time = time()
    UPIMAPI().upimapi()
    print(f'UPIMAPI analysis finished in {strftime("%Hh%Mm%Ss", gmtime(time() - start_time))}')
