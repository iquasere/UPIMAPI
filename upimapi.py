"""
UPIMAPI - UniProt Id Mapping through API

By João Sequeira

Mar 2020
"""

from uniprot_support import UniprotSupport
from progressbar import ProgressBar
from io import StringIO
import pandas as pd
import argparse, time, os, urllib.request, urllib.parse, urllib.error

upmap = UniprotSupport()

class UPIMAPI:
    
    def __init__ (self, **kwargs):
        self.__dict__ = kwargs
        
    '''
    Input:
        filename: str - FASTA filename
    Output:
        dict {header : sequence}
    '''
    def parse_fasta(self, filename):
        lines = [line.rstrip('\n') for line in open(filename)]
        i = 0
        sequences = dict()
        while i < len(lines):
            if lines[i].startswith('>'):
                name = lines[i][1:]
                sequences[name] = ''
                i += 1
                while i < len(lines) and not lines[i].startswith('>'):
                    sequences[name] += lines[i]
                    i += 1
        return sequences
    
    '''
    Input:
        filename: str - BLAST filename
    Output:
        pandas.DataFrame
    '''
    def parse_blast(self, blast):
        result = pd.read_csv(blast, sep='\t', header = None)
        result.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                          'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                          'evalue', 'bitscore']
        return result
    
    '''
    Input:
        ids: list of UniProt IDs to query
        original_database: database from where the IDs are
        database_destination: database to where to map (so far, only works with 'ACC'
        output_format: format of response to get
        columns: names of UniProt columns to get info on
        databases: names of databases to cross-reference with
    Output:
        Returns the content of the response from UniProt
    '''
    def uniprot_request(self, ids, original_database = 'ACC+ID', database_destination = '',
                        output_format = 'tab', columns = list(), databases = list()):
        base_url = 'https://www.uniprot.org/uploadlists/'
        
        params = {
            'from':original_database,
            'format':output_format,
            'query':' '.join(ids),
            'columns':upmap.string4mapping(columns = columns, databases = databases)
        }
        
        if database_destination == '' or original_database == 'ACC+ID':
            params['to'] = 'ACC'
        
        data = urllib.parse.urlencode(params).encode("utf-8")
        request = urllib.request.Request(base_url, data)
        response = urllib.request.urlopen(request)
        return response.read().decode("utf-8")
    
    '''
    Input:
        ids: list of UniProt IDs to query
        original_database: database from where the IDs are
        database_destination: database to where to map (so far, only works with 'ACC'
        chunk: INT, number of IDs to send per request
        sleep: INT, number of seconds to wait between requests
        columns: list - names of UniProt columns to get info on
        databases: list - names of databases to cross-reference with
    Output:
        pd.DataFrame will be returned with the information about the IDs queried.
    '''
    def get_uniprot_information(self, ids, original_database = 'ACC+ID',
                                database_destination = '', step = 1000, sleep = 30,
                                columns = list(), databases = list()):
        pbar = ProgressBar()
        print('Retrieving UniProt information from ' + str(len(ids)) + ' IDs.')
        result = pd.DataFrame()
        for i in pbar(range(0, len(ids), int(step))):
            j = i + step if i + step < len(ids) else len(ids)
            try:
                data = self.uniprot_request(ids[i:j], original_database, 
                            database_destination, columns = columns, 
                            databases = databases)
                if len(data) > 0:
                    uniprotinfo = pd.read_csv(StringIO(data), sep = '\t')
                    k = (len(uniprotinfo.columns) - (30 if (len(columns + databases) == 0 
                        or (columns == [''] and databases == ['']))
                        else (len(columns) + len(databases))))                      # Removes the "yourlist:" and "isomap:" columns
                    uniprotinfo = uniprotinfo[uniprotinfo.columns.tolist()[:-k]]
                    result = pd.concat([result, uniprotinfo[uniprotinfo.columns.tolist()]])
                time.sleep(sleep)
            except:
                print('Mapping failed at some point!')
                return result
        return result

    '''
    Input:
        ids: list of UniProt IDs to query
        chunk: INT, number of IDs to send per request
        sleep: INT, number of seconds to wait between requests
    Output:
        str object containing the fasta sequences and headers
        of the proteis belonging to the IDs queried will be returned
    '''
    def get_uniprot_fasta(self, ids, step = '1000', sleep = 30):
        pbar = ProgressBar()
        print('Building FASTA from ' + str(len(ids)) + ' IDs.')
        result = str()
        for i in pbar(range(0, len(ids), int(step))):
            j = i + step if i + step < len(ids) else len(ids)
            data = self.uniprot_request(ids[i:j], original_database = 'ACC+ID', 
                        database_destination = '', output_format = 'fasta')
            if len(data) > 0:
                result += data
            time.sleep(sleep)
        return result
    
    def recursive_uniprot_fasta(self, all_ids, output, max_iter = 5, step = '1000'):
        if os.path.isfile(output):
            print(output + ' was found. Will perform mapping for the remaining IDs.')
            ids_done = self.parse_fasta(output).keys()
        else:
            print(output + ' not found. Will perform mapping for all IDs.')
            ids_done = list()
        tries = 0
        ids_unmapped_output = '/'.join(output.split('/')[:-1]) + '/ids_unmapped.txt'
        ids_missing = list(set(all_ids) - set(ids_done))
        
        tries = 0
        ids_done = ([ide.split('|')[1] for ide in self.parse_fasta(output).keys()]
                    if os.path.isfile(output) else list())
        while len(ids_done) < len(all_ids) and tries < max_iter:
            print('Checking which IDs are missing information.')
            pbar = ProgressBar()
            ids_missing = list(set([ide for ide in pbar(all_ids) if ide not in ids_done]))
            print('Information already gathered for ' + str(len(ids_done)) + 
                  ' ids. Still missing for ' + str(len(ids_missing)) + '.')
            uniprotinfo = self.get_uniprot_information(ids_missing, step = step,
                                                       output_format = 'fasta')
            with open(output, 'a') as file:
                file.write(uniprotinfo)
            ids_done = [ide.split('|')[1] for ide in self.parse_fasta(output).keys()]
            tries += 1
        if len(ids_done) == len(all_ids):
            print('Results for all IDs are available at ' + output)
        else:
            ids_unmapped_output = '/'.join(output.split('/')[:-1]) + '/ids_unmapped.txt'
            handler = open(ids_unmapped_output, 'w')
            handler.write('\n'.join(ids_missing))
            print('Maximum iterations were made. Results related to {} IDs were not obtained. IDs with missing information are available at {} and information obtained is available at {}'.format(
                    str(len(ids_missing)), ids_unmapped_output, output))
        
    def recursive_uniprot_information(self, ids, output, max_iter = 5, excel = False,
                                      columns = list(), databases = list(), step = '1000'):
        if os.path.isfile(output):
            print(output + ' was found. Will perform mapping for the remaining IDs.')
            result = pd.read_csv(output, sep = '\t', low_memory=False).drop_duplicates()
            ids_done = list(set(result['Entry'].tolist() + result['Entry name'].tolist()))
        else:
            print(output + ' not found. Will perform mapping for all IDs.')
            result = pd.DataFrame()
            ids_done = list()
        tries = 0
        ids_unmapped_output = '/'.join(output.split('/')[:-1]) + '/ids_unmapped.txt'
        ids_missing = list(set(ids) - set(ids_done))
        last_ids_missing = None
        
        print('IDs present in uniprotinfo file: ' + str(len(ids_done)))
        print('IDs missing: ' + str(len(ids_missing)))
        
        while len(ids_missing) > 0 and tries < max_iter and ids_missing != last_ids_missing:
            print('Information already gathered for {} ids. Still missing for {}.'.format(
                    str(len(ids_done)), str(len(ids_missing))))
            last_ids_missing = ids_missing
            uniprotinfo = self.get_uniprot_information(ids_missing, step = step,
                                columns = columns, databases = databases)
            if len(uniprotinfo) > 0:
                ids_done += list(set(uniprotinfo['Entry'].tolist() + uniprotinfo['Entry name'].tolist()))
                result = pd.concat([result, uniprotinfo], ignore_index = True)
            ids_missing = list(set(ids) - set(ids_done))
            
            if len(ids_missing) > 0:
                if last_ids_missing == ids_missing:
                    print('Could not map additional IDs for this mapping. There were probably some outdated IDs. For more questions, please contact through https://github.com/iquasere/UPIMAPI/issues'"")
                else:
                    print('Failed to retrieve information for some IDs. Retrying request.')
                    tries += 1
        
        if not excel:
            result.to_csv(output, sep = '\t', index = False)
        else:
            output += '.xlsx'
            result.to_excel(output, index = False)
        
        if len(ids_missing) == 0:
            print('Results for all IDs are available at ' + output)
        else:
            open(ids_unmapped_output, 'w').write('\n'.join(ids_missing))
            print('Maximum iterations were made. Results related to {} IDs were not obtained. IDs with missing information are available at {} and information obtained is available at {}'.format(
                    str(len(ids_missing)), ids_unmapped_output, output))
        
    '''
    Input:
    Output:
    '''
    def get_ids(self, inpute, blast = False, full_id = False):
        if blast:
            ids = self.parse_blast(inpute)['sseqid']
            ids = [ide for ide in ids if ide != '*']                            # removes the non identified
        else:
            ids = open(inpute).read().split('\n')
        if full_id:
            ids = [ide.split('|')[1] for ide in ids]
        return ids
    
    def upimapi(self):
        parser = argparse.ArgumentParser(description = "UniProt Id Mapping through API",
                             epilog = """A tool for retrieving information from UniProt.""")
        parser.add_argument("-i", "--input", required = True,
                            help = """Input filename - can be a list of IDs (one per line)
                            or a BLAST TSV file - if so, specify with the --blast parameter""")
        parser.add_argument("-o", "--output", help = "filename of output",
                            default = "./uniprotinfo.tsv")
        parser.add_argument("--excel", help = "Will produce output in EXCEL format (default is TSV)",
                            action = "store_true", default = False)
        parser.add_argument("-anncols", "--annotation-columns", default = '',
                            help = "List of UniProt columns to obtain information from")
        parser.add_argument("-anndbs", "--annotation-databases", default = '',
                            help = "List of databases to cross-check with UniProt information")
        parser.add_argument("--blast", help = "If input file is in BLAST TSV format",
                            action = "store_true", default = False)
        parser.add_argument("--full-id", help = "If IDs are in 'full' format: tr|XXX|XXX",
                            action = "store_true", default = False)
        parser.add_argument("--fasta", help = "Output will be generated in FASTA format",
                            action = "store_true", default = False)
        parser.add_argument("--step", help = "How many IDs to submit per request to the API (default is 1000)",
                            default = '1000')
        
        args = parser.parse_args()
        
        # Get the IDs
        ids = self.get_ids(args.input, blast = args.blast, full_id = args.full_id)
        
        # Get UniProt information
        if not args.fasta:
            columns = args.annotation_columns.split(',') if args.annotation_columns != '' else list()
            databases = args.annotation_databases.split(',') if args.annotation_databases != '' else list()
            
            self.recursive_uniprot_information(ids, args.output, columns = columns,
                                              databases = databases, excel = args.excel,
                                              step = args.step)
        else:
            self.recursive_uniprot_fasta(ids, args.output, step = args.step)
        
if __name__ == '__main__':
    
    UPIMAPI().upimapi()
    