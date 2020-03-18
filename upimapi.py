"""
UPIMAPI - UniProt Id Mapping through API

By João Sequeira

Mar 2020
"""

from uniprot_support import UniprotSupport
from progressbar import ProgressBar
from io import StringIO
import pandas as pd
import numpy as np
import argparse, subprocess, time, os, sys, urllib.request, urllib.parse, urllib.error

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
                        output_format = 'tab', columns = None, databases = None):
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
        output_format: format of response to get
        database_destination: database to where to map (so far, only works with 'ACC'
        chunk: INT, number of IDs to send per request
        sleep: INT, number of seconds to wait between requests
        columns: names of UniProt columns to get info on
        databases: names of databases to cross-reference with
    Output:
        If output format is 'tab', a pd.DataFrame will be returned with the information
        about the IDs queried.
        If output format is 'fasta', a Str containing the fasta sequences and headers
        of the proteis belonging to the IDs queried will be returned.
    '''
    def get_uniprot_information(self, ids, original_database = 'ACC+ID', output_format = 'tab',
                                database_destination = '', chunk = 1000, sleep = 30,
                                columns = None, databases = None):
        pbar = ProgressBar()
        print('Retrieving UniProt information from ' + str(len(ids)) + ' IDs.')
        if output_format == 'tab':
            result = pd.DataFrame()
            for i in pbar(range(0, len(ids), chunk)):
                j = i + chunk if i + chunk < len(ids) else len(ids)
                try:
                    data = self.uniprot_request(ids[i:j], original_database, 
                                database_destination, columns = columns, 
                                databases = databases)
                    if len(data) > 0:
                        uniprotinfo = pd.read_csv(StringIO(data), sep = '\t')
                        result = pd.concat([result, uniprotinfo[uniprotinfo.columns.tolist()[:-1]]])    # last column is uniprot_list
                    time.sleep(sleep)
                except:
                    return result
        elif output_format == 'fasta':
            result = str()
            for i in pbar(range(0, len(ids), chunk)):
                j = i + chunk if i + chunk < len(ids) else len(ids)
                data = self.uniprot_request(ids[i:j], original_database, 
                            database_destination, output_format = output_format)
                if len(data) > 0:
                    result += data
                time.sleep(sleep)
        return result
    
    def recursive_uniprot_fasta(self, output, fasta = None, blast = None, 
                                entries = None, max_iter = 5):
        if fasta is not None:
            fasta = self.parse_fasta(fasta)
            all_ids = list(set(fasta.keys()))
        elif blast is not None:
            blast = self.parse_blast(blast)
            all_ids = list(set([ide.split('|')[1] if ide != '*' else ide 
                                for ide in blast.sseqid]))
        elif entries is not None:
            all_ids = list(set(entries))
        else:
            print('Must specify either fasta or blast!')
            return
        i = 0
        ids_done = ([ide.split('|')[1] for ide in self.parse_fasta(output).keys()]
                    if os.path.isfile(output) else list())
        while len(ids_done) < len(all_ids) and i < max_iter:
            print('Checking which IDs are missing information.')
            pbar = ProgressBar()
            ids_missing = list(set([ide for ide in pbar(all_ids) if ide not in ids_done]))
            print('Information already gathered for ' + str(len(ids_done)) + 
                  ' ids. Still missing for ' + str(len(ids_missing)) + '.')
            uniprotinfo = self.get_uniprot_information(ids_missing, 
                                                       output_format = 'fasta')
            with open(output, 'a') as file:
                file.write(uniprotinfo)
            ids_done = [ide.split('|')[1] for ide in self.parse_fasta(output).keys()]
            i += 1
        if len(ids_done) == len(all_ids):
            print('Results for all IDs are available at ' + output)
        else:
            ids_unmapped_output = '/'.join(output.split('/')[:-1]) + '/ids_unmapped.txt'
            handler = open(ids_unmapped_output, 'w')
            handler.write('\n'.join(ids_missing))
            print(str(i) + ' iterations were made. Results related to ' + str(len(ids_missing)) + 
                  ' IDs were not obtained. IDs with missing information are available' +
                  ' at ' + ids_unmapped_output + ' and information obtained is available' +
                  ' at ' + output)
        
    def recursive_uniprot_information(self, ids, output, max_iter = 5,
                                      columns = None, databases = None):
        if os.path.isfile(output):
            print(output + ' was found. Will perform mapping for the remaining IDs.')
            result = pd.read_csv(output, sep = '\t', low_memory=False).drop_duplicates()
            ids_done = list(set(result['Entry']))
        else:
            print(output + ' not found. Will perform mapping for all IDs.')
            result = pd.DataFrame()
            ids_done = list()
        tries = 0
        ids_unmapped_output = '/'.join(output.split('/')[:-1]) + '/ids_unmapped.txt'
        ids_missing = list(set(ids) - set(ids_done))
        
        print('IDs present in blast file: ' + str(len(ids)))
        print('IDs present in uniprotinfo file: ' + str(len(ids_done)))
        print('IDs missing: ' + str(len(ids_missing)))
        
        while len(ids_missing) > 0 and tries < max_iter:
            print('Information already gathered for {} ids. Still missing for {}.'.format(
                    str(len(ids_done)), str(len(ids_missing))))
            uniprotinfo = self.get_uniprot_information(ids_missing,
                                columns = columns, databases = databases)
            ids_done += list(set(uniprotinfo['Entry']))
            result = pd.concat([result, uniprotinfo])
            ids_missing = list(set(ids) - set(ids_done))
            if len(ids_missing) > 0:
                print('Failed to retrieve information for some IDs. Retrying request.')
                tries += 1
            
        result.to_csv(output, sep = '\t', index = False)
        
        if len(ids_missing) == 0:
            print('Results for all IDs are available at ' + output)
        else:
            open(ids_unmapped_output, 'w').write('\n'.join(ids_missing))
            print('Maximum iterations were made. Results related to ' + str(len(ids_missing)) + 
                  ' IDs were not obtained. IDs with missing information are available' +
                  ' at ' + ids_unmapped_output + ' and information obtained is available' +
                  ' at ' + output)
        
    '''
    Input:
        pathway: a String row from UniProt's 'Pathway' column
    Output:
        returns List of pathways included in that row
    '''
    def split(self, pathway):
        pathway = pathway.split('. ')
        return [path for path in pathway if path != '']
    
    '''
    Input:
        ec: a String row from UniProt's 'EC number' column
    Output:
        returns List of EC numbers included in that row
    '''
    def split_ec(self, ec):
        ec = ec.split('; ')
        return [ec_number for ec_number in ec if ec_number != '']
    
    '''
    Input:
        pathway: a String row from UniProt's 'Pathway' column
    Output:
        Reeives a pd.DataFrame object and breaks the list elements in the 'Pathway'
        column, multiplicating the rows with several elements of 'Pathway' into
        individual 'Pathway' elements, each with its row
    '''    
    def using_repeat(self, df, column = 'Pathway'):
        lens = [len(item) for item in df[column]]
        dictionary = dict()
        for column in df.columns:
            dictionary[column] = np.repeat(df[column].values,lens)
        dictionary[column] = np.concatenate(df[column].values)
        return pd.DataFrame(dictionary) 
    
    '''
    Input:
        uniprotinfo: information from UniProt ID mapping
        blast: blast file from DIAMOND annotation
        output: basename for EXCEL files to output
    Output:
        Two EXCEL files formated for Krona plotting with taxonomic and functional 
        information.
        This function is very useful if wanting to use UniProt 'Pathway' information,
        as it uses the three previous functions to organize the information from
        that column into the three functional levels of UniProt Pathways.
        This function was very cool
    '''
    def uniprotinfo_to_excel(self, uniprotinfo, blast, output):
        blast = self.parse_blast(blast)
        uniprotdf = pd.read_csv(uniprotinfo, sep = '\t', index_col = 0).drop_duplicates()
        pbar = ProgressBar()
        blast['Coverage'] = [float(ide.split('_')[5]) for ide in pbar(blast.qseqid)]
        pbar = ProgressBar()
        blast['ID'] = [ide.split('|')[-1] for ide in pbar(blast.sseqid)]
        uniprotdf = pd.merge(uniprotdf, blast, left_index = True, right_on = 'ID')
        tax_columns = ['Taxonomic lineage (SUPERKINGDOM)','Taxonomic lineage (PHYLUM)',
                       'Taxonomic lineage (CLASS)','Taxonomic lineage (ORDER)',
                       'Taxonomic lineage (FAMILY)','Taxonomic lineage (GENUS)',
                       'Taxonomic lineage (SPECIES)','Coverage']
        taxdf = uniprotdf[tax_columns].groupby(tax_columns[:-1])['Coverage'].sum().reset_index()
        taxdf.to_excel(output + '_taxonomic_krona.xlsx', index = False)
        print('Saved taxonomy')
        func_columns = ['Pathway','Protein names','EC number']
        funcdf = uniprotdf[uniprotdf.Pathway.notnull()][func_columns + ['Coverage']]
        funcdf.Pathway = funcdf.Pathway.apply(self.split)
        funcdf = self.using_repeat(funcdf)
        pathways = pd.DataFrame([(path.split('; ') + [np.nan] * (3 - len(path.split('; ')))) 
                                    for path in funcdf.Pathway], index = funcdf.index)
        pathways.columns = ['Superpathway','Pathway','Subpathway']
        del funcdf['Pathway']
        funcdf = pd.concat([pathways, funcdf], axis = 1)
        funcdf = funcdf[['Superpathway','Pathway','Subpathway','EC number',
                         'Protein names','Coverage']]
        funcdf = funcdf.groupby(funcdf.columns.tolist()[:-1])['Coverage'].sum().reset_index()
        funcdf.to_excel(output + '_functional_krona.xlsx', index = False)
        print('Saved pathways')
           
    '''
    Input:
        tsv: filename of TSV file to be inputed. Must have the format 
        value\tcategorie1\tcategorie2\t..., with no header
        output: filename of HTML krona plot to output
    Output:
        A krona plot will be created at output if it has been specified, else
        at tsv.replace('.tsv','.html')
    '''
    def create_krona_plot(self, tsv, output = None):
        if output is None:
            output = tsv.replace('.tsv','.html')
        bashCommand = 'perl Krona/KronaTools/scripts/ImportText.pl {} -o {}'.format(tsv, output)
        subprocess.run(bashCommand.split(), stdout=sys.stdout, check = True)
        
    '''
    Input:
    Output:
    '''
    def get_ids(self, inpute, blast = False, entry_name = False):
        if blast:
            ids = upimapi.parse_blast(inpute)['sseqid']
        else:
            ids = open(inpute).read().split('\n')
        if entry_name:
            ids = [ide.split('|')[1] for ide in ids]
        return ids
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description = "UniProt Id Mapping through API",
                             epilog = """A tool for retrieving information from UniProt.""")
    parser.add_argument("-i", "--input", required = True,
                        help = """Input filename - can be a list of IDs (one per line)
                        or a BLAST TSV file - if so, specify with the --blast parameter""")
    parser.add_argument("-o", "--output", help = "filename of output",
                        default = "./uniprotinfo.tsv")
    parser.add_argument("--excel", help = "Will produce output in EXCEL format (default is TSV)",
                        action = "store_true", default = False)
    parser.add_argument("-anncols", "--annotation-columns", default = None,
                        help = "List of UniProt columns to obtain information from")
    parser.add_argument("-anndbs", "--annotation-databases", default = None,
                        help = "List of databases to cross-check with UniProt information")
    parser.add_argument("--blast", help = "If input file is in BLAST TSV format",
                        action = "store_true", default = False)
    parser.add_argument("--entry_name", help = "If IDs are in 'Entry name' format: tr|XXX|XXX",
                        action = "store_true", default = False)
    parser.add_argument("--fasta", help = "Output will be generated in FASTA format",
                        action = "store_true", default = False)
    
    args = parser.parse_args()
    upimapi = UPIMAPI()    
    
    # Get the IDs
    ids = upimapi.get_ids(args.input, blast = args.blast, entry_name = args.entry_name)
    print(ids)
    ids = [ide for ide in ids if ide != '*']                                    # removes the non identified that happen if input is blast
    
    # Get UniProt information
    upimapi.recursive_uniprot_information(ids, args.output, columns = args.annotation_columns,
                                          databases = args.annotation_databases)
    