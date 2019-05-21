import pandas as pd
import time
import urllib.request
from bs4 import BeautifulSoup
import sys
import re  # regular expression package
from datetime import datetime as dt
import os
from multiprocessing import Pool
if sys.version_info[0] < 3:
    from StringIO import StringIO
else:
    from io import StringIO


# get accession numbers and more for conserved domains from NCBI CDD
def get_cdd_information_from_gi(gi_list, filters=[], run_by_batch=250):
    # make temporary file to hold sub runs
    temp_file_name = 'tmp/tmp-cdd-information-{}.csv'.format(dt.now())
    # make tmp dir if it doesn't exist
    try:
        os.makedirs("tmp/")
    except FileExistsError:
        pass
    # if not a list, typecast single one to list
    if isinstance(gi_list, int) or isinstance(gi_list, str):
        gi_list = [gi_list]
    # make sure it isn't a pandas series
    elif not isinstance(gi_list, list):
        gi_list = list(gi_list)


    # find unique values
    gi_list = list(set(gi_list))

    # Check if results have been returned or if None returned every time
    results_returned = False
    # get information in chunks
    for i in range(0, len(gi_list), run_by_batch):
        sub_gis = []
        # make sure genes gets correct indicies
        if i + run_by_batch < len(gi_list):
            sub_gis = gi_list[i:(i+run_by_batch)]
        else:
            sub_gis = gi_list[i:]

        print("\n* Conserved_domains now processing {}-{} GIs of {}".format(
              i+1, (i+len(sub_gis)), len(gi_list)))

        query_results = query_ncbi_cdd(sub_gis)

        # check if something was returned, if not return None
        if query_results is not None:
            query_results = add_gi_to_ncbi_query_results(query_results)

            # accession id portion
            df_acc_id = query_results[['gi', 'Accession', 'Short name',  'E-Value']]
            df_acc_id = df_acc_id.rename(index=str, columns={'Short name': 'cdd_name',
                                                             'Accession': 'accession'})  # , inplace=True)
            # print(df_acc_id)

            # write to temporary csv file
            with open(temp_file_name, 'a') as f:
                if not results_returned:
                    df_acc_id.to_csv(f, index=False)
                else:
                    df_acc_id.to_csv(f, index=False, header=False)
            results_returned = True

    # TODO: Change to file checking and returning empty df
    if not results_returned:
        print('Conserved domains found no information for list of GIs')
        return None

    # read in final df and return it
    final_df = pd.read_csv(temp_file_name, dtype=str)
    # filter output if user gave filters
    if len(filters) > 0:
        final_df = filter_output(final_df, filters)
    return final_df


# Returns the resulting file from ncbi cdd query
# TODO: Return error if gi_list is longer than 4000
def query_ncbi_cdd(gi_list):
    # if not a list, typecast single one to list
    if not isinstance(gi_list, list):
        gi_list = [gi_list]

    # gi is joined string
    gi = '%0A'.join([str(x) for x in gi_list])

    # amount of time to wait before moving on
    max_wait_time = 20 + 2*len(gi_list)

    # get the id to search if query has finished yet
    url = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?queries="\
          + gi + "&useid1=true&tdata=hits"
    success = False
    try:
        contents = urllib.request.urlopen(url).read().decode("utf-8")

        # get the query id from it's location on the content returned
        query_id = contents.split('cdsid\t')[1].split('\n')[0]
        # print('\nquery_id', query_id)
        url_for_checking = 'https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb'\
                           '.cgi?cdsid=' + query_id
        # check the url eveery 2 seconds until there is a success
        for t in range(int(max_wait_time/2)):
            time.sleep(2)
            cdd_results = urllib.request.urlopen(url_for_checking).read().decode('utf-8')
            if 'success' in cdd_results:
                success = True
                print('NCBI query of {} returned in {} seconds'.format(len(gi_list),
                                                                        t*2))
                break
    # if the web page doesn't work, report what value for
    except urllib.error.URLError:
        print("URLError when downloading query for GI: {}".format(gi))
        return(None)

    # the job didn't finish in time, print ouput
    if success == False:
        print("For GI {}, query did not return after{} seconds. Aborting."\
              "".format(max_wait_time, gi_list))
        # return empty list so not a nontype
        return(None)

    # make pandas table from the resulting string
    df = pd.read_csv(StringIO(cdd_results.split("\n\n")[1]), sep="\t")

    return df


# get the GI number from an ncbi 'Query' column
def add_gi_to_ncbi_query_results(ncbi_results):
    # example query results
    # Q#1 - 116863(Warning: this sequence record may be obsolete or preliminary)
    # Q#1 - 116863

    gi_values = []
    for line in ncbi_results['Query']:
        gi = line.split(' ')[2]
        gi = gi.split('(')[0]
        gi_values.append(gi)

    ncbi_results['gi'] = gi_values

    return ncbi_results


def get_cdd_descriptions(accession_list, num_processes=5):
    # if not a list, typecast single one to list
    if isinstance(accession_list, int) or isinstance(accession_list, str):
        accession_list = [accession_list]
    # make sure it isn't a pandas series
    elif not isinstance(accession_list, list):
        accession_list = list(accession_list)

    accession_list = list(set(accession_list))
    start_time = time.time()

    # divide into three pocesses
    # num_processes = 5
    desct_list = []
    with Pool(num_processes) as p:
        desct_list = p.map(get_single_cdd_description, accession_list)

    print('Processes {} run time {}'.format(num_processes, time.time()-start_time))
    print('len desct list: ', len(desct_list), ' len acc: ', len(accession_list))

    desct_df = pd.DataFrame({'accession': accession_list, 'cdd_description':
                             desct_list})

    return desct_df


# Uses accession number from conserved domains to get description from
# webpage
def get_single_cdd_description(cdd_accession):
    description = ""
    # get webpage using accession id for conserved domain
    url = 'https://www.ncbi.nlm.nih.gov/Structure/cdd/' + cdd_accession
    try:
        webpage = BeautifulSoup(urllib.request.urlopen(url).read(),
                                features='html.parser')


        # find the description based on looking at layout of html
        # may not work in all cases
        desct = webpage.find("div", {"id": "dscpt"})
        if desct is not None:
            span_descript = str(desct.find_all("span")[-1])

            # remove tags <something here> using regular expressions
            description = re.sub('<.*?>', ' ', span_descript).strip()
        else:
            description = ''

    # catch if problem with the webpage
    except urllib.error.URLError:
        print("URLError for cdd accession website: {}".format(cdd_accession))
    except ConnectionResetError as e:
        print(e)

    return description


# Select columns in filters list and removes others
def filter_output(df, filters):
    current_columns = list(df.columns.values)
    if len(filters) > 0:
        bad_filters = []
        # find unwanted columns and drop them
        for new in filters:
            if new in current_columns:
                current_columns.remove(new)
            else:
                bad_filters.append(new)
        for col in current_columns:
            df.drop(col, axis=1, inplace=True)
        if len(bad_filters) > 0:
            print('Not known options to filter: {}'.format(', '.join(bad_filters)))

    return df
