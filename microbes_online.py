import mysql.connector
from mysql.connector.errors import Error as MySQLError
import pandas as pd
from datetime import datetime as dt
import os
import time


# TODO: check database for gene, if doesn't have the query database
def get_microbes_online_df(gene_file, fields, run_by_batch=250):
    # make tmp dir if it doesn't exist
    try:
        os.makedirs("tmp/")
    except FileExistsError:
        pass
    all_genes = file_as_list(gene_file)
    # empty df to hold all gene information
    temp_file_name = 'tmp/tmp-microbes-online-{}-{}.csv'.format(gene_file,
                                                                dt.now())
    # get information in chunks
    connection = get_mysql_connection()
    for i in range(0, len(all_genes), run_by_batch):
        genes = []
        # make sure genes gets correct indicies
        if i + run_by_batch < len(all_genes):
            genes = all_genes[i:(i+run_by_batch)]
        else:
            genes = all_genes[i:]

        print("\n* Microbes_online now processing {}-{} genes of {}".format(
              i+1, (i+len(genes)), len(all_genes)))

        # make query and run
        query = make_query(genes)
        query_result, field_names = run_query(connection, query)

        # return df with all possible result combinations
        raw_df = postprocess_query_result(query_result, field_names)

        # write to temporary csv file
        with open(temp_file_name, 'a') as f:
            if i == 0:
                raw_df.to_csv(f, index=False)
            else:
                raw_df.to_csv(f, index=False, header=False)

    connection.close()
    final_df = pd.read_csv(temp_file_name, dtype=str)
    final_df = filter_output(final_df, fields)
    return final_df


# function to read in files as lists
def file_as_list(file_name):
    with open(file_name) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    return content


# function to generate query
def make_query(genes):
    genes = ["'" + g + "'" for g in genes]
    genes = ", ".join(genes)
    query = " ".join(file_as_list('query.txt'))
    query = query.replace("INSERT_GENE_LIST", genes)
    return query


def get_mysql_connection():
    # catch System error 104: connection reset by peer
    success = False
    num_iter = 0
    max_connect_errors = 1  # found in database variable
    while not success:
        try:
            connection = mysql.connector.connect(
                host='pub.microbesonline.org',
                user='guest',
                password='guest',
                db='genomics')
            success = True
            print("Connected to database at: pub.microbesonline.org")
            return connection
        except MySQLError as e:
            num_iter = num_iter + 1
            sleep_time = 60 + 2**num_iter
            msg = 'microbes online error at {}. Waiting {} seconds before\
                   trying to reconnect.'.format(dt.now(), sleep_time)
            print('\n'.join([msg, str(e)]))
            with open('error-messages-microbes-online.txt', 'a') as f:
                f.write('\n'.join(['\n\n', msg, str(e)]))
            time.sleep(sleep_time)
            if num_iter >= max_connect_errors:
                raise RuntimeError('Too many iterations for microbes\
                                    online query.')


# Connect to the database and run query
def run_query(connection, query):
    # Run the query
    cursor = connection.cursor()
    cursor.execute(query)
    result = cursor.fetchall()
    field_names = [i[0] for i in cursor.description]
    return(result, field_names)


def extract_synonym(df, synonym):
    data = df[df['synonym_description'] == synonym][['name', 'synonym']]
    data = data.drop_duplicates()
    data = data.rename({'synonym': synonym}, axis=1)
    df = pd.merge(df, data)
    return(df)


def postprocess_query_result(result, field_names):
    print("Post-processing query output.")
    # clean up duplicate rows and write to CSV
    df = pd.DataFrame(result, columns=field_names)
    df = extract_synonym(df, 'GI')
    df = extract_synonym(df, 'NCBI accession number')
    df = extract_synonym(df, 'NCBI GeneID')
    df.drop('synonym', axis=1, inplace=True)
    df.drop('synonym_description', axis=1, inplace=True)
    df = df.rename(index=str, columns={'GI': 'gi',
                        'NCBI accession number': 'ncbi_accession_number',
                        'NCBI GeneId': 'ncbi_gene_id'})

    # return raw data frame, with all columns
    return(df)


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
            print('Microbes online cannot filter: {}'.format(
                  ', '.join(bad_filters)))
    else:
        print("Nothing to filter. Returning full results.")

    return df
