import microbes_online as mo
import conserved_domains as cdd
import pandas as pd
import datetime as dt
import time
import sys

try:
    gene_file = sys.argv[1]
except:
    print('Need a file to read genes from.')
    sys.exit(1)

try:
    output_file = sys.argv[2]
except:
    print('Need an output file name.')

# Merge columns of information that must stay together when reducing
# information based on how many unique values there are
def merge_columns(df, column_names):
    for name in column_names.copy():
        try:
            temp = df[name]
        except KeyError:
            print('\tmerge_columns: {} not in dataframe, removing'.format(name))
            # prevent second error
            time.sleep(.1)
            column_names.remove(name)

    merged_column_name = '|||'.join(column_names)
    print(column_names)
    merged_columns = []

    # merge columns that have values and leaves ones that don't blank
    for i, row in df.iterrows():
        joined_row = ''
        for name in column_names:
            if row[name] is None:
                row[name] = ''
            else:
                row[name] = str(row[name])
        joined_row = '|||'.join(row[column_names])
        merged_columns.append(joined_row)

    # add new column with merged columns
    df[merged_column_name] = merged_columns

    # remove old columns
    df.drop(column_names, axis=1, inplace=True)
    return(df, column_names)


# make prevously merged columns into separate fields again
def unmerge_columns(df, column_names):
    merged_column_name = '|||'.join(column_names)
    new_columns = pd.DataFrame([x.split('|||') for x in df[merged_column_name]],
                               columns=column_names)
    df = pd.concat([df, new_columns], axis=1)
    df.drop(merged_column_name, axis=1, inplace=True)
    return(df)


def reshape_data(df):
    # merge columns where values must stay together
    GO_columns = ['go_id',
                  'go_name',
                  'go_type']
    interpro_columns = ['ipr_id',
                        'ipr_name']
    cdd_columns = ['accession','cdd_name',
                   'e-value', 'cdd_description']
    # convert E-Vlaue to strings
    df, go_names = merge_columns(df, GO_columns)
    df, ipr_names = merge_columns(df, interpro_columns)
    df, cdd_names = merge_columns(df, cdd_columns)

    # find unique list of genes
    genes = df['name'].unique()
    new_rows = []
    columns = df.columns.values

    # look through all genes and find uniqe values for each gene
    for gene in genes:
        # get all rows and convert NaN to empty strings
        rows = df[df['name'] == gene].fillna(value='')
        num_rows = 0
        for column in columns:
            # compare old columns to new ones and take maximum values
            num_unique_vals = len(set([x for x in rows[column]]))
            num_rows = max(num_rows, num_unique_vals)

        # build new part of data frame for gene
        for i in range(num_rows):
            row_values = []
            for column in columns:
                col_vals = rows[column].unique()
                # only include each vlaue once, unless it's the name
                if i < len(col_vals):
                    row_values.append(col_vals[i])
                elif column == 'name':
                    row_values.append(gene)
                else:
                    row_values.append('')
            new_rows.append(row_values)
    # write modifed data to df
    df = pd.DataFrame(new_rows, columns=columns)

    # make merged fields individual again
    df = unmerge_columns(df, go_names)
    df = unmerge_columns(df, ipr_names)
    df = unmerge_columns(df, cdd_names)

    return(df)


def create_interpro_link(df):
    # url for interpro website information
    url = 'https://www.ebi.ac.uk/interpro/entry/'

    # an empty list to hold hyperlinks
    interpro_links = []
    for ipr in df['ipr_id']:
        link = ''
        # check if there is an interpro id, if so make hyperlink
        if isinstance(ipr, str):
            link = '=HYPERLINK("' + url + ipr + '")'
        interpro_links.append(link)

    # add hyperlink list to df and return
    df['ipr_link'] = interpro_links
    return(df)


def read_tmp_csv(name):
    df = pd.read_csv('tmp/' + name, dtype=str)

if __name__ == '__main__':
    fields = ['name', 'locus_id', 'organism', 'gene_name', 'gene_description',
            'cdd_name', 'e-value', 'cdd_description', 'ipr_id', 'ipr_name',
            'ipr_link', 'fun_code', 'fun_code_description',
            'fun_code_group', 'cog_info_id', 'cog_description',
            'tigr_description', 'go_id', 'go_name', 'go_type',
            'ncbi_accession_number', 'ncbi_gene_id', 'gi', 'accession']

    start_time = time.time()
    print("Application started at {}".format(dt.datetime.now().time()))

    print('\n*** Microbes Online {}'.format(''.join(['*' for x in range(10)])))

    # 1. send fileds to microbes_online and get df with that information
    mo_df = mo.get_microbes_online_df(gene_file, fields)
    mo_time = time.time()
    print("Return microbes online dataframe: {0:.4f} sec\n".format(mo_time-start_time))

    #mo_df = pd.read_csv('tmp/tmp-microbes-online-laura-genes-1000.txt-2019-04-14 13:47:12.312666.csv', dtype=str)
    #mo_time = time.time()


    # 2. Create an interpo link from 'iprId'
    print('\n*** Creating interpro links {}'.format(''.join(['*' for x in range(10)])))
    mo_df = create_interpro_link(mo_df)
    ipr_time = time.time()
    print("Create interpro links: {0:.4f} sec\n".format(ipr_time - mo_time))


    # 3. Get conserved domain information using 'GI'
    print('\n*** Conserved Domains {}'.format(''.join(['*' for x in range(10)])))
    gi_list = mo_df['gi']
    cdd_df = cdd.get_cdd_information_from_gi(gi_list, fields)
    cdd_time = time.time()
    print("Get conserved domain info: {0:.4f} sec\n".format(cdd_time-ipr_time))

    #cdd_df = pd.read_csv('tmp/tmp-cdd-information-2019-04-14 13:49:20.305633.csv', dtype=str)

    # 3.5 Get cdd description using accession information
    print('\n*** Conserved Domains description {}'.format(''.join(['*' for x in range(10)])))
    acc_list = cdd_df['accession']
    desct_df = cdd.get_cdd_descriptions(acc_list)
    desct_time = time.time()
    print("Get conserved domain description: {0:.4f} sec\n".format(desct_time-cdd_time))


    # 4. join cdd and mo df
    print('\n*** Joining dataframes {}'.format(''.join(['*' for x in range(10)])))
    df = pd.merge(mo_df, cdd_df, how='left')
    df = pd.merge(df, desct_df, how='left')
    merge_time = time.time()
    print("Merge cdd and mo data frames: {0:.4f} sec\n".format(merge_time-desct_time))


    # 5. Reformat data to desired output
    print('\n*** Reformat dataframe {}'.format(''.join(['*' for x in range(10)])))
    df = reshape_data(df)
    reshape_time = time.time()
    print("Reshape data: {0:.4f} sec\n".format(reshape_time-merge_time))

    #TODO: rename lous_id to vimss id
    # 6. Write to csv
    print('\n*** Print to csv file {}'.format(''.join(['*' for x in range(10)])))
    # or compare fields to current columns
    df = df.reindex(columns=fields)
    df.to_csv(output_file, index=False)

    print("Application finished at: {}".format(dt.datetime.now()))
    print("Total runtime:: {0:.4f} sec".format((time.time()-start_time)))
