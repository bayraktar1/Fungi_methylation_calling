import pandas as pd
from collections import defaultdict
import pyarrow.feather as feather
from argparse import ArgumentParser

parser = ArgumentParser(
    description='Calculate weighted methylation score per window from bismark, megalodon and deepsignalplant output')
parser.add_argument('--source', type=str, required=True, help='Tool name: bismark, modkit, deepsignalplant')
parser.add_argument('--file', type=str, required=True, help='Path to input file')
parser.add_argument('--chromosome_sizes', required=True, type=str, help='Path to file with chromosomse sizes:chromosome size')
parser.add_argument('--window_size', required=True, type=int, help='Size of windows')
parser.add_argument('--output', type=str, required=True, help='Name of the output file')
parser.add_argument('--score_calculation', required=True, type=str, help='fraction_score, weighted_meth, fungi_weighted_meth')
parser.add_argument('--min_cov', required=True, type=int, help='Minimal coverage for base')
args = parser.parse_args()


def read_bismark(file_name):
    """
    Extracts relevant columns from Bismark output and puts them into a dataframe
    :param file_name: Name of the file
    :return: Dataframe
    """
    file_header = ['Chr', 'position', 'count_methylated', 'count_unmethylated']
    df = pd.read_csv(file_name, sep='\t', usecols=[0, 1, 3, 4],
                     names=file_header)
    df['Chr'] = df['Chr'].astype(str)
    return df


def read_modkit(file_name):
    """
    Extracts relevant columns from modkit and puts them into a dataframe
    :param file_name: Name of the file
    :return: Dataframe
    """
    file_header = ['Chr', 'position', 'coverage', 'methylated_percentage']
    df = pd.read_csv(file_name, sep='\t', usecols=[0, 1, 9, 10],
                     names=file_header)
    df['count_methylated'] = df['coverage'] / 100 * df['methylated_percentage']
    df['count_unmethylated'] = df['coverage'] - df['count_methylated']
    df.drop(columns=['coverage', 'methylated_percentage'], inplace=True)
    df['Chr'] = df['Chr'].astype(str)
    return df


def read_deepsignalplant(file_name):
    """
    Extracts relevant columns from Deepsignalplant and puts them into a dataframe
    :param file_name: Name of the file
    :return: Dataframe
    """
    file_header = ['Chr', 'position', 'count_methylated', 'count_unmethylated']
    df = pd.read_csv(file_name, sep='\t', usecols=[0, 1, 6, 7],
                     names=file_header)
    df['Chr'] = df['Chr'].astype(str)
    return df


def read_file(source, file_name):
    """
    Calls file reading functions
    :param source: Name of the tool that produced the file
    :param file_name: Name of the file
    :return: Dataframe
    """
    match source:
        case 'bismark':
            df = read_bismark(file_name)
            return df
        case 'modkit':
            df = read_modkit(file_name)
            return df
        case 'deepsignalplant':
            df = read_deepsignalplant(file_name)
            return df
        case _:
            raise ValueError(f'Unknown file format {source} exiting...')


def read_chromosome_file(file):
    """
    Reads the chromosome file
    :param file: chromosome file in format: chromosome length
    :return: dictionary with lengths if chromosomes
    """
    chrom_dict = {}
    with open(file) as f:
        for line in f:
            (key, val) = line.split()
            # chrom_dict[int(key)] = val
            chrom_dict[key] = val
    return chrom_dict


def fungi_calculate_weighted_score(single_window):
    """
    Calculated the weighted methylation score adjusted for fungi
    :param single_window: Dataframe
    :return: Score for given window
    """
    # this is the sum of all methylated READS in the window
    sum_c = single_window['count_methylated'].sum()
    # this is the sum of all unmethylated READS in the window
    sum_t = single_window['count_unmethylated'].sum()
    # This is the number of all READS in the window
    sum_c_t = sum_c + sum_t

    # print(f'NUMBER OF READS IN WINDOW {sum_c_t}')
    # print(f'METHYLATED READS {sum_c}')
    weighted_methyl_score = sum_c / sum_c_t if sum_c_t != 0 else 0

    # if sum_c_t == 0:
    #     weighted_methyl_score = 0
    #     print('sum is zero')
    # else:
    #     weighted_methyl_score = sum_c / sum_c_t
    #     print('NOT zero')

    return weighted_methyl_score


def calculate_weighted_score(single_window):
    """
    Calculated the weighted methylation score
    :param single_window: Dataframe
    :return: Score for given window
    """
    # this assumes bimodal distribution and high methylation site levels
    # we don't have this in fungi
    total_reads = single_window['count_methylated'].sum() + single_window['count_unmethylated'].sum()
    mod_bool_list = []
    for mod in single_window['count_methylated']:
        mod_bool = True if mod >= 3 else False
        mod_bool_list.append(mod_bool)
    single_window['mod'] = mod_bool_list

    # the total number of reads that support methylation on methylated reads
    # only
    total_methylated_reads_from_meth_C = 0
    for num, mod in zip(single_window['count_methylated'], single_window['mod']):
        if mod:
            total_methylated_reads_from_meth_C += num

    score = 0 if total_reads == 0 else total_methylated_reads_from_meth_C / total_reads

    return score


def calculate_fraction_score(single_window):
    """
    Calculates the fraction score
    :param single_window: Dataframe
    :return: Score for given window
    """
    # fractional methylation, calculated as the number of methylated
    # cytosines divided by all cytosine positions
    # For fractional methylation, a cytosine was considered methylated
    # if it was at least 5% methylated from all the reads covering that
    # cytosine or in our case > 3 methylated reads

    mod_mask = (single_window['count_methylated'] >= 3)
    total_c = len(single_window[mod_mask].index)
    mod_c = len(single_window[mod_mask].index)
    score = 0 if total_c == 0 else mod_c / total_c
    return score


def create_window(dataframe, chromosome, chrom_size, window_size):
    """
    Divides a chromosome into windows of x kb and calculates a score per window
    :param dataframe: Contains called cytosines on chromosome
    :param chromosome: Name of the chromosome
    :param chrom_size: Size of the chromosome
    :param window_size: Size of the window
    :return: Windows of chromosomes in a dataframe with scores
    """
    window_values = defaultdict(list)
    min_value = 1
    max_value = chrom_size

    num_windows = (max_value - min_value) // window_size + 1

    start_window = min_value
    end_window = start_window + window_size

    for i in range(num_windows):

        if end_window > chrom_size:
            end_window = chrom_size
        if start_window > end_window:
            continue

        range_mask = (dataframe['position'].values >= start_window) & \
                     (dataframe['position'].values < end_window)
        current_window_cytosines = dataframe.loc[range_mask]

        window_values['Number'].append(i)
        window_values['start'].append(start_window)
        window_values['end'].append(end_window)
        window_values['Chr'].append(chromosome)

        match args.score_calculation:
            case 'fraction_score':
                score = calculate_fraction_score(current_window_cytosines.copy())
            case 'weighted_meth':
                score = calculate_weighted_score(current_window_cytosines.copy())
            case 'fungi_weighted_meth':
                score = fungi_calculate_weighted_score(current_window_cytosines.copy())

        window_values['score'].append(score)

        start_window += window_size + 1
        end_window += window_size + 1

    return pd.DataFrame(window_values)


if __name__ == "__main__":
    df_raw = read_file(args.source, args.file)

    # not filtering on min methylated reads for fungi
    filter_mask = (
            df_raw['count_methylated'] + df_raw['count_unmethylated'] > args.min_cov)

    df_filtered = df_raw[filter_mask]

    chrom_sizes = read_chromosome_file(args.chromosome_sizes)
    chr_names = df_raw['Chr'].unique()
    chr_window_list = []

    for name in chr_names:
        chr_mask = (df_filtered['Chr'] == name)
        size = chrom_sizes[str(name)]
        windows_df = create_window(df_filtered[chr_mask], name,
                                   int(size), args.window_size)
        chr_window_list.append(windows_df)

    final_df = pd.concat(chr_window_list)
    final_df['file'] = args.file
    feather.write_feather(final_df, args.output)
