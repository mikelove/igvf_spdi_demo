# this script takes as input a variant list that looks like "1_25253604_hg38_G_A"
# easily modifiable for other specifications, adjust `split_char` and `info_idx`
# note that this script assumes input is 1-based position
# note skript expects csv and given column of the user string and adds SPDI column to the dataframe
# note if you run the script with call-spdi-batch flag it will expect read and write permissions in the current directory
# example call: python make_spdi_list.py --input-file=/test_data/user_variant_pattern_100.csv --output-file=test_data/user_variant_pattern_100_spdi.csv --column-separator='\t' --string-separator="_" --column-name="variant_string" --indices="0,1,3,4" --spdi-batch-processing-output="test_data/spdi_for_batch_processing.txt" --spdi-batch-output="test_data/spdi_batch_output.txt"

import click
import numpy as np
import pandas as pd
import subprocess # if call-spdi-batch flag is used

# No scientific notation
np.set_printoptions(suppress=True)

# dict of chr number to refseq chromosome number
chrom_refseq_dict = {
    "1": "NC_000001.11",
    "2": "NC_000002.12",
    "3": "NC_000003.12",
    "4": "NC_000004.12",
    "5": "NC_000005.10",
    "6": "NC_000006.12",
    "7": "NC_000007.14",
    "8": "NC_000008.11",
    "9": "NC_000009.12",
    "10": "NC_000010.11",
    "11": "NC_000011.10",
    "12": "NC_000012.12",
    "13": "NC_000013.11",
    "14": "NC_000014.9",
    "15": "NC_000015.10",
    "16": "NC_000016.10",
    "17": "NC_000017.11",
    "18": "NC_000018.10",
    "19": "NC_000019.10",
    "20": "NC_000020.11",
    "21": "NC_000021.9",
    "22": "NC_000022.11",
    "X": "NC_000023.11",
    "Y": "NC_000024.10"}

@click.command()
@click.option('--input-file', type=click.Path(exists=True), required=True)
@click.option('--output-file', type=click.Path(readable=True), required=True)
@click.option('--column-separator', default='-', help='column-separator: Character that separates the columns.')
@click.option('--string-separator', default='\t', help='Character that separates the fields (chromosome, position, reference and alternative).')
@click.option('--column-name', help='Column name of the user defined string including information like chromosome, position, reference and alternative seperated by user defined separator.')
@click.option('--indices', default='0,1,3,4', help='Comma-separated list of indices for chromosome, position, ref, and alt.')
@click.option('--call-spdi-batch', default=True, help='Check SPDI directly with spdi_batch.py (expects write permissions in the current directory).')
@click.option('--spdi-batch-path', default="spdi_batch.py", help='Path to the spdi_batch.py script (default you run the script within the ifvf_spdi_demo directory).')
@click.option('--spdi-batch-processing-output', default="test_data/spdi_for_batch_processing.txt", help='Path to the output file of the created SPDI column if you run with --call-spdi-batch flag.')
@click.option('--spdi-batch-output', default="test_data/spdi_batch_output.txt", help='Path to the output file of spdi_batch.py script.')

def cli(input_file, output_file, column_separator, string_separator, column_name, indices, call_spdi_batch, spdi_batch_path, spdi_batch_processing_output, spdi_batch_output):
    def variant_string_2_spdi_column(user_string, separator="-", indices=[0,1,2,3]):
        """
        Returns SPDI identifier for given variant
        SPDI: refseq_chromosome:pos:ref:alt
        @params: indices: list of indices for chromosome, position, ref and alt in user_string
        @params: separator: separator used in the string defined by the user to seperate chromosome, position, reference and alternative
        @params: user_string: user defined string of chromosome, position, reference and alternative and other information seperated by separator
        """
        if len(user_string.split(separator)) < 4:
            raise ValueError('Not enough indices given: Expected chrom, pos, ref, alt position in user_string')
        split_string = user_string.split(separator)
        zero_based_position = int(split_string[indices[1]]) # - 1 # (input: 1-based => 0-based)
        return f'{chrom_refseq_dict[split_string[indices[0]]]}:{zero_based_position}:{split_string[indices[2]]}:{split_string[indices[3]]}'


    def process_variants(input_file, column_separator, string_separator, column_name, indices):
        """
        Reads the user specified file and adds SPDI column to the dataframe assumes 1-based input
        """
        # read user specified file
        df = pd.read_csv(input_file, sep=column_separator, engine='python')
        # add SPDI column
        df['SPDI'] = df[column_name].apply(lambda x: variant_string_2_spdi_column(x, separator=string_separator, indices=indices))
        return df


    def call_spdi_batch(spdi_batch_path, spdi_batch_output):
        """
        Calls the skript spdi_batch.py with the generated SPDI file and write to a file
        example call: python spdi_batch.py -i spdi_100.txt -t SPDI | grep "NC_" > spdi_batch_output.txt (see https://github.com/mikelove/igvf_spdi_demo for more details)
        """
        batch_output = open(spdi_batch_output, "w")
        p1 = subprocess.Popen(
            ["python", spdi_batch_path, "-i", spdi_batch_processing_output, "-t", "SPDI"],
            stdout=subprocess.PIPE
        )
        p2 = subprocess.Popen(
            ["grep", "NC_"],
            stdin=p1.stdout,
            stdout=batch_output
        )
        # Close the pipe in the first process
        p1.stdout.close()
        p1.wait()
        p2.wait()
        batch_output.close()
        return True


    def check_spdi_batch(spdi_batch_output):
        """
        Check if the output from the spdi_batch script mentioned warnings and shows them to the user
        """
        # check if in the second column are warnings
        spdi_batch_output_df = pd.read_csv(spdi_batch_output, sep="\t", header=None)
        spdi_batch_output_df.columns = ['given_SPDI', 'canonical_SPDI']
        warnings = spdi_batch_output_df[spdi_batch_output_df["canonical_SPDI"].str.contains("WARNING")]
        if warnings.shape[0] > 0:
            print("Warnings in SPDI batch processing:")
            print(warnings)
            print("Please check your input additional informations can be found here: https://github.com/mikelove/igvf_spdi_demo")
            return False
        return True

    # process given comma-separated list
    processed_indices = [int(i) for i in indices.split(',')]
    # add SPDI column to the dataframe
    speedy_df = process_variants(input_file, column_separator, string_separator, column_name, processed_indices)

    # write to file
    speedy_df.to_csv(output_file, index=False)

    if call_spdi_batch:
        # write SPDI column to file for batch processing
        speedy_df[['SPDI']].to_csv(spdi_batch_processing_output, index=False)
        if call_spdi_batch(spdi_batch_path, spdi_batch_output):
            print('NOTE: spdi batch script was called.\nOutput will be checked...')
        if check_spdi_batch(spdi_batch_output):
            print('NOTE: No warnings detected using the spdi_batch script')
        print('NOTE: Removing Temporary Files')
        # remove temporary file spdi_batch_output
        subprocess.run(["rm", spdi_batch_output])
        subprocess.run(["rm", spdi_batch_processing_output])

if __name__ == '__main__':
    cli()

