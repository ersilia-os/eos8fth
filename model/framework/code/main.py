# imports
import os
import csv
import sys
import pandas as pd 
import shutil
import tempfile
from run_predictions import get_predictions, make_dictn, get_consensus
from rdkit import Chem

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# current file directory
root = os.path.dirname(os.path.abspath(__file__))

def combine_consensus_files(output_folder, combined_file):
    files = os.listdir(output_folder)
    consensus_files = [f for f in files if f.endswith('-consensus.csv')]

    combined_df = None

    for file in consensus_files:
        file_path = os.path.join(output_folder, file)
        df = pd.read_csv(file_path)
        consensus_column = df['Consensus']

        filename = file.replace('-consensus.csv', '')

        if combined_df is None:
            combined_df = pd.DataFrame({'Consensus': consensus_column})
        else:
            combined_df[filename] = consensus_column

    combined_df.to_csv(combined_file, index=False)


# my model
def my_model(csv_file):
    temp_results_folder = tempfile.mkdtemp()
    try:
        #run model
        get_predictions(temp_results_folder, csv_file)

        #adapt output
        combined_file = 'consensus_files.csv'
        combine_consensus_files(temp_results_folder, combined_file)
        return combined_file

    finally:
        if os.path.exists(temp_results_folder):
            shutil.rmtree(temp_results_folder)


# read SMILES from .csv file, assuming one column with header
with open(input_file, "r") as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    smiles_list = [r[0] for r in reader]

# run model

output_consensus = my_model(input_file)

# write output in a .csv file
with open(output_file, "w") as f:
    writer = csv.writer(f)
    writer.writerow(["value"])  # header

    #adapt output again
    consensus_df = pd.read_csv(consensus_file)
    outputs = consensus_df['Consensus'].tolist()

    for o in outputs:
        writer.writerow([o])
