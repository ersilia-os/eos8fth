# imports
import os
import csv
import sys
import pandas as pd 
import tempfile
from run_predictions import get_predictions, make_dictn, get_consensus
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# current file directory
root = os.path.dirname(os.path.abspath(__file__))

def combine_consensus_files(outputs, combined_file):
    # Get a list of all files in the output folder
    files = os.listdir(output_folder)

    # Filter files to only include those ending with '-consensus.csv'
    consensus_files = [f for f in files if f.endswith('-consensus.csv')]

    # Initialize a DataFrame to store the combined data
    combined_df = None

    for file in consensus_files:
        file_path = os.path.join(output_folder, file)
        # Read each output file and extract the 'Consensus' column
        df = pd.read_csv(file_path)
        consensus_column = df['Consensus']

        # Extract the filename without the '-consensus.csv' suffix
        filename = file.replace('-consensus.csv', '')

        # Add the 'Consensus' column to the combined DataFrame with the filename as the column header
        if combined_df is None:
            combined_df = pd.DataFrame({'Consensus': consensus_column})
        else:
            combined_df[filename] = consensus_column

    # Save the combined DataFrame to a new CSV file
    combined_df.to_csv(combined_file, index=False)


# my model
def my_model(smiles_list):
    temp_results_folder = tempfile.mkdtemp()
    try:
        # Call the existing script's function to generate predictions
        get_predictions(temp_results_folder, csv_file)

        # Call the new function to combine the consensus columns
        combined_file = 'consensus_files.csv'
        combine_consensus_files(temp_results_folder, combined_file)
        return combined_file

    finally:
        # Cleanup: Remove the temporary directory and its contents
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

    consensus_df = pd.read_csv(consensus_file)
    # Extract the values from the "Consensus" column
    outputs = consensus_df['Consensus'].tolist()

    for o in outputs:
        writer.writerow([o])
