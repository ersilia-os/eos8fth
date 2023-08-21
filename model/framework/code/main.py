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
temp_folder = tempfile.mkdtemp()
print(temp_folder)

# read SMILES from .csv file, assuming one column with header
with open(input_file, "r") as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    smiles_list = [r[0] for r in reader]

data_with_header = [['SMILES']] + [[smiles] for smiles in smiles_list]

with open(os.path.join(temp_folder,"tmp_input.csv"), "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerows(data_with_header)
csv_file = os.path.join(temp_folder, "tmp_input.csv")
get_predictions(temp_folder, temp_folder, csv_file)
files = os.listdir(temp_folder)
consensus_files = [f for f in files if f.endswith('-consensus.csv')]

combined_df = pd.DataFrame()

for file in consensus_files:
    file_path = os.path.join(temp_folder, file)
    df = pd.read_csv(file_path)
    consensus_column = df['Consensus']
    exp_name = file.split("-tmp_input-consensus.csv")[0]
    combined_df[exp_name] = consensus_column

combined_df.to_csv(output_file, index=False)
        
#if os.path.exists(temp_folder):
    #shutil.rmtree(temp_folder)
