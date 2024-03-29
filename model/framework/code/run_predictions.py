######################################################################################
# Author: Govinda KC, Computational Science, UTEP                                    #
# Code developed in Sirimulla Research Group (http://sirimullaresearchgroup.com/)    #
# University of Texas at El Paso, Tx, USA                                            #
# Last modified: 01/01/2021                                                          #
######################################################################################

import tempfile, shutil
import glob, os, pandas as pd, numpy as np
import argparse, pickle, json
from get_features import FeaturesGeneration
from smile_standardization import StandardSmiles

#paths
model_path = os.path.dirname(os.path.abspath(__file__))

# SMILES standardization
def standardize(temp_dir, csv_file):
    df = pd.read_csv(csv_file)
    smiles = df['SMILES']
    smiles_standard = []
    for i in range(len(smiles)):
        sd = StandardSmiles()
        stand_smi = sd.preprocess_smi(smiles[i])
        smiles_standard.append(stand_smi)
    df.insert(1, 'SMILES_stand', smiles_standard)
    df['SMILES_stand'].replace('', np.nan, inplace=True)
    df.dropna(subset=['SMILES_stand'], inplace=True)
    df.reset_index(drop=True, inplace=True)
    df.to_csv(temp_dir+'/'+'temp_file_stand.csv', index=False)
    return df
# Features generation
def automate(temp_dir, _file):
    stand_df = standardize(temp_dir, _file)
    features_dictn = dict()
    dictn = json.load(open(os.path.join(model_path, "dictn_models_fp.json"), 'r'))
    fg = FeaturesGeneration()
    pharmacophore = fg.get_fingerprints(stand_df, 'dummy_name', 'tpatf', 'dummy_split', 'dummpy_numpy_folder')
    
    for k, v in dictn.items():
        local_dictn = dict()
        features = fg.get_fingerprints(stand_df, k, v, 'dummy_split', 'dummpy_numpy_folder')
        features_rdkit = fg.get_fingerprints(stand_df, k, 'rdkDes', 'dummy_split', 'dummpy_numpy_folder')
        local_dictn['fingerprint']=features
        local_dictn['rdkDescriptor']=features_rdkit
        local_dictn['pharmacophore']=pharmacophore
        features_dictn[k] = local_dictn
    
    return features_dictn

# Get fingerprint type of best fingerprint models
def make_dictn():
    file_path = os.path.join(model_path,"..","..", "checkpoints/models_tuned_best/*.pkl")
    files = glob.glob(file_path)
    fingerprints, models = [], []
    for f in files:
        filename = os.path.splitext(os.path.basename(f))[0]
        filename2 = filename.split('-')[0]
        if filename2 == 'rdkDes' or filename2 == 'tpatf' or filename2 == 'volsurf':
            continue
        fingerprints.append(filename2)
        model_name = "-".join(filename.split('-')[1:])[:-31]
        models.append(model_name)
    dictn = dict(zip(models, fingerprints))
    return dictn

# Get consensus label
def get_consensus(df):
    consensus_label = []
    for i in range(len(df)):
        sum_proba1 = df['fingerprint'][i] + df['pharmacophore'][i] + df['rdkDescriptor'][i]
        consensus_proba1 = sum_proba1/3
        consensus_label.append(consensus_proba1)
    df.insert(len(df.columns), 'Consensus', consensus_label)
    return df

# Prediction 
def get_predictions(temp_dir, results, csv_file):
    dictn_fp = make_dictn()
    filename = os.path.splitext(os.path.basename(csv_file))[0]
    features_dictn = automate(temp_dir, csv_file)
    model_types = ['3CL', 'ACE2', 'AlphaLISA', 'CoV1-PPE', 'CoV1-PPE_cs', 'CPE', 'cytotox', 'hCYTOX', 'MERS-PPE', 'MERS-PPE_cs', 'TruHit'] 
    
    for m in model_types:
        fptype = dictn_fp[m]
        exact_fpnames = [fptype, 'tpatf', 'rdkDes']
        fpnames = ['fingerprint', 'pharmacophore', 'rdkDescriptor']
        predictions = []
        
        df = pd.read_csv(temp_dir+'/'+'temp_file_stand.csv')

        for fp, fp_name in zip(fpnames, exact_fpnames):
            data = features_dictn[m][fp]
            X_true= data
            models_tuned_dir = os.path.join(model_path,"..","..", "checkpoints/models_tuned_best/")
            pickle_filename = '{}-{}-balanced_randomsplit7_70_15_15.pkl'.format(fp_name, m)
            pickle_path = os.path.join(models_tuned_dir, pickle_filename)
            
            model = pickle.load(open(pickle_path, 'rb'))
            y_pred = model.predict_proba(X_true)[:,1]
            predictions.append(y_pred)
        
        df.insert(len(df.columns), 'fingerprint', predictions[0])
        df.insert(len(df.columns), 'pharmacophore', predictions[1])
        df.insert(len(df.columns), 'rdkDescriptor', predictions[2])
        
        df_consensus = get_consensus(df)
        df_consensus.to_csv(results+'/'+m+'-'+filename+'-consensus.csv', index=False)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="SARS-CoV-2 Activities Prediction")

    parser.add_argument('--csvfile', action='store', dest='csvfile', required=True, \
                        help='csv file with a SMILES column')
    

    parser.add_argument('--results', action='store', dest='results', required=True,\
                        help='Results folder path')    

    args = vars(parser.parse_args())

    csv_file = args['csvfile']
    results = args['results']
    temp_dir = tempfile.mkdtemp() 
    get_predictions(temp_dir, results, csv_file)
    shutil.rmtree(temp_dir)
    print('Done')
