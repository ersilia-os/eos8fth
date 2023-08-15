import sys,os,glob
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import MACCSkeys, AllChem
from rdkit.Chem.AtomPairs import Pairs, Torsions
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import rdMolDescriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import AllChem, Descriptors
import pandas as pd
import pickle
import tempfile, os
import shutil
from sklearn.preprocessing import LabelEncoder
from config import fpFunc_dict
import argparse
from sklearn.impute import SimpleImputer

descriptor_names =  ['MaxEStateIndex', 'MinEStateIndex', 'MaxAbsEStateIndex', 'MinAbsEStateIndex', 'qed', 'MolWt', 'HeavyAtomMolWt', 'ExactMolWt', 'NumValenceElectrons', 'NumRadicalElectrons', 'MaxPartialCharge', 'MinPartialCharge', 'MaxAbsPartialCharge', 'MinAbsPartialCharge', 'FpDensityMorgan1', 'FpDensityMorgan2', 'FpDensityMorgan3', 'BalabanJ', 'BertzCT', 'Chi0', 'Chi0n', 'Chi0v', 'Chi1', 'Chi1n', 'Chi1v', 'Chi2n', 'Chi2v', 'Chi3n', 'Chi3v', 'Chi4n', 'Chi4v', 'HallKierAlpha', 'Ipc', 'Kappa1', 'Kappa2', 'Kappa3', 'LabuteASA', 'PEOE_VSA1', 'PEOE_VSA10', 'PEOE_VSA11', 'PEOE_VSA12', 'PEOE_VSA13', 'PEOE_VSA14', 'PEOE_VSA2', 'PEOE_VSA3', 'PEOE_VSA4', 'PEOE_VSA5', 'PEOE_VSA6', 'PEOE_VSA7', 'PEOE_VSA8', 'PEOE_VSA9', 'SMR_VSA1', 'SMR_VSA10', 'SMR_VSA2', 'SMR_VSA3', 'SMR_VSA4', 'SMR_VSA5', 'SMR_VSA6', 'SMR_VSA7', 'SMR_VSA8', 'SMR_VSA9', 'SlogP_VSA1', 'SlogP_VSA10', 'SlogP_VSA11', 'SlogP_VSA12', 'SlogP_VSA2', 'SlogP_VSA3', 'SlogP_VSA4', 'SlogP_VSA5', 'SlogP_VSA6', 'SlogP_VSA7', 'SlogP_VSA8', 'SlogP_VSA9', 'TPSA', 'EState_VSA1', 'EState_VSA10', 'EState_VSA11', 'EState_VSA2', 'EState_VSA3', 'EState_VSA4', 'EState_VSA5', 'EState_VSA6', 'EState_VSA7', 'EState_VSA8', 'EState_VSA9', 'VSA_EState1', 'VSA_EState10', 'VSA_EState2', 'VSA_EState3', 'VSA_EState4', 'VSA_EState5', 'VSA_EState6', 'VSA_EState7', 'VSA_EState8', 'VSA_EState9', 'FractionCSP3', 'HeavyAtomCount', 'NHOHCount', 'NOCount', 'NumAliphaticCarbocycles', 'NumAliphaticHeterocycles', 'NumAliphaticRings', 'NumAromaticCarbocycles', 'NumAromaticHeterocycles', 'NumAromaticRings', 'NumHAcceptors', 'NumHDonors', 'NumHeteroatoms', 'NumRotatableBonds', 'NumSaturatedCarbocycles', 'NumSaturatedHeterocycles', 'NumSaturatedRings', 'RingCount', 'MolLogP', 'MolMR', 'fr_Al_COO', 'fr_Al_OH', 'fr_Al_OH_noTert', 'fr_ArN', 'fr_Ar_COO', 'fr_Ar_N', 'fr_Ar_NH', 'fr_Ar_OH', 'fr_COO', 'fr_COO2', 'fr_C_O', 'fr_C_O_noCOO', 'fr_C_S', 'fr_HOCCN', 'fr_Imine', 'fr_NH0', 'fr_NH1', 'fr_NH2', 'fr_N_O', 'fr_Ndealkylation1', 'fr_Ndealkylation2', 'fr_Nhpyrrole', 'fr_SH', 'fr_aldehyde', 'fr_alkyl_carbamate', 'fr_alkyl_halide', 'fr_allylic_oxid', 'fr_amide', 'fr_amidine', 'fr_aniline', 'fr_aryl_methyl', 'fr_azide', 'fr_azo', 'fr_barbitur', 'fr_benzene', 'fr_benzodiazepine', 'fr_bicyclic', 'fr_diazo', 'fr_dihydropyridine', 'fr_epoxide', 'fr_ester', 'fr_ether', 'fr_furan', 'fr_guanido', 'fr_halogen', 'fr_hdrzine', 'fr_hdrzone', 'fr_imidazole', 'fr_imide', 'fr_isocyan', 'fr_isothiocyan', 'fr_ketone', 'fr_ketone_Topliss', 'fr_lactam', 'fr_lactone', 'fr_methoxy', 'fr_morpholine', 'fr_nitrile', 'fr_nitro', 'fr_nitro_arom', 'fr_nitro_arom_nonortho', 'fr_nitroso', 'fr_oxazole', 'fr_oxime', 'fr_para_hydroxylation', 'fr_phenol', 'fr_phenol_noOrthoHbond', 'fr_phos_acid', 'fr_phos_ester', 'fr_piperdine', 'fr_piperzine', 'fr_priamide', 'fr_prisulfonamd', 'fr_pyridine', 'fr_quatN', 'fr_sulfide', 'fr_sulfonamd', 'fr_sulfone', 'fr_term_acetylene', 'fr_tetrazole', 'fr_thiazole', 'fr_thiocyan', 'fr_thiophene', 'fr_unbrch_alkane', 'fr_urea']
descriptor_funcs = {desc_name: getattr(Descriptors, desc_name) for desc_name in descriptor_names}

class FeaturesGeneration:
    def __init__(self):
        self.fingerprints = []
    def get_fingerprints(self, df, model, fp_name, split, numpy_folder):

        smiles_list = df['SMILES_stand'].to_list()
        
        not_found = []
        for smi in smiles_list:
            try: 
                m = Chem.MolFromSmiles(smi)
            
                can_smi = Chem.MolToSmiles(m, True)
        
                # Calculate only the specific 200 descriptors
                descriptor_values = []
                for descriptor_name in descriptor_names:
                    descriptor_func = descriptor_funcs.get(descriptor_name)
                    if descriptor_func:
                        descriptor_value = descriptor_func(m)
                        descriptor_values.append(descriptor_value)
                    else:
                        descriptor_values.append(np.nan)
                
                self.fingerprints.append(descriptor_values)
            except:
                # Handle exceptions and add placeholder values
                self.fingerprints.append([np.nan] * len(descriptor_names))

                
                if fp_name == 'tpatf':
                    add = [np.nan for i in range(self.fingerprints[0].shape[1])]
                elif fp_name == 'rdkDes':
                    add = [np.nan for i in range(len(self.fingerprints[0]))]
                else:
                    add = [np.nan for i in range(len(self.fingerprints[0]))]
                tpatf_arr = np.array(add, dtype=np.float32)
                self.fingerprints.append(tpatf_arr) 
                
                pass
            
        self.fingerprints = np.array(self.fingerprints)  

        if fp_name == 'rdkDes':
            X = np.array(self.fingerprints)
            ndf = pd.DataFrame.from_records(X)
            ndf.isnull().sum().sum()
            r, _ = np.where(df.isna())
            ndf.isnull().sum().sum()

            for col in ndf.columns:
                ndf[col].fillna(ndf[col].mean(), inplace=True)
            ndf.isnull().sum().sum()
            X = ndf.iloc[:,0:].values
            fp_array = ( np.asarray((X), dtype=object) )
            X = X.astype(np.float32)
            X = np.nan_to_num(X)
            scalers_dir = os.path.abspath('eos8fth/model/checkpoints/scalers')
            pickle_filename = '{}-rdkDes_scaler.pkl'.format(model)
            pickle_path = os.path.join(scalers_dir, pickle_filename)
            rdkDes_scaler = pickle.load(open(pickle_path, 'rb'))
            X = rdkDes_scaler.transform(X)

        else:
            fp_array = ( np.asarray((self.fingerprints), dtype=object) )
            X = np.vstack(fp_array).astype(np.float32)
            imp_median = SimpleImputer(missing_values=np.nan, strategy='median')
            imp_median.fit(X)  
            X = imp_median.transform(X)
        
#        Y = df['Label'].values
#        Y = Y.reshape(Y.shape[0],1)
#        Y = np.vstack(Y).astype(np.float32)
        final_array = X #np.concatenate((X, Y), axis=1)
        self.fingerprints = []
        return final_array