import pickle
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import os

class MoleculeProcessor:
    def __init__(self, smiles, dtx):
        self.mol = Chem.MolFromSmiles(smiles)
        self.dtx = dtx
    
    def chain_length(self, ch=30):
        mysr = 'C(F)(F)'
        mylst = []
        for n in range(1, ch):
            a = self.mol.HasSubstructMatch(Chem.MolFromSmarts(''.join(mysr * n)))
            mylst.append(a)
        return mylst.index(False)

    def get_morgan_fingerprint(self, radius = 3, nBits=1024):
        mgrn_df = pd.DataFrame([np.array(AllChem.GetMorganFingerprintAsBitVect(self.mol,radius,nBits))] )
        mgrn_df.columns = ['mrgn_%d'%i for i in mgrn_df.columns]
        mgrn_df.index = [self.dtx]
        return mgrn_df

class CategoryManager:
    def __init__(self, categories):
        self.categories = categories
        
    def make_df(self, smiles, dtx, cat='others'):
        '''
    Construct the appropriate input df to facilitate RF predictions
        '''
        if cat not in self.categories:
            raise ValueError(f"Category {cat} is not recognized.") 
        mol_processor = MoleculeProcessor(smiles, dtx)
        df1 = mol_processor.get_morgan_fingerprint()
        df1['category'] = cat
        df1['chain_length'] = mol_processor.chain_length()
        df1['category'] = df1['category'].astype('category')
        return df1
    
class CategoryPrediction:
    def __init__(self, rev_dict='reverse_dict2.pkl', rf_model='final_model_v2.sav'):

        self.rev_dict = pickle.load(open(rev_dict, 'rb'))
        self.rf_model = pickle.load(open(rf_model, 'rb'))

    def make_prediction(self, df):
        '''
        Make prediction using the RF model and convert back into the original terminal categories
        '''
        pred = self.rf_model.predict(df)[0]
        return self.rev_dict[pred]


categories = ['Aromatic PFASs',
 'HFCs',
 'Other PFASs',
 'Other PFASs, cyclic',
 'PASF-based substances',
 'PFAA precursors',
 'PFAA precursors, cyclic',
 'PFAAs',
 'PFAAs, cyclic',
 'PolyFCA derivatives',
 'Polyfluoroalkanes',
 'Polyfluoroalkyl acids',
 'Polyfluoroalkyl acids, cyclic',
 'Si PFASs',
 'n:2 fluorotelomer-based substances',
 'others',
 'others, cyclic',
 'unclassified']





