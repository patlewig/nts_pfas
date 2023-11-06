import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import PandasTools
from rdkit.Chem import AllChem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
import os
import sys
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import fcluster
from scipy.cluster.hierarchy import cophenet
from scipy.cluster.hierarchy import inconsistent





def mk_fp(df):

    '''
    Function to create a Morgan FP  of length 1024 and radius 3. Input file is expecting a dtxsid and smiles column in a df. 
    Expects dtxsid as identifier and smiles as SMILES.
    Returns df of index as dtxsid values and columns as 1024 morgan FP
    '''

    MOLS = dict(zip(df['dtxsid'], df['smiles']))
    MOLS = {k:Chem.MolFromSmiles(v) for k,v in MOLS.items()}
    MOLS = {i:j for i,j in MOLS.items() if j}
    FP0 = pd.DataFrame([np.array(AllChem.GetMorganFingerprintAsBitVect(i,3,1024)) for i in MOLS.values()])
    FP0.index = MOLS.keys()
    FP0.columns = ['mrgn_%d'%i for i in FP0.columns]
    return FP0

def distance_matrix(df):
    '''
    Function to create a pairwise square distance matrix using the Jaccard index
    '''
    D_mgrn = pd.DataFrame(squareform(pdist(df, 'jaccard')), columns = df.index, index = df.index)

    return D_mgrn

def medoid_calc(D):
    '''
    Calculate the medoid of a dataframe
    '''

    a = np.argmin(D.sum(axis = 1))
    return D.index[a]


def nn_neighbours(chem, df, n = 10):
    '''
    Calculate nearest neighbours for a chemical from an existing distance matrix
    '''
    ids = {i:x for i,x in enumerate(df.index)}
    index_chem = df[chem].values
    temp_chem = np.argpartition(index_chem, n)
    chem_nn = temp_chem[:n]
    mydf = pd.DataFrame(zip([ids[x] for x in chem_nn], index_chem[chem_nn]), columns = ['dtxsid', 'dist'])
    return mydf

def ecdf(data):
    '''
    Calculate the ecdf of a 1d-array
    Returns the x, y needed to plot the empirical ecdf
    '''
    n = len(data)
    x = np.sort(data)
    y = np.arange(1, n+1)/n
    return x, y

def plug_in(D, ptiles = 95, label = None):
    '''
    computes the dictionary of cluster information and summary metrics 
    '''
    summary = {}
    Zm = linkage(squareform(D), 'ward')
    summary['C'], summary['Coph_dists'] = cophenet(Zm, squareform(D))
    summary['5PCTILE'] = np.percentile(squareform(D), 5)
    summary['10PCTILE'] = np.percentile(squareform(D), 10)
    summary['Z'] = Zm
    summary['label'] = label
    summary['std'] = squareform(D).std()
    summary['50PCTILE'] = np.percentile(squareform(D), 50)
    summary['max'] = squareform(D).max()
    summary['ecdfx'], summary['ecdfy'] = ecdf(squareform(D))
    summary['diff_dist'] = squareform(D).max() - squareform(D).min()
    
    return summary
    
    
def prep(df):
    '''
    Preparatory wrangling for between cluster distances
    '''
    df = df.where(np.triu(np.ones(df.shape)).astype(np.bool))
    df = df.stack().reset_index()
    df = df[df['level_0'] != df['level_1']]
    return df[0].values
    
def plug_in_stack(df, label = None):
    summary = {}
    summary['label'] = label
    summary['dist'] = prep(df)
    summary['ecdfx'], summary['ecdfy'] = ecdf(summary['dist'])
    summary['diff_dist'] = summary['dist'].max() - summary['dist'].min()
    summary['5PCTILE'] = np.percentile(summary['dist'], 5)
    summary['10PCTILE'] = np.percentile(summary['dist'], 10)
    summary['50PCTILE'] = np.percentile(summary['dist'], 50)
    summary['std'] = summary['dist'].std()
    summary['max'] = summary['dist'].max()
    return summary
    
 

