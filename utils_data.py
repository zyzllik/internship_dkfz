from sklearn import mixture

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from scipy.sparse import csr_matrix
import torch

from tqdm import tqdm

import anndata as ad
import gc


# Select indexes randomly (for train-val-test split)
def select_idx(idx_initial, percentage):

    # Number of data points to be selected
    n_selected = np.round(idx_initial.shape[0]*percentage).astype(int)
    # Select percentage% indexes randomly
    rng = np.random.default_rng()
    chosen_idx = rng.choice(idx_initial, n_selected, replace=False)
    # Not selected indexes
    excluded_idx = np.setdiff1d(idx_initial, chosen_idx)
    
    return chosen_idx, excluded_idx

# Split for one genotype
def split_sample(onehot_genotype, p_train, p_val):
    idx_all = np.where(onehot_genotype==1)[0]
    train_idx, other_idx = select_idx(idx_all, p_train)
    val_idx, test_idx = select_idx(other_idx, p_val/(1-p_train))
    return train_idx, val_idx, test_idx

# Split all data
def dataset_split(hom_alt, hom_ref, het, expression, lib_size, p_train=0.7, p_val=0.2): # 70% 20% 10% split
    
    # Split data for each genotype separatly
    train_hom_alt, val_hom_alt, test_hom_alt = split_sample(hom_alt, p_train, p_val)
    train_hom_ref, val_hom_ref, test_hom_ref = split_sample(hom_ref, p_train, p_val)
    train_het, val_het, test_het = split_sample(het, p_train, p_val)
    
    # Get union of indexes in train/validation/test over all genotypes
    train_idx = np.concatenate([train_hom_alt, train_hom_ref, train_het])
    val_idx = np.concatenate([val_hom_alt, val_hom_ref, val_het])
    test_idx = np.concatenate([test_hom_alt, test_hom_ref, test_het])
    
    # Create datasets
    train_data = [torch.tensor(hom_alt[train_idx]), 
                  torch.tensor(hom_ref[train_idx]), 
                  torch.tensor(het[train_idx]), 
                  torch.tensor(lib_size[train_idx]), 
                  torch.tensor(expression[train_idx])]
    val_data = [torch.tensor(hom_alt[val_idx]), 
                torch.tensor(hom_ref[val_idx]), 
                torch.tensor(het[val_idx]), 
                torch.tensor(lib_size[val_idx]), 
                torch.tensor(expression[val_idx])]
    test_data = [torch.tensor(hom_alt[test_idx]), 
                 torch.tensor(hom_ref[test_idx]), 
                 torch.tensor(het[test_idx]), 
                 torch.tensor(lib_size[test_idx]), 
                 torch.tensor(expression[test_idx])]
    
    return train_data, val_data, test_data

# Get data for one eQTL
def get_one_eqtl(selected_eqtl, data, lib_size, raw = False):
    data_gene = data[data.obs.old_cell_label==selected_eqtl['Cell type'], selected_eqtl['Gene Ensembl ID']]
    df_gene = data_gene.to_df()
    lib_size = lib_size[data.obs.old_cell_label==selected_eqtl['Cell type']]

    samples = ['hom_ref_samples', 'het_samples', 'hom_alt_samples']
    for sample in samples:
        search_str = selected_eqtl[sample]
        search_str = search_str.replace(',', '$|^') # Transform to regex
        search_str = '^'+search_str + '$'
        mask = data_gene.obs.individual.str.contains(search_str, regex=True)
        df_gene[sample] = mask.astype(int)

    expression = np.array(df_gene[selected_eqtl['Gene Ensembl ID']])
    hom_alt = np.array(df_gene['hom_alt_samples'])
    hom_ref = np.array(df_gene['hom_ref_samples'])
    het = np.array(df_gene['het_samples'])

    train_data, val_data, test_data = dataset_split(hom_alt, hom_ref, het, expression, lib_size, p_train=0.7, p_val=0.2)

    train_x = train_data[:4]
    train_y = train_data[4]
    val_x = val_data[:4]
    val_y = val_data[4]
    test_x = test_data[:4]
    test_y = test_data[4]
    
    return train_x, train_y, val_x, val_y, test_x, test_y

# Get top n_top eQTLs according to stat
def get_top_eqtls(n_top, stat = 'qvalue', raw = False):
    # Load data
    print('Loading expression data..')
    data = ad.read_h5ad('/home/e860a/chernova/my_onek1k_data/onek1k.norm_pflofpf.h5ad')
    print('Expression data loaded')

    if raw:
        print('Getting the raw counts')
        data.X = data.layers['counts']
        lib_size = np.asarray(csr_matrix.sum(data.X, axis=1)).reshape(data.X.shape[0],)
    else:
        lib_size = np.zeros(data.X.shape[0])
        
    eqtl = pd.read_csv('/home/e860a/chernova/my_onek1k_data/eQTLs_patients.csv')

    # Load TF data
    tf_df = pd.read_csv('/home/e860a/chernova/my_onek1k_data/TF_list.csv')
    tf_df = tf_df[tf_df['Is TF?']=='Yes']

    # Just select eQTLs acting on TFs
    tf_ids = np.intersect1d(eqtl['Gene Ensembl ID'], tf_df['Ensembl ID'])
    search_str = '|'.join(tf_ids)
    tf_eqtl = eqtl[eqtl['Gene Ensembl ID'].str.contains(search_str, regex=True)]
    tf_eqtl_sort = tf_eqtl.sort_values(stat)

    selected_eqtls = tf_eqtl_sort.iloc[:n_top, :]
    
    all_model_data = {}
    print('Getting data for single SNPs')
    for eqtl_idx in tqdm(range(n_top)):
        selected = selected_eqtls.iloc[eqtl_idx, :]
        model_data = get_one_eqtl(selected, data, lib_size, raw=raw)
        name = selected['SNP'] + '_' + selected['Cell type'] + '_' + selected['Gene ID']
        all_model_data[name] = model_data
    del(data)
    gc.collect()
    return all_model_data