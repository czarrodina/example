#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from pathlib import Path
import pandas as pd
import os
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist, squareform

import sys


# In[ ]:


# generates neighborhood m/z intensities for each cell in each core
# by utilizing a guassian kernal that averages the intensity of the neighbors
# based upon the distance from middle of the center cell to the middle
# of the neighboring cells

# input files containing MSI/IMC features
msifile = "/mnt/d/tma_mask_msi_imc_groups_coords_ind_v2.csv"
dfmsi = pd.read_csv(msifile)


# In[ ]:


# distance for creating the neighborhoods
cutoff = 320


# In[ ]:


# identifies the MSI features
mcols = dfmsi.columns[dfmsi.columns.str.contains('m_')]
pcols= dfmsi.columns[dfmsi.columns.str.contains('p_')]
gcols= dfmsi.columns[dfmsi.columns.str.contains('g_')]
msicols = [*mcols, *pcols, *gcols]


# In[ ]:


# identifies the IMC columns
imc = ['imc_CXCR3','imc_COLLAGEN','imc_P504S','imc_34BE12','imc_CD146',
       'imc_ERG','imc_PSMA','imc_CD11c','imc_P63','imc_CD16',	
        'imc_CD163','imc_CD45RA','imc_CD15','imc_PTEN','imc_CD56',
       'imc_CD45','imc_CD44','imc_CD14','imc_FoxP3','imc_CD4','imc_CD11b',
        'imc_CD74','imc_CD68','imc_CCR6','imc_CD20','imc_CD8',
       'imc_MBD-1','imc_C-MYC','imc_PD-1','imc_NKX3.1','imc_EFCAB4B',
        'imc_Ki67','imc_CD206','imc_CD3','imc_HLA-DR','imc_DCR3',
       'imc_CD45RO','imc_CD7','imc_CD49A','imc_CD10','imc_SMA',
        'imc_DNA1','imc_DNA2','imc_E-Cadherin','imc_ATPaseCalculin','imc_VCL',
       'imc_PanKeratin','imc_CD57']

# shrinks the dataset size to make it more RAM friendly
other = ['axis_major_length','axis_minor_length', 'area','eccentricity']
dfmsi.drop(imc, axis=1, inplace=True)
dfmsi.drop(other, axis=1, inplace=True)


# In[ ]:


# stores all of neighborhood intensities
allint = []
numc = 0

# for each unique core in the TMA
for core in dfmsi['core'].unique():
    # subselect the core features
    coremsi = dfmsi[dfmsi['core']==core]
    allx = dfmsi[msicols]

    # gets X/Y coordinates for each cell
    disdf = pd.DataFrame()
    disdf['X'] = coremsi['x']
    disdf['Y'] = coremsi['y']

    # creates of matrix of distances between cells
    dist = pdist(disdf.values, metric='euclidean')
    distm = squareform(dist)

    # sets cells outside the neighborhood to a distance of 0
    distance_df = pd.DataFrame(distm, index=disdf.index, columns=disdf.index)
    distance_df[distance_df > cutoff] = 0

    corecols = []

    # for each cell in the core, create a neighborhood
    # apply a guassian kernal to the neighbor intensities
    # and then sum up the values to get the intensity for the neighborhood
    for index, row in distance_df.iterrows():
        neigh = row[row != 0]
        temp = allx.iloc[neigh.index]
        kernel = np.exp(-neigh/neigh.std())

        kern_df = temp.mul(kernel, axis=0)
        gauss = kern_df.sum(axis=0)

        corecols.append(gauss)

    # stores neighborhood intensities for all the cells in the core
    coredf = pd.DataFrame(corecols, columns = msicols)
    coredf['Object'] = dfmsi['Object']
    coredf['core'] = core

    allint.append(coredf)
    numc+=1
    print (str(numc), end=" ") 


# In[ ]:


# creates the neighborhood dataframe for all cores, renames the features to indicate distance
# and then saves the resulting dataframe
allintdf = pd.concat(allint)


# In[ ]:


allintdf.columns = allintdf.columns.str.replace('m_', 'mint' + str(cutoff) + '_').str.replace('p_', 'pint' + str(cutoff) + '_').str.replace('g_', 'gint' + str(cutoff) + '_')


# In[ ]:


allintdf


# In[ ]:


allintdf.to_csv("/mnt/d/msi_neighborhood_matrixes/tma_mask_prostate_neighborhood_msi_" + str(cutoff) + "um_v2.csv",index=False)


# In[ ]:





# In[ ]:




