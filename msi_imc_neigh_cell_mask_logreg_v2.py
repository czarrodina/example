#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, balanced_accuracy_score, ConfusionMatrixDisplay
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
import os
import csv
import gc
import multiprocessing as mp
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
import time
from scipy.stats import ttest_ind
import warnings
np.random.seed(9999)
warnings.filterwarnings("ignore")


# In[ ]:


# runs logistic regression on specified Y matrixes for MSI, IMC, and neighorhood features
# due to lack of effect that parameter had during grid search, this run uses a C of 1 and L1 of 0
# as these parameters have help to minimize runtimes
# expected steady state memory usage is 250GB, max during file loading can be above 300GB
# initial data loading iteratively loads, merges, and deletes extra data to help minimize the
# max and longterm memory usage


# In[ ]:


# name of run
fname = 'mask_logreg_original_msi_imc_neigh'
# Y value for classification
yvarn = 'cell_gleason_class'
# classes to exclude
excval = 'Mask2'
# class to use for creating k-folds
splitvar = 'MainGleason'


# In[ ]:


# files containing msi, imc, metadata, and neighborhoods
allf = "../tma_mask_msi_imc_groups_coords_ind_cat_log_v4.csv"
msin = "../tma_mask_prostate_neighborhood_msi_80um_v2.csv"
imcn = "../tma_mask_prostate_neighborhood_imc_80um_cat_log_v4.csv"
msin2 = "../tma_mask_prostate_neighborhood_msi_160um_v2.csv"
imcn2 = "../tma_mask_prostate_neighborhood_imc_160um_cat_log_v4.csv"
msin3 = "../tma_mask_prostate_neighborhood_msi_320um_v2.csv"
imcn3 = "../tma_mask_prostate_neighborhood_imc_320um_cat_log_v4.csv"


# In[ ]:


# reads msi, imc, and 80um neighborhood data
alldf = pd.read_csv(allf)
msidf = pd.read_csv(msin)
imcdf = pd.read_csv(imcn)


# In[ ]:


# drops redundant metadata across neighborhood files (used for initial QC)
msidf = msidf.drop(['Object','core'],axis =1)
imcdf = imcdf.drop(['Object','core'],axis =1)


# In[ ]:


# merges neighborhood then deletes initial dataframes (helps save memory)
nmerge = pd.merge(msidf, imcdf, left_index=True, right_index=True, how='inner')
del msidf
del imcdf


# In[ ]:


# repeats process for 160um neighborhoods
msidf2 = pd.read_csv(msin2)
imcdf2 = pd.read_csv(imcn2)
msidf2 = msidf2.drop(['Object','core'],axis =1)
imcdf2 = imcdf2.drop(['Object','core'],axis =1)


# In[ ]:


nmerge2 = pd.merge(msidf2, imcdf2, left_index=True, right_index=True, how='inner')
del msidf2
del imcdf2


# In[ ]:


# merges 80um and 160um neighborhood features, then deletes extra dataframe
nmerge = pd.merge(nmerge, nmerge2, left_index=True, right_index=True, how='inner')
del nmerge2


# In[ ]:


# repeats the process for 320um neighborhoods
msidf3 = pd.read_csv(msin3)
imcdf3 = pd.read_csv(imcn3)
msidf3 = msidf3.drop(['Object','core'],axis =1)
imcdf3 = imcdf3.drop(['Object','core'],axis =1)


# In[ ]:


nmerge3 = pd.merge(msidf3, imcdf3, left_index=True, right_index=True, how='inner')
del msidf3
del imcdf3


# In[ ]:


nmerge = pd.merge(nmerge, nmerge3, left_index=True, right_index=True, how='inner')
del nmerge3


# In[ ]:


# merges all the neighborhoods with msi/imc and metadata
dfmsi = pd.merge(alldf, nmerge, left_index=True, right_index=True, how='inner')


# In[ ]:


del alldf
del nmerge


# In[ ]:


# drops some of unused metadata
other = ['axis_major_length','axis_minor_length','y','x','area','eccentricity']
dfmsi.drop(other, axis=1, inplace=True)


# In[ ]:


print(dfmsi.shape)


# In[ ]:


# creates a data frame to modify classes, extra dataframe helps save on memory and processor time
dfshort = pd.DataFrame()
dfshort['indication'] = dfmsi['indication']


# In[ ]:


# generates score class
score = []
for index, row in dfshort.iterrows():
    if row['indication'] == 'Benign':
        score.append(0)

    else:
        vals = row['indication'].split('+')
        score.append(int(vals[0]) + int(vals[1]))
        


# In[ ]:


dfshort['totscore'] = score
score = []


# In[ ]:


# generates group class
group = []
for index, row in dfshort.iterrows():
    if row['indication'] == 'Benign':
        group.append(0)
    elif row['indication'] == '3+4':
        group.append(2)
    elif row['indication'] == '4+3':
        group.append(3)
    elif row['totscore'] == 8:
        group.append(4)
    elif row['totscore'] > 8:
        group.append(5)
    elif row['totscore'] < 7:
        group.append(1)
    else:
        group.append(10)


# In[ ]:


dfshort['groupscore'] = group
group = []


# In[ ]:


# adds the group and score classes to the main dataframe, then deletes extra dataframe
dfmsi ['groupscore'] = dfshort['groupscore']
dfmsi['totscore'] = dfshort['totscore']


# In[ ]:


del dfshort


# In[ ]:


# QC for extra class generation
dfmsi['groupscore'].value_counts()


# In[ ]:


# removes excluded class
dfmsi = dfmsi[dfmsi[yvarn] != excval]
dfmsi.reset_index(drop=True, inplace = True)


# In[ ]:


# identifies all features for X matrix
mcols = dfmsi.columns[dfmsi.columns.str.contains('m_')]
pcols= dfmsi.columns[dfmsi.columns.str.contains('p_')]
gcols= dfmsi.columns[dfmsi.columns.str.contains('g_')]
mncols80 = dfmsi.columns[dfmsi.columns.str.contains('mint80_')]
pncols80= dfmsi.columns[dfmsi.columns.str.contains('pint80_')]
gncols80= dfmsi.columns[dfmsi.columns.str.contains('gint80_')]
mncols160 = dfmsi.columns[dfmsi.columns.str.contains('mint160_')]
pncols160= dfmsi.columns[dfmsi.columns.str.contains('pint160_')]
gncols160= dfmsi.columns[dfmsi.columns.str.contains('gint160_')]
mncols320 = dfmsi.columns[dfmsi.columns.str.contains('mint320_')]
pncols320= dfmsi.columns[dfmsi.columns.str.contains('pint320_')]
gncols320= dfmsi.columns[dfmsi.columns.str.contains('gint320_')]
imc = dfmsi.columns[dfmsi.columns.str.contains('imc_')]
imc80 = dfmsi.columns[dfmsi.columns.str.contains('imc80_')]
imc160 = dfmsi.columns[dfmsi.columns.str.contains('imc160_')]
imc320 = dfmsi.columns[dfmsi.columns.str.contains('imc320_')]

allcols = [*mcols, *pcols, *gcols, *mncols80, *pncols80, *gncols80, *mncols160, *pncols160, *gncols160,
          *mncols320, *pncols320, *gncols320,*imc, *imc80, *imc160, *imc320]


# In[ ]:


# QC to make sure contains all desired features
print(allcols)


# In[ ]:


# creates X matrix, class matrix, clinician annotation list, and cell list
allx = dfmsi[allcols]
ally = dfmsi[yvarn]
ayfull = dfmsi['group']
objfull = dfmsi['core'] + '_' + dfmsi['Object'].astype(str)

# scales X feature matrix, helps improve runtimes
X_scaler = StandardScaler()
X_scaler.fit(allx)
allx = X_scaler.transform(allx)


# In[ ]:


# creates metadata for splitting cores into folds
# all cells in the core got into the same split
test_train_metadata = pd.DataFrame()
test_train_metadata['core'] = dfmsi['core']
test_train_metadata['score'] = dfmsi[splitvar]


# In[ ]:


# creates list of cores and their class of interest
# test/train splits attempts to balance main class of the core
setv = []

for core in test_train_metadata['core'].unique():
    temp = test_train_metadata[test_train_metadata['core'] == core]
    mc = temp['score'].value_counts().index[0]
    setv.append([core,mc])

tma_gleason = pd.DataFrame(setv, columns = ['TMAcore','score'])


# In[ ]:


# deletes original data to same memory
del dfmsi
gc.collect()


# In[ ]:


# 3 cross folds
cv = 3
skf = StratifiedKFold(n_splits = cv, shuffle = True)


# In[ ]:


# Generate a train,test iterable
train_test_index = []
for i, (train_index, test_index) in enumerate(skf.split(tma_gleason, tma_gleason["score"])):
    #print(f"{i+1}-th fold, {train_index=}")
    TMAcore_train_index = []
    for j in train_index:
        TMAcore_name = tma_gleason.iloc[j]["TMAcore"]
        # print(TMAcore_name, tma_gleason["Dominant_Gleason_Class"][j])
        TMAcore_train_index.append(test_train_metadata[test_train_metadata["core"] == TMAcore_name].index.to_list())
    TMAcore_train_index = [x for xs in TMAcore_train_index for x in xs]
    #print(f"{i+1}-th fold, {test_index=}")
    TMAcore_test_index = []
    for k in test_index:
        TMAcore_name = tma_gleason.iloc[k]["TMAcore"]
        # print(TMAcore_name, tma_gleason["Dominant_Gleason_Class"][k])
        TMAcore_test_index.append(test_train_metadata[test_train_metadata["core"] == TMAcore_name].index.to_list())
    TMAcore_test_index = [x for xs in TMAcore_test_index for x in xs]
    train_test_index.append((TMAcore_train_index, TMAcore_test_index))


# In[ ]:


# grid search code, not used in this version
'''
logENet = LogisticRegression(penalty = "l2", max_iter = 1000, warm_start = True, solver = "saga", n_jobs = -1)
param_grid = [
  {"C": np.linspace(0.05,1.5,num = 3)}
 ]
clf = GridSearchCV(logENet, param_grid, cv = train_test_index, scoring = "balanced_accuracy", n_jobs = 3, verbose = 3)
clf.fit(allx,ally)
'''


# In[ ]:


'''
# Show the performance statistics of various models
clf_results = pd.DataFrame(clf.cv_results_)
clf_results.sort_values(by = "rank_test_score", inplace = True)
clf_results.to_csv("tma_lr_grid_msi_imc_80_320_neigh_score.csv", index = False)
clf_best_estimator = clf.best_estimator_
'''


# In[ ]:





# In[ ]:


# iterates through each fold, running logistic regression using best parameters
# in this case its C 1 and L1 of 0 due to lack of effect on accuracy
# the run with the highest accuracy is saved for a feature analysis
Y_pred_collection = []
Y_collection = []
objcol = []
yfc = []
maxacc = 0
accl = []
for i in range(cv):
    #LogReg_forCV = LogisticRegression(penalty = "l2", max_iter = 1000, warm_start = False, solver = "saga",C = clf_best_estimator.get_params().get("C"), n_jobs = -1, class_weight = 'balanced')
    LogReg_forCV = LogisticRegression(penalty = "l2", max_iter = 1000, warm_start = False, solver = "saga",C = 1, n_jobs = -1, class_weight = 'balanced')

    train_index, test_index = train_test_index[i]
    # get test_train data
    X_test = allx[test_index]
    X_train = allx[train_index]
    Y_test = ally[test_index]
    Y_train = ally[train_index]
    LogReg_forCV.fit(X_train, Y_train)
    Y_pred = LogReg_forCV.predict(X_test)
    Y_pred_collection.append(Y_pred)
    Y_collection.append(Y_test)
    Y_full = ayfull[test_index]
    objsub = objfull[test_index]
    yfc.append(Y_full)
    objcol.append(objsub)
    accuracy = balanced_accuracy_score(Y_test, Y_pred)
    if accuracy > maxacc:
        clf_best_estimator = LogReg_forCV
        maxacc = accuracy
    accl.append(accuracy)
    print(maxacc)
Y_pred_collection = np.concatenate(Y_pred_collection)
Y_collection =np.concatenate(Y_collection)
yfc = np.concatenate(yfc)
objcol = np.concatenate(objcol)


# In[ ]:


# creates a confusion matrix to display class accuracy
disp = ConfusionMatrixDisplay(confusion_matrix(Y_collection, Y_pred_collection, normalize = "true"), display_labels=clf_best_estimator.classes_)
disp.plot()
plt.savefig(fname + '_Confusion_Matrix.png', bbox_inches='tight', dpi = 300)
plt.show()
plt.close()


# In[ ]:


dfacc = pd.DataFrame(accl)
dfacc.to_csv(fname + "_Generalized_all_features.csv", index = False)


# In[ ]:


# From the best estimator and a list of features, create a dataframe of processed data
# saves feature importance per class and the sum of absolute value across classes
Clusters_Coef_collection = []
column_names = ["Predictors"]
Clusters_Coef_collection.append(allcols)
for i in range(clf_best_estimator.coef_.shape[0]):
    Clusters_Coef_collection.append(clf_best_estimator.coef_[i])
    column_names.append(clf_best_estimator.classes_[i])
for i in range(clf_best_estimator.coef_.shape[0]):
    Clusters_Coef_collection.append(abs(clf_best_estimator.coef_[i]))
    column_names.append(f"abs({clf_best_estimator.classes_[i]})")
Clusters_Coef_collection.append(np.sum(abs(clf_best_estimator.coef_), axis = 0))
column_names.append("sum(abs(Coef))")


# In[ ]:


Clusters_Coef = pd.DataFrame(Clusters_Coef_collection).transpose()
Clusters_Coef.columns = column_names
Clusters_Coef.to_csv(os.path.join(fname + "_Generalized_Coef.csv"), index = False)
Clusters_Coef.head()


# In[ ]:





# In[ ]:


# saves class and prediction information for each cell
allydf = pd.DataFrame()


# In[ ]:


allydf['Gleason Score'] = Y_collection
allydf['Full Gleason'] = yfc
allydf['Prediction'] = Y_pred_collection
allydf['Cell'] = objcol


# In[ ]:


allydf.to_csv(os.path.join(fname + "_ally_prediction.csv"), index = False)


# In[ ]:


np.sum(accl)/3


# In[ ]:


# saves accuracy information for each run
dfa = pd.DataFrame()
dfa['fold'] = ['1','2','3']
dfa['accuracy'] = accl
dfa.to_csv(os.path.join(fname + "_fold_accuracy.csv"), index = False)

