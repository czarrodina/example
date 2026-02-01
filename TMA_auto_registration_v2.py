#!/usr/bin/env python
# coding: utf-8

# In[1]:


import itk
from itkwidgets import compare, checkerboard
import numpy as np
import SimpleITK as sitk
from tifffile import imread
from tifffile import imwrite
import matplotlib.pyplot as plt
import SimpleITK as sitk
import xml.etree.ElementTree as ET
from tifffile import TiffFile
from ome_types import from_xml
import geojson
from shapely.geometry import Point, Polygon
from skimage.draw import polygon
import pandas as pd
import PIL

import os
import shutil
import glob


# In[2]:


# registration code is split into two parts, the first is code for all the functions, and the second runs those function 


# In[3]:


# Note: there are two optical image registrations, one inverts the image and one does not
# The example below does not. It all depends on if the image created
# is side scanning or not (side scanning will need to be inverted and the background removed as 
# the the non-tissue will be reversed and the tissue will be color shifted)


# In[4]:


# Note2: this was specifically designed for the structure of the current run, as file management was the bigger issue than 
# making the methods work for any set of TMAs  (future MSI projects should include my recommendations
# for data aquisition, microscopy, and data storage, to help streamline future code)
# each cores can have 50+ files at the end, with optional masks and temporary masked images
# generating these files requires a lot of information, time, and knowledge


# In[5]:


# location of the imc and H&E files
dirimc = r"D:\tma_registrations\tmas_masks_he\PR1211_181"
# location of the msi files
dirmsi = r"D:\tma_registrations\tma_10_4_24_groups\PR1121"


# In[6]:


# creates an overview plot for the results of the registration
def plot_reg(result, moving, fixed, fname):
    plt.ioff()
    fig, axs = plt.subplots(1,3, sharey=False, figsize=[30,30])
    plt.figsize=[100,100]
    axs[0].imshow(result)
    axs[0].set_title('Result', fontsize=30)
    axs[1].imshow(moving)
    axs[1].set_title('Moving', fontsize=30)
    axs[2].imshow(fixed)
    axs[2].set_title('Fixed', fontsize=30)
    plt.savefig(fname, bbox_inches = 'tight')
    plt.close()


# In[7]:


# creates a list for moving around all the necessary images for the MSI - IMC registration
# I found this much easier than to have dozens of parameters in each method
# It does obfuscate the data though
def msi_imc_files(mflist):
    num = 0
    flist = {}
    for file in mflist:
        if len(file) > 0:
            flist[mlorder[num]] = file[0]
        elif 'mask' in mlorder[num]:
            flist[mlorder[num]] = ''
        else:
            return {}
        num+=1
    return flist 


# In[8]:


# creates a list for moving around all the necessary images for the H&E - IMC registration
def he_imc_files(filelist):
    num = 0
    flist = []
    for file in filelist:
        if len(file) > 0:
            flist.append(file[0])
        elif num ==3:
            flist.append('')
        else:
            return []
        num +=1
    return flist 


# In[9]:


# if the mask exists, then create a temporary image that only includes the masked areas in the folder where the image file is located
def geomask(maskfile, imagefile):
    imf = imread(imagefile)
    
    with open(maskfile, 'r') as f:
        hm = geojson.load(f)
        hmarr = np.array(hm['features'][0]['geometry']['coordinates'])
        shmarr = np.squeeze(hmarr)
        
        img = np.zeros((imf.shape[1],imf.shape[0]), 'uint8')
        rr, cc = polygon(shmarr[:,0], shmarr[:,1], img.shape)
        img[rr,cc] = 1
        imgt = img.transpose()
       
        imset = []
        for i in range (0,imf.shape[2]):
            imf2 = imf[:,:,i]
            immasked = imf2 * imgt
            imset.append(immasked)

        imarr  = np.array(imset)
        immasked = np.transpose(imarr, (1, 2, 0))
        imwrite(imagefile.replace('.tiff','_temp.tiff'), immasked, metadata ={'PhysicalSizeX':"1", 'PhysicalSizeXUnit':"µm", 'PhysicalSizeY':"1", 'PhysicalSizeYUnit':"µm"}, ome=True,  bigtiff=True)

    return imagefile.replace('.tiff','_temp.tiff')


# In[10]:


# inverts the image and then removes purely white pixels (basically background)
# useful for side scanned images
def invert(image):
    test = PIL.Image.open(image)
    test_n = np.array(test)
    t_inv = PIL.ImageOps.invert(test)
    t_inv_n = np.array(t_inv)
    t_inv_mask = set_to_zero_if_greater(t_inv_n, 254)

    return t_inv_mask


# In[11]:


# helps for easy removal if high intensity pixels (useful for both H&E and side scanned MSI)
def set_to_zero_if_greater(arr, value):
    """Sets elements to 0 if all x, y, z components are greater than the given value."""

    # Create a boolean mask where True indicates elements that meet the condition
    mask = (arr[:,:, 0] > value) & (arr[:,:, 1] > value) & (arr[:,:, 2] > value)

    # Set elements to 0 where the mask is True
    arr[mask] = 0

    return arr


# In[12]:


# registers the optical images to the IMC, this is the version for regular optical images
def register_msi_imc(mdict):
    # list of registration images and files
    glyregf = []
    pepregf = []
    metregf = []

    # loads images into itk format
    fixed_image = itk.imread(mdict['imcfile'], itk.F)
    moving_gly = itk.imread(mdict['mgopt'], itk.F)
    moving_pep = itk.imread(mdict['mpopt'], itk.F)
    moving_met = itk.imread(mdict['mmopt'], itk.F)

    # if the images have a mask, apply the mask
    if mdict['imcmask'] != '':
        imcu = geomask(mdict['imcmask'], mdict['imcfile'])
        fixed_image = itk.imread(imcu, itk.F)

    if mdict['mglymask'] != '':
        mglyu = geomask(mdict['mglymask'], mdict['mgopt'])
        moving_gly = itk.imread(mglyu, itk.F)

    if mdict['mpepmask'] != '':
        mpepu = geomask(mdict['mpepmask'], mdict['mpopt'])
        moving_pep = itk.imread(mpepu, itk.F)

    if mdict['mmetmask'] != '':
        mmetu = geomask(mdict['mmetmask'], mdict['mmopt'])
        moving_met = itk.imread(mmetu, itk.F)

    # loads the landmarks files and generates the file needed for itk in the core directory
    imcpland = get_landmarks(mdict['imcpland'])
    imcgland = get_landmarks(mdict['imcgland'])
    imcmland = get_landmarks(mdict['imcmland'])

    mgland = get_landmarks(mdict['mgland'])
    mpland = get_landmarks(mdict['mpland'])
    mmland = get_landmarks(mdict['mmland'])

    ufile = open(mdict['imcpland'].replace('.txt','_temp.txt'), 'w')
    ufile.write(imcpland)
    ufile.close()

    ufile = open(mdict['imcgland'].replace('.txt','_temp.txt'), 'w')
    ufile.write(imcgland)
    ufile.close()

    ufile = open(mdict['imcmland'].replace('.txt','_temp.txt'), 'w')
    ufile.write(imcmland)
    ufile.close()

    ufile = open(mdict['mgland'].replace('.txt','_temp.txt'), 'w')
    ufile.write(mgland)
    ufile.close()

    ufile = open(mdict['mpland'].replace('.txt','_temp.txt'), 'w')
    ufile.write(mpland)
    ufile.close()

    ufile = open(mdict['mmland'].replace('.txt','_temp.txt'), 'w')
    ufile.write(mmland)
    ufile.close()
    
    # registers the glycan optical image to the IMC
    result_image1, result_transform_parameters1 = itk.elastix_registration_method(fixed_image, moving_gly,
        fixed_point_set_file_name=mdict['imcgland'].replace('.txt','_temp.txt'),
        moving_point_set_file_name=mdict['mgland'].replace('.txt','_temp.txt'),
        log_to_console=True,
        parameter_object=rab_param)

    glyregf.append(fixed_image)
    glyregf.append(moving_gly)
    glyregf.append(result_image1)
    glyregf.append(result_transform_parameters1)  

    # registers the peptide optical image to the IMC
    result_image1, result_transform_parameters1 = itk.elastix_registration_method(fixed_image, moving_pep,
        fixed_point_set_file_name=mdict['imcpland'].replace('.txt','_temp.txt'),
        moving_point_set_file_name=mdict['mpland'].replace('.txt','_temp.txt'),
        log_to_console=True,
        parameter_object=rab_param)

    pepregf.append(fixed_image)
    pepregf.append(moving_pep)
    pepregf.append(result_image1)
    pepregf.append(result_transform_parameters1)

    # registers the mets optical image to the IMC
    result_image1, result_transform_parameters1 = itk.elastix_registration_method(fixed_image, moving_met,
        fixed_point_set_file_name=mdict['imcmland'].replace('.txt','_temp.txt'),
        moving_point_set_file_name=mdict['mmland'].replace('.txt','_temp.txt'),
        log_to_console=True,
        parameter_object=rab_param)

    metregf.append(fixed_image)
    metregf.append(moving_met)
    metregf.append(result_image1)
    metregf.append(result_transform_parameters1)

    return glyregf, pepregf, metregf


# In[13]:


# registers the optical images to the IMC, this is the version for the side scanning
# as an initial inversion is done on the optical images
def register_msi_imc_invert(mdict):
    # list of registration images and transformation files
    glyregf = []
    pepregf = []
    metregf = []

    # creates inversions of the RBG optical images
    mginv = invert(mdict['mgopt'])
    imwrite(mdict['mgopt'].replace('.tiff','_invert.tiff'), mginv, metadata ={'PhysicalSizeX':"1", 'PhysicalSizeXUnit':"µm", 'PhysicalSizeY':"1", 'PhysicalSizeYUnit':"µm"}, ome=True,  bigtiff=True)
    mpinv = invert(mdict['mpopt'])
    imwrite(mdict['mpopt'].replace('.tiff','_invert.tiff'), mpinv, metadata ={'PhysicalSizeX':"1", 'PhysicalSizeXUnit':"µm", 'PhysicalSizeY':"1", 'PhysicalSizeYUnit':"µm"}, ome=True,  bigtiff=True)
    mminv = invert(mdict['mmopt'])
    imwrite(mdict['mmopt'].replace('.tiff','_invert.tiff'), mminv, metadata ={'PhysicalSizeX':"1", 'PhysicalSizeXUnit':"µm", 'PhysicalSizeY':"1", 'PhysicalSizeYUnit':"µm"}, ome=True,  bigtiff=True)

    # loads those images in the itk forma
    fixed_image = itk.imread(mdict['imcfile'], itk.F)
    moving_gly = itk.imread(mdict['mgopt'].replace('.tiff','_invert.tiff'), itk.F)
    moving_pep = itk.imread(mdict['mpopt'].replace('.tiff','_invert.tiff'), itk.F)
    moving_met = itk.imread(mdict['mmopt'].replace('.tiff','_invert.tiff'), itk.F)

    # checks to see if there are masks for the optical images, and applis them if there are
    if mdict['imcmask'] != '':
        imcu = geomask(mdict['imcmask'], mdict['imcfile'])
        fixed_image = itk.imread(imcu, itk.F)

    if mdict['mglymask'] != '':
        mglyu = geomask(mdict['mglymask'], mdict['mgopt'].replace('.tiff','_invert.tiff'))
        moving_gly = itk.imread(mglyu, itk.F)

    if mdict['mpepmask'] != '':
        mpepu = geomask(mdict['mpepmask'], mdict['mpopt'].replace('.tiff','_invert.tiff'))
        moving_pep = itk.imread(mpepu, itk.F)

    if mdict['mmetmask'] != '':
        mmetu = geomask(mdict['mmetmask'], mdict['mmopt'].replace('.tiff','_invert.tiff'))
        moving_met = itk.imread(mmetu, itk.F)

    # loads the landmarks files and generates the file needed for itk in the core directory
    imcpland = get_landmarks(mdict['imcpland'])
    imcgland = get_landmarks(mdict['imcgland'])
    imcmland = get_landmarks(mdict['imcmland'])

    mgland = get_landmarks(mdict['mgland'])
    mpland = get_landmarks(mdict['mpland'])
    mmland = get_landmarks(mdict['mmland'])

    ufile = open(mdict['imcpland'].replace('.txt','_temp.txt'), 'w')
    ufile.write(imcpland)
    ufile.close()

    ufile = open(mdict['imcgland'].replace('.txt','_temp.txt'), 'w')
    ufile.write(imcgland)
    ufile.close()

    ufile = open(mdict['imcmland'].replace('.txt','_temp.txt'), 'w')
    ufile.write(imcmland)
    ufile.close()

    ufile = open(mdict['mgland'].replace('.txt','_temp.txt'), 'w')
    ufile.write(mgland)
    ufile.close()

    ufile = open(mdict['mpland'].replace('.txt','_temp.txt'), 'w')
    ufile.write(mpland)
    ufile.close()

    ufile = open(mdict['mmland'].replace('.txt','_temp.txt'), 'w')
    ufile.write(mmland)
    ufile.close()
    
    # registers the glycan optical image to the IMC# glycans
    result_image1, result_transform_parameters1 = itk.elastix_registration_method(fixed_image, moving_gly,
        fixed_point_set_file_name=mdict['imcgland'].replace('.txt','_temp.txt'),
        moving_point_set_file_name=mdict['mgland'].replace('.txt','_temp.txt'),
        log_to_console=True,
        parameter_object=rab_param)

    glyregf.append(fixed_image)
    glyregf.append(moving_gly)
    glyregf.append(result_image1)
    glyregf.append(result_transform_parameters1)  

    # registers the peptide optical image to the IMC
    result_image1, result_transform_parameters1 = itk.elastix_registration_method(fixed_image, moving_pep,
        fixed_point_set_file_name=mdict['imcpland'].replace('.txt','_temp.txt'),
        moving_point_set_file_name=mdict['mpland'].replace('.txt','_temp.txt'),
        log_to_console=True,
        parameter_object=rab_param)

    pepregf.append(fixed_image)
    pepregf.append(moving_pep)
    pepregf.append(result_image1)
    pepregf.append(result_transform_parameters1)

    # registers the mets optical image to the IMC
    result_image1, result_transform_parameters1 = itk.elastix_registration_method(fixed_image, moving_met,
        fixed_point_set_file_name=mdict['imcmland'].replace('.txt','_temp.txt'),
        moving_point_set_file_name=mdict['mmland'].replace('.txt','_temp.txt'),
        log_to_console=True,
        parameter_object=rab_param)

    metregf.append(fixed_image)
    metregf.append(moving_met)
    metregf.append(result_image1)
    metregf.append(result_transform_parameters1)

    return glyregf, pepregf, metregf


# In[14]:


# registers the H&E to the IMC
def register_he_imc(hifiles):
    impfiles = []
    imask = False
    imcf = imread(hifiles[0])

    # loads the IMC mask if it exists and then applies it
    if hifiles[3] != '':
        imask = True
        with open(hifiles[3], 'r') as f:
            hm = geojson.load(f)
        hmarr = np.array(hm['features'][0]['geometry']['coordinates'])
        shmarr = np.squeeze(hmarr)
        hmp = Polygon(shmarr)
        
        img = np.zeros((imcf.shape[1],imcf.shape[0]), 'uint8')
        rr, cc = polygon(shmarr[:,0], shmarr[:,1], img.shape)
        img[rr,cc] = 1
        imgt = img.transpose()
       
        imcset = []
        for i in range (0,imcf.shape[2]):
            imcf2 = imcf[:,:,i]
            imcmasked = imcf2 * imgt
            imcset.append(imcmasked)
            
        imcarr  = np.array(imcset)
        imcmasked = np.transpose(imcarr, (1, 2, 0))
        imwrite(hifiles[3].replace('.geojson','_temp.tiff'), imcmasked, metadata ={'PhysicalSizeX':"1", 'PhysicalSizeXUnit':"µm", 'PhysicalSizeY':"1", 'PhysicalSizeYUnit':"µm"}, ome=True,  bigtiff=True)

    # loads the IMC file
    fixed_image = itk.imread(hifiles[0], itk.F)
    if imask:
        fixed_image = itk.imread(hifiles[3].replace('.geojson','_temp.tiff'), itk.F)

    # loads and normalizes the H&E mask (was necessary as sometimes the masks are floats and sometimes ints)
    # Final mask is either 0s or 1s
    hei = imread(hifiles[1])
    dfh = imread(hifiles[2])
    dfh = dfh/dfh.max()
    dfh = dfh.astype(np.uint8)

    # H&E image and registration files
    heset = []

    # applies H&E mask and saves a temp file
    for i in range (0,hei.shape[2]):
        heii = hei[:,:,i]
        heimasked = heii * dfh
        heset.append(heimasked)
    
    hesetarr  = np.array(heset)
    hemasked = np.transpose(hesetarr, (1, 2, 0))
    
    imwrite(hifiles[1].replace('.tiff','_temp.tiff'), hemasked, metadata ={'PhysicalSizeX':"1", 'PhysicalSizeXUnit':"µm", 'PhysicalSizeY':"1", 'PhysicalSizeYUnit':"µm"}, ome=True,  bigtiff=True)

    moving_image = itk.imread(hifiles[1].replace('.tiff','_temp.tiff'), itk.F)

    # loads the landmarks and creates itk files
    imcland = get_landmarks(hifiles[5])
    heland = get_landmarks(hifiles[4])

    ufile = open(hifiles[5].replace('.txt','_temp.txt'), 'w')
    ufile.write(imcland)
    ufile.close()
    
    hfile = open(hifiles[4].replace('.txt','_temp.txt'), 'w')
    hfile.write(heland)
    hfile.close()

    # registers the H&E to the IMC
    result_image1, result_transform_parameters1 = itk.elastix_registration_method(fixed_image, moving_image,
        fixed_point_set_file_name=hifiles[5].replace('.txt','_temp.txt'),
        moving_point_set_file_name=hifiles[4].replace('.txt','_temp.txt'),
        log_to_console=True,
        parameter_object=rab_param)

    impfiles.append(fixed_image)
    impfiles.append(moving_image)
    impfiles.append(result_image1)
    impfiles.append(result_transform_parameters1)

    impfiles.append(imread(hifiles[1].replace('.tiff','_temp.tiff')))

    return impfiles


# In[15]:


# generates landmark format for itk
def get_landmarks(lfile):
    x = []
    y = []
    lvals = ""
    
    with open(lfile, 'r') as file:
        # Read each line in the file
        for line in file:
            # Print each line
            tab = line.split("\t")
            x.append(int(float(tab[0])))
            y.append(int(float(tab[1])))
        lvals = "index\n" + str(len(x)) + '\n'
    
        for i in range (0, len(x)):
            lvals += str(x[i]) + ' ' + str(y[i]) + '\n'
    return lvals


# In[16]:


# Rigid-Affine-Bspline registration
# uses landmarking for the R and A
def rab_registration():
    parameter_object = itk.ParameterObject.New()
    parameter_map_rigid = parameter_object.GetDefaultParameterMap('rigid')
    parameter_map_rigid['Registration'] = ['MultiMetricMultiResolutionRegistration']
    original_metric = parameter_map_rigid['Metric']
    parameter_map_rigid['Metric'] = [original_metric[0], 'CorrespondingPointsEuclideanDistanceMetric']
    #parameter_map_rigid['Metric0Weight'] = ['0.0']
    #parameter_map_rigid['Metric1Weight'] = ['1.0']
    parameter_object.AddParameterMap(parameter_map_rigid)
    parameter_map_affine = parameter_object.GetDefaultParameterMap('affine')
    parameter_map_affine['Registration'] = ['MultiMetricMultiResolutionRegistration']
    original_metric = parameter_map_rigid['Metric']
    parameter_map_affine['Metric'] = [original_metric[0], 'CorrespondingPointsEuclideanDistanceMetric']
    parameter_object.AddParameterMap(parameter_map_affine)
    parameter_object.SetParameter(1, "MaximumNumberOfIterations", "1000")
    parameter_map_bspline = parameter_object.GetDefaultParameterMap('bspline')
    parameter_object.AddParameterMap(parameter_map_bspline)
    parameter_object.SetParameter(2, "MaximumNumberOfIterations", "1000")
    parameter_object.SetParameter('ResampleInterpolator','FinalNearestNeighborInterpolator')

    return parameter_object


# In[17]:


# Rigid-Bspline registration, not used in the current version as performed worse than RAB
# Rigid-Affine wasn't made because while testing it cause large distortions if 
# too much rotation was needed or there was too much difference in shape/size
def rb_registration():
    parameter_object = itk.ParameterObject.New()
    parameter_map_rigid = parameter_object.GetDefaultParameterMap('rigid')
    parameter_map_rigid['Registration'] = ['MultiMetricMultiResolutionRegistration']
    original_metric = parameter_map_rigid['Metric']
    parameter_map_rigid['Metric'] = [original_metric[0], 'CorrespondingPointsEuclideanDistanceMetric']
    parameter_object.AddParameterMap(parameter_map_rigid)
    parameter_map_bspline = parameter_object.GetDefaultParameterMap('bspline')
    parameter_object.AddParameterMap(parameter_map_bspline)
    parameter_object.SetParameter(1, "MaximumNumberOfIterations", "1000")
    parameter_object.SetParameter('ResampleInterpolator','FinalNearestNeighborInterpolator')
    
    return parameter_object


# In[18]:


# Applies the transform to another image/layer
# This one is for when the image and the transform are done at the same resolution
def get_transform(image, transform_params):
    comb = []
    imgs = 1

    # if there is more than 1 layer in the image, apply the transform to each layer
    if len(image.shape) >2:
        imgs = image.shape[2]
    
        for i in range (0,imgs):
            imagec = image[:,:,i]
            mit = itk.GetImageFromArray(imagec)
        
            transformix_object = itk.TransformixFilter.New(mit)
            transformix_object.SetTransformParameterObject(transform_params)
            
            # Update object (required)
            transformix_object.UpdateLargestPossibleRegion()
            
            # Results of Transformation
            result_image_transformix = transformix_object.GetOutput()
        
            np_view = itk.array_view_from_image(result_image_transformix)
        
            comb.append(np_view)
    # otherwise apply the transform to the whole image
    else:
        mit = itk.GetImageFromArray(image)
    
        transformix_object = itk.TransformixFilter.New(mit)
        transformix_object.SetTransformParameterObject(transform_params)
        
        # Update object (required)
        transformix_object.UpdateLargestPossibleRegion()
        
        # Results of Transformation
        result_image_transformix = transformix_object.GetOutput()
    
        np_view = itk.array_view_from_image(result_image_transformix)
    
        comb.append(np_view)
    return comb


# In[19]:


# Applies the transform to another image/layer
# This one is for when the image is 20um and the transform was done at 1um
def get_transform_msi(image, transform_params):
    comb = []
    imgs = 1

    # if there is more than 1 layer in the image, apply the transform to each layer
    if len(image.shape) >2:
        imgs = image.shape[2]
    
        for i in range (0,imgs):
            imagec = image[:,:,i]
            # expands 20umx20um into 1umx1um by duplicating pixels
            imagec = imagec.repeat(20, axis=0).repeat(20, axis=1)
            mit = itk.GetImageFromArray(imagec)
        
            transformix_object = itk.TransformixFilter.New(mit)
            transformix_object.SetTransformParameterObject(transform_params)
            
            # Update object (required)
            transformix_object.UpdateLargestPossibleRegion()
            
            # Results of Transformation
            result_image_transformix = transformix_object.GetOutput()
        
            np_view = itk.array_view_from_image(result_image_transformix)
        
            comb.append(np_view)
    # otherwise apply the transform to the whole image
    else:
        mit = itk.GetImageFromArray(image)
    
        transformix_object = itk.TransformixFilter.New(mit)
        transformix_object.SetTransformParameterObject(transform_params)
        
        # Update object (required)
        transformix_object.UpdateLargestPossibleRegion()
        
        # Results of Transformation
        result_image_transformix = transformix_object.GetOutput()
    
        np_view = itk.array_view_from_image(result_image_transformix)
    
        comb.append(np_view)
    return comb


# In[ ]:





# In[20]:


# initializes the Rigid-Affine-Bspline registration object
rab_param = rab_registration()


# In[21]:


# list of masks that you will register
# Since focal four was named F4, you will get issues with cores named F4 
# (in the future just make sure that your cores and your masks can't be named the same thing; you can fix this by moving/renaming the base image for this core)
pannot = ['_2.tiff','_3.tiff','_4.tiff','_5.tiff', '_B.tiff', '_S.tiff','_W.tiff','_3G.tiff','_3S.tiff','_3M.tiff','_F4.tiff']


# In[22]:


# metadata list that passes information between functions
# opt = optical image
# land is landmark
# mask is masks; not all masks have to be present as they are mostly used for removing damaged tissue
# if it starts with m its MSI, otherwise m is mets, g is glycans, p is peptides
mlorder = ['mgopt','mpopt','mmopt','mgland','mpland','mmland','imcgland','imcpland','imcmland','imcfile','mglymask','mpepmask','mmetmask','imcmask']


# In[23]:


# iterates through all cores to perform H&E-IMC registration
with os.scandir(dirimc) as entries:
    for entry in entries:
        if entry.is_dir():
            # checkpoint to see if the registration was already performed
            # if you want to reregister a ficore and this exists, just delete it
            # or change the 1 to something like 2
            # is useful so you dont have to reregister everything if/when something fails
            if len(glob.glob(entry.path + "/" + 'complete.out')) <1:
                # 3 channel IMC
                imcfile = glob.glob(entry.path + "/" + '*imc3d.tiff')  
                # H&E file that has background brightness removed
                hefile = glob.glob(entry.path + "/" + '*HE_masked.tiff') 
                # core mask for the H&E
                hemask = glob.glob(entry.path + "/" + '*W.tiff') 
                # mask for IMC, if present
                imcmask = glob.glob(entry.path + "/" + '*imc3d_mask.geojson')
                # landmark files
                heland = glob.glob(entry.path + "/" + '*HE_masked_landmarks.txt')
                imcland = glob.glob(entry.path + "/" + '*imc3d_landmarks.txt')

                # identifies which files are present (as not all masks have to be)
                filelist = [imcfile,hefile,hemask,imcmask,heland,imcland]
    
                hifiles = he_imc_files(filelist)
    
                if len(hifiles) > 0:
                    # registers H&E to 3 channel IMC
                    regfiles = register_he_imc(hifiles)
                    timage1 = get_transform(regfiles[4],regfiles[3])
                    u_reshaped = np.transpose(timage1, (1, 2, 0))
                    imwrite(hefile[0].replace('.tiff','_rab_reg.tiff'), u_reshaped, metadata ={'PhysicalSizeX':"1", 'PhysicalSizeXUnit':"µm", 'PhysicalSizeY':"1", 'PhysicalSizeYUnit':"µm"}, bigtiff=True, ome=True, dtype=np.uint8)
                    plot_reg(regfiles[2], regfiles[1], regfiles[0], hefile[0].replace('.tiff','_rab_reg_image.jpeg'))

                    # for each of the masks in the H&E transform them
                    for item in pannot:
                        fitem = glob.glob(entry.path + "/" + '*' + item)
                        if len(fitem) > 0:
                            mf = imread(fitem[0])
                            mf = mf/mf.max()
                            mf = mf.astype(np.uint8)
                            timagea = get_transform(mf,regfiles[3])
                            u_reshaped = np.transpose(timagea, (1, 2, 0))
                            imwrite(fitem[0].replace('.tiff','_rab_reg.tiff'), u_reshaped, metadata ={'PhysicalSizeX':"1", 'PhysicalSizeXUnit':"µm", 'PhysicalSizeY':"1", 'PhysicalSizeYUnit':"µm"}, bigtiff=True, ome=True, dtype=np.uint8)

                    # creates an empty file that denotes that transformation was completed
                    with open(entry.path + "/" + 'complete.out', 'a') as file:
                        file.write("completed.\n")


# In[28]:


# iterates through each core and performs MSI-IMC registration
# uses the presence of the IMC as to whether or not the core is used
with os.scandir(dirimc) as entries:
    for entry in entries:
        if entry.is_dir():
            # MSI and IMC had different IDs and this maps between them
            name = entry.name.split('PR1211_181_')[1]
            dirn = glob.glob(dirmsi + "/" + 'PR1121_180_' + name)[0]
            # checkpoint to see if the registration was already performed
            # if you want to reregister a ficore and this exists, just delete it
            # or change the 1 to something like 2
            # is useful so you dont have to reregister everything if/when something fails
            if len(glob.glob(dirn + "/" + 'complete.out')) <1:
                # 3-channel IMC image
                imcfile = glob.glob(entry.path + "/" + '*imc3d.tiff')  
                # mask file, if present
                imcmask = glob.glob(entry.path + "/" + '*imc3d_mask.geojson')

                # optical images at 1um
                mgopt = glob.glob(dirn + "/*glycans*optical_1um.tiff")
                mpopt = glob.glob(dirn + "/*peptides*optical_1um.tiff")
                mmopt = glob.glob(dirn + "/*mets*optical_1um.tiff")
                # optical image landmarks
                mgland = glob.glob(dirn + "/*glycans*optical_landmarks.txt")
                mpland = glob.glob(dirn + "/*peptides*optical_landmarks.txt")
                mmland = glob.glob(dirn + "/*mets*optical_landmarks.txt")

                # msi images at 20 um
                mgly = glob.glob(dirn + "/*glycans_msi.tiff")
                mpep = glob.glob(dirn + "/*peptides_msi.tiff")
                mmet = glob.glob(dirn + "/*mets_msi.tiff")

                # imc landmarks
                imcgland = glob.glob(dirn + "/*glycans_landmarks.txt")
                imcpland = glob.glob(dirn + "/*peptides_landmarks.txt")
                imcmland = glob.glob(dirn + "/*mets_landmarks.txt")

                # msi mask, if present
                mglymask = glob.glob(dirn + "/*glycans_mask.geojson")
                mpepmask = glob.glob(dirn + "/*peptides_mask.geojson")
                mmetmask = glob.glob(dirn + "/*mets_mask.geojson")

                # identifies which files are present (as not all masks have to be)
                mflist = [mgopt,mpopt,mmopt,mgland,mpland,mmland,imcgland,imcpland,imcmland,imcfile,mglymask,mpepmask,mmetmask,imcmask]
                mdict = msi_imc_files(mflist)
    
                if len(mdict) > 1:
                    glyregf, pepregf, metregf = register_msi_imc(mdict)
    
                    # registers the glycan optical image to IMC
                    timage1 = get_transform(glyregf[1],glyregf[3])
                    u_reshaped = np.transpose(timage1, (1, 2, 0))
                    u_reshaped = u_reshaped.astype(np.uint8)
                    imwrite(mgopt[0].replace('.tiff','_rab_reg.tiff'), u_reshaped, metadata ={'PhysicalSizeX':"1", 'PhysicalSizeXUnit':"µm", 'PhysicalSizeY':"1", 'PhysicalSizeYUnit':"µm"}, bigtiff=True, ome=True, dtype=np.uint8)
                    plot_reg(glyregf[2], glyregf[1], glyregf[0], mgopt[0].replace('.tiff','_rab_reg_image.jpeg'))
    
                    if len(mgly) > 0:
                        # obtains m/z metadata information for glycans
                        # not necessary to do every time as they will be the same tiff to tiff
                        # if you ever use imzmls you will have to
                        mgnames = []
                        with TiffFile(mgly[0]) as tif:
                            ome_xml = tif.ome_metadata
                        
                            ome = from_xml(ome_xml)
                        
                            # Access metadata elements
                            for channel in ome.images[0].pixels.channels:
                                mgnames.append(channel.name)

                        # transforms glycan MSI images based on optical image transformation
                        mgf = imread(mgly[0])
                        mgf = np.transpose(mgf, (1, 2, 0))
                        timagea = get_transform_msi(mgf,glyregf[3])
                        # exists to format to numpy instead of itk
                        u_reshaped = np.transpose(timagea, (0, 1, 2))
                        imwrite(mgly[0].replace('.tiff','_rab_reg.tiff'), u_reshaped, metadata ={'PhysicalSizeX':"1", 'PhysicalSizeXUnit':"µm", 'PhysicalSizeY':"1", 'PhysicalSizeYUnit':"µm", "Channel": {"Name": mgnames}}, bigtiff=True, ome=True, dtype=np.float32, compression ='zlib')
    
                    # registers peptide optical image to IMC
                    timage1 = get_transform(pepregf[1],pepregf[3])
                    u_reshaped = np.transpose(timage1, (1, 2, 0))
                    u_reshaped = u_reshaped.astype(np.uint8)
                    imwrite(mpopt[0].replace('.tiff','_rab_reg.tiff'), u_reshaped, metadata ={'PhysicalSizeX':"1", 'PhysicalSizeXUnit':"µm", 'PhysicalSizeY':"1", 'PhysicalSizeYUnit':"µm"}, bigtiff=True, ome=True, dtype=np.uint8)
                    plot_reg(pepregf[2], pepregf[1], pepregf[0], mpopt[0].replace('.tiff','_rab_reg_image.jpeg'))
    
                    if len(mpep) > 0:
                        # obtains m/z metadata information for peptides
                        # not necessary to do every time as they will be the same tiff to tiff
                        # if you ever use imzmls you will have to
                        mpnames = []
                        with TiffFile(mpep[0]) as tif:
                            ome_xml = tif.ome_metadata
                        
                            ome = from_xml(ome_xml)
                        
                            # Access metadata elements
                            for channel in ome.images[0].pixels.channels:
                                mpnames.append(channel.name)

                        # transforms peptide MSI images based on optical image transformation
                        mpf = imread(mpep[0])
                        mpf = np.transpose(mpf, (1, 2, 0))
                        timagea = get_transform_msi(mpf,pepregf[3])
                        # exists to format to numpy instead of itk
                        u_reshaped = np.transpose(timagea, (0, 1, 2))
                        imwrite(mpep[0].replace('.tiff','_rab_reg.tiff'), u_reshaped, metadata ={'PhysicalSizeX':"1", 'PhysicalSizeXUnit':"µm", 'PhysicalSizeY':"1", 'PhysicalSizeYUnit':"µm", "Channel": {"Name": mpnames}}, bigtiff=True, ome=True, dtype=np.float32, compression ='zlib')
                       
                    # transforms mets optical image to IMC
                    timage1 = get_transform(metregf[1],metregf[3])
                    u_reshaped = np.transpose(timage1, (1, 2, 0))
                    u_reshaped = u_reshaped.astype(np.uint8)
                    imwrite(mmopt[0].replace('.tiff','_rab_reg.tiff'), u_reshaped, metadata ={'PhysicalSizeX':"1", 'PhysicalSizeXUnit':"µm", 'PhysicalSizeY':"1", 'PhysicalSizeYUnit':"µm"}, bigtiff=True, ome=True, dtype=np.uint8)
                    plot_reg(metregf[2], metregf[1], metregf[0], mmopt[0].replace('.tiff','_rab_reg_image.jpeg'))
                  
                    if len(mmet) > 0:
                        # obtains m/z metadata information for mets
                        # not necessary to do every time as they will be the same tiff to tiff
                        # if you ever use imzmls you will have to
                        mmnames = []
                        with TiffFile(mmet[0]) as tif:
                            ome_xml = tif.ome_metadata
                        
                            ome = from_xml(ome_xml)
                        
                            # Access metadata elements
                            for channel in ome.images[0].pixels.channels:
                                mmnames.append(channel.name)
                                
                        # transforms mets MSI images based on optical image transformation
                        mmf = imread(mmet[0])
                        mmf = np.transpose(mmf, (1, 2, 0))
                        timagea = get_transform_msi(mmf,metregf[3])
                        # exists to format to numpy instead of itk
                        u_reshaped = np.transpose(timagea, (0, 1, 2))
                        imwrite(mmet[0].replace('.tiff','_rab_reg.tiff'), u_reshaped, metadata ={'PhysicalSizeX':"1", 'PhysicalSizeXUnit':"µm", 'PhysicalSizeY':"1", 'PhysicalSizeYUnit':"µm", "Channel": {"Name": mmnames}}, bigtiff=True, ome=True, dtype=np.float32, compression ='zlib')
                       
                    # creates an empty file that denotes that transformation was completed
                    with open(dirn + "/" + 'complete.out', 'a') as file:
                        file.write("completed.\n")


# In[ ]:




