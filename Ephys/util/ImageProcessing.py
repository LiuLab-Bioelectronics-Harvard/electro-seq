#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from os import path
import sys
sys.path.append("./")
from features_ECG import *
import numpy as np
import widgets as wd
import pywt
import matplotlib.pyplot as plt
import pandas as pd

def contact_cell_idx_finding(meta_file, save_pkl,plot,search_method,data_folder,cell_mask_folder,
                             gene_typing_cluster,cardiomypcyte_cluster,radius):
    if path.exists(save_pkl):
        cell_recorded = wd.load_obj(save_pkl)
    else:
        cell_recorded = wd.contact_cell_finding(data_folder,cell_mask_folder,meta_file,save_pkl,
                                                plot,search_method,gene_typing_cluster,cardiomypcyte_cluster,radius)
    stage_contact_idx = []
    cell_meta = pd.read_csv(meta_file, header=0)
    for i in range(len(cell_recorded['position'])):
        position = cell_recorded['position'][i]
        cell_id = cell_recorded['cell_id'][i]
        if not position:
            continue
        if cell_id:
            for sub_cell_id in cell_id:
                # sub_cell_id = cell_id[0]
                position_idx = np.array(cell_meta['position'] == position[0])
                cell_id_idx = np.array(cell_meta['Unnamed: 0'] == sub_cell_id)
                cell_idx = np.logical_and(position_idx, cell_id_idx)
                idx = cell_meta.index[cell_idx].tolist()
                stage_contact_idx.append(idx)
    return stage_contact_idx

