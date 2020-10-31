#!/usr/bin/env python
# coding: utf-8

# # Calculate CS shiftx2

# In[ ]:


import MDAnalysis as mda
import sys
import shutil
import pandas as pd
sys.path.append("..")
from src.features.build_features import get_chemical_shifts
import os
import numpy as np
# from spc_imports import *
# set_up_plt()


# In[ ]:


raw_data_dir = '../data/raw/'
interim_data_dir = '../data/interim/'
processed_data_dir = '../data/processed/'
external_data_dir = '../data/external/'


# In[ ]:


states = { 
#               '5VKH_lb': {'begin': 0,
#                       'end': 1.e+20},
#     '3FB5_lb' : {'begin': 400000.,
#                       'end': 1000000.},
#           '5VK6_lb': {'begin': 0.,
#                       'end': 350000.},
          '5VKE_lb': {'begin': 0,
                      'end': 1000000.}
         }


# ## How to run this Notebook
# 
# This notebook must be used after processing the trajectories with `calculate_cs.ipynb`
# 
# This notebook must be run with the environment created with: `environment_for_shiftx2.yml`. This is because shiftx2 is hard to install and use otherwise.
# 
# Shiftx2 was written in python2. The ominia channel version of shiftx2 for python3.5 is essentially only a wrapper calling shiftx2 with python2 from python3.5.
# In order to not have any problem of python3 being picked up instead of python2, it is best to:
# ```bash
# cd /path/to/anaconda3/envs/nmr_assign_state_for_shiftx2/share/shiftx2
# find . -iname "[a-zA-Z]*.py" -exec sed -i '1s/python/python2/' {} \;
# find . -iname "[a-zA-Z]*.py" -exec sed -i 's/"python/"python2/' {} \;
# ```
# 
# 
# 
# 

# ### Get Chemical Shifts

# In[ ]:


for method in ['shiftx2']:
    for state in states.keys():
        univ = mda.Universe(interim_data_dir + state+ '/protein_sk1_pbc.pdb',
                            interim_data_dir + state+'/protein_sk1_pbc.xtc')
        df = get_chemical_shifts(univ, '../data/processed/',method=method,pH=4.,skip=1, protein_selection='protein and not resid 26 121')
        df.to_pickle(processed_data_dir+state+'/CS_'+method+'_'+state+'.pkl')


# In[ ]:




