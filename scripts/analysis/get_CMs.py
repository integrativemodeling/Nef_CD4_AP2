import IMP
import IMP.pmi
import IMP.pmi.analysis
import sys
import os

sys.path.append('/home/ignacia/SOFTW/PMI_analysis/pyext/src/')
from contact_maps import CMTable

import glob
import numpy as np
import random

 

#####################################################
# calculate contact frequency
#####################################################
CM = CMTable(GSMs_dir = 'analys/GSMs_-1/',
             clustering_dir = 'analys/clustering_cl-1/',
             out_dir = 'analys/clustering_cl-1/CMs/',
             cluster = 0,
             number_of_models = 10000,
             selection = ['Nef','CD4mut','AP2alpha1','AP2mu2','AP2sigma','AP2beta2'],    
             cutoff = 14.0,
             nproc = 20)
CM.add_XLs("Interlinks_Nef_AP2_20190723_renamed_renumbered_nonambiguos.csv")
CM.compute_contact_maps()

#CM.get_close_contacts()
CM.plot_contact_maps()

