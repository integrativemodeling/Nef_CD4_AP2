#!/usr/bin/env python

'''
Different distance metrics to compare
structures in a ensemble
'''

import os
import random
import math
import itertools
import numpy as np

import sys

sys.path.append('/home/ignacia/SOFTW/PMI_analysis/pyext/src/')
from compute_distance_metrics import *

sel = int(sys.argv[1])
print(sel)

alignment = [('Nef',21,202),
             ('CD4mut',394,419),
             ('AP2mu2',1,135),
             ('AP2sigma',1,142),
             ('AP2beta2',1,615),
             ('AP2beta2_s',90,615),
             ('AP2alpha2',1,641)]



alignment = [('AP2beta2', 91,583),
             ('AP2beta2', 29,43),
             ('AP2beta2', 49,61),
             ('AP2beta2', 64,78),
             ('AP2beta2', 81,86),
             ('AP2alpha2', 9,619), 
             ('Nef', 33,202)]


alignment_sel = [alignment[sel]]
print('Running (alignment): ', alignment_sel)

D = get_distance_metrics('analys/clustering_cl-1/',
                        0,
                        align_to = alignment_sel,
                        out_dir_name = 'RMSDs_5K_%s_%s'%(alignment_sel[0][0],alignment_sel[0][1]),
                        number_of_models = 5000)

D.compute_RMSD_all_versus_centroid()
del(D)

# nohup sh -c 'for i in {0..4}; do python get_RMSD_all.py $i > rmsd_$i.log; done > log' & 



