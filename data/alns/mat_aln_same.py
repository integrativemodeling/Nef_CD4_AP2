import pandas as pd
import numpy as np
import os, sys
from sys import exit
import pprint
import random
import copy
import glob
import Bio
from Bio import SeqIO

seq_file = '../Nef_AP2.fasta'
seqs = SeqIO.parse(seq_file,'fasta')
for seq in seqs:
    s = seq.seq
    name = seq.name
    aln = []
    for i, rt in enumerate(s):
        aln.append([rt, name, i+1, rt, name, i+1, 9.0])
    np.savetxt('mat_'+name+'_'+name+'.align', np.array(aln), fmt='%s')

exit()

