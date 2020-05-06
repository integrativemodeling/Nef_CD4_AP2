############################################
# Modeling the Nef-CD4_AP2 complex
#
# iecheverria - Sali Lab - UCSF
############################################
import IMP
import RMF
import IMP.atom
#import IMP.rmf
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.em
import IMP.pmi.restraints.occams
from IMP.pmi.io.crosslink import CrossLinkDataBaseKeywordsConverter

import random
import numpy as np
import glob
from sys import exit
from sys import argv

###################### SYSTEM SETUP #####################
include_Occams = True

mdl = IMP.Model()

###############################
# Species 1
###############################
top_spec1 = top_dir+'top_Nef_AP2_4rbs.dat'
reader_spec1 = IMP.pmi.topology.TopologyReader(top_spec1,
                                               pdb_dir = '../data/',
                                               fasta_dir = '../data/')

bs_spec1 = IMP.pmi.macros.BuildSystem(mdl,
                                      resolutions=[1])
bs_spec1.add_state(reader_spec1)

hier,  dof = bs_spec1.execute_macro(max_rb_trans=3.0,
                                          max_rb_rot=0.03)
mols_S1 = bs_spec1.get_molecules()[0]

##############################
# Connectivity
##############################
output_objects = [] # keep a list of functions that need to be reported
sample_objects = []
rmf_restraints = []

crs = []
for molname in mols_S1:
    for mol in mols_S1[molname]:
        copy_n = IMP.atom.Copy(mol.get_hierarchy()).get_copy_index()
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
        cr.set_label(mol.get_name()+'.'+str(copy_n))
        cr.add_to_model()
        output_objects.append(cr)
        crs.append(cr)


##############################
# Excluded Volume
##############################
evr1 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=mols_S1.values(),
                                                               resolution=10)
evr1.add_to_model()
evr1.set_weight(1.0)
output_objects.append(evr1)

##############################
# Cross-links
##############################
# INITIALIZE DB    

cldbkc=IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
cldbkc.set_protein1_key("Protein1")
cldbkc.set_protein2_key("Protein2")
cldbkc.set_residue1_key("AbsPos1")
cldbkc.set_residue2_key("AbsPos2")
cldbkc.set_unique_id_key("Id")
cldbkc.set_psi_key("Score")

# XLs RESTRAINT
cldb=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
cldb.create_set_from_file("../data/Interlinks_Nef_AP2_20190116_unique_modeling.csv")

xl1 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=hier,
                                                                            CrossLinkDataBase=cldb,
                                                                            resolution=1.0,
                                                                            length=21.0,
                                                                            slope=0.05)
xl1.add_to_model()
xl1.set_weight(1.0)

rmf_restraints.append(xl1)
output_objects.append(xl1)
dof.get_nuisances_from_restraint(xl1)

##############################
# Occams Spring restaint
##############################
if include_Occams:
    occ = IMP.pmi.restraints.occams.OccamsRestraint(hier,
                                                    hier,
                                                    'equiv_assis_self.dat',
                                                    '../data/alns',
                                                    sample_sys_1 = False,
                                                    sigma_init=8.0,
                                                    slope=0.0005,
                                                    psi_nuisances_are_optimized=True,
                                                    sigma_nuisances_are_optimized=True)
    occ.add_to_model()
    occ.set_weight(2.0)
    occ.write_distances('')
    rmf_restraints.append(occ)
    output_objects.append(occ)
    
    print('nuisances:', occ.get_particles_to_sample())
    dof.get_nuisances_from_restraint(occ)

    socc = occ.get_output()
    print(socc)
        
##############################
# Shuffle
##############################    
IMP.pmi.tools.shuffle_configuration(hier,
                                    max_translation=50)
#                                    bounding_box=((-100,-100,-100),(100,100,150)))
dof.optimize_flexible_beads(200)

socc = occ.get_output()
print(socc)
print(xl1.get_output()['CrossLinkingMassSpectrometryRestraint_Data_Score'])
  
############################# SAMPLING ##############################
# Run replica exchange Monte Carlo sampling
rex=IMP.pmi.macros.ReplicaExchange0(mdl,
                                    root_hier=hier,                          
                                    crosslink_restraints=rmf_restraints,           
                                    monte_carlo_sample_objects=dof.get_movers(),
                                    replica_exchange_maximum_temperature=3.0,
                                    global_output_directory="output/",
                                    output_objects=output_objects,
                                    monte_carlo_steps=10,
                                    number_of_frames=50000,
                                    number_of_best_scoring_models=10)

rex.execute_macro()
exit()

