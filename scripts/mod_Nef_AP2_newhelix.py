############################################
# Modeling the Nef-CD4_AP2 complex
#
# iecheverria - Sali Lab - UCSF
# ignacia@salilab.org
############################################
import IMP
import RMF
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.occams
from IMP.pmi.io.crosslink import CrossLinkDataBaseKeywordsConverter
import IMP.pmi.mmcif

import random
import numpy as np
import sys

###################### SYSTEM SETUP #####################

include_Occams = True

mdl = IMP.Model()

###############################
# Setup system
###############################
top_spec1 = 'top_Nef_AP2_newhelices.dat'
reader_spec1 = IMP.pmi.topology.TopologyReader(top_spec1,
                                               pdb_dir ='../data/',
                                               fasta_dir ='../data/')

bs_spec1 = IMP.pmi.macros.BuildSystem(mdl,
                                      resolutions=[1])
    
##############################
# Generate mmcif file
##############################

if '--mmcif' in sys.argv:
    print('Hola')
    # Record the modeling protocol to an mmCIF file
    po = IMP.pmi.mmcif.ProtocolOutput(open('IM_Nef-CD4-AP2.cif', 'w'))
    po.system.title = ('Structural Basis of CD4 Downregulation by HIV-1 Nef')
    bs_spec1.system.add_protocol_output(po)
    
    # Add publication
    #po.system.citations.append(ihm.Citation.from_pubmed_id(0000))

##############################
# Build state
##############################

bs_spec1.add_state(reader_spec1)

hier,  dof = bs_spec1.execute_macro(max_rb_trans=3.0,
                                    max_rb_rot=0.03)
mols_S1 = bs_spec1.get_molecules()[0]


# Write coordinates
out = IMP.pmi.output.Output()
out.init_rmf("all_ini.rmf3", [hier])
out.write_rmf("all_ini.rmf3")
out.close_rmf("all_ini.rmf3")

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

# Add distance restraint between CD4 and Nef
dr = IMP.pmi.restraints.basic.DistanceRestraint(root_hier = hier,
                                                tuple_selection1=(222,222,"Nef"),
                                                tuple_selection2=(394,394,"CD4mut"),
                                                distancemin=3,
                                                distancemax=5)
dr.add_to_model()

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
cldb.create_set_from_file("../data/Interlinks_Nef_AP2_20190723_renamed_renumbered_nonambiguos.csv")

xl1 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=hier,
                                                                            CrossLinkDataBase=cldb,
                                                                            resolution=1.0,
                                                                            length=21.0,
                                                                            slope=0.02)
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
    occ.set_weight(1.5)
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
dof.optimize_flexible_beads(200)
  
############################# SAMPLING ##############################
# Run replica exchange Monte Carlo sampling


num_frames = 50000
if '--mmcif' in sys.argv or '--test' in sys.argv:
    num_frames=5

rex=IMP.pmi.macros.ReplicaExchange0(mdl,
                                    root_hier=hier,                          
                                    crosslink_restraints=rmf_restraints,           
                                    monte_carlo_sample_objects=dof.get_movers(),
                                    replica_exchange_maximum_temperature=3.0,
                                    global_output_directory="output/",
                                    output_objects=output_objects,
                                    monte_carlo_steps=10,
                                    number_of_frames=num_frames,
                                    number_of_best_scoring_models=10)

rex.execute_macro()

##############################
# Generate mmcif
##############################  
if '--mmcif' in sys.argv:
    import ihm.cross_linkers
    import ihm.dumper
    import ihm.format
    import ihm.location
    import ihm.representation
    import ihm.startmodel
    import ihm.dataset
    import ihm.protocol
    import ihm.analysis
    import ihm.model
    import ihm.restraint
    import ihm.geometry
    
    fname = '../data/Interlinks_Nef_AP2_20190723_renamed_renumbered_nonambiguos.csv'
        
    s = po.system
    print("restraint datasets:", [r.dataset for r in s.restraints])
    # Datasets for XL-MS restraint
    for r in s.restraints:
        if isinstance(r, ihm.restraint.CrossLinkRestraint):
            r.linker = ihm.cross_linkers.dsso
            print("XL-MS dataset at:", r.dataset.location.path)
            print("Details:", r.dataset.location.details)
    
    # Correct number of output models to account for multiple runs
    protocol = po.system.orphan_protocols[-1]
    protocol.steps[-1].num_models_end = 2007800

    # Get last protocol in the file
    protocol = po.system.orphan_protocols[-1]
    # State that we filtered the 200000 frames down to one cluster of
    # 9999 models:
    analysis = ihm.analysis.Analysis()
    protocol.analyses.append(analysis)
    analysis.steps.append(ihm.analysis.ClusterStep(
                            feature='RMSD', num_models_begin=200000,
                            num_models_end=9999))
    
    # Create an ensemble for the cluster
    e = po._add_simple_ensemble(analysis.steps[-1],
                                name="Cluster 0", num_models=9999,
                                drmsd=8.3, num_models_deposited=1,
                                localization_densities={}, ensemble_file=None)
    
    # Add the model from RMF
    rh = RMF.open_rmf_file_read_only('../results/clustering/cluster.0/cluster_center_model.rmf3')
    IMP.rmf.link_hierarchies(rh, [hier])
    IMP.rmf.load_frame(rh, RMF.FrameID(0))
    del rh
    model = po.add_model(e.model_group)

    # Add localization densities
    # Look up the ihm.AsymUnit corresponding to a PMI component name
    for asym in po.asym_units:
        name = asym.split('.')[0]
        fname = f'../results/clustering/cluster.0/LPD_{name}_orie.mrc'
        print('fname', fname)
        loc = ihm.location.OutputFileLocation(fname)
        den = ihm.model.LocalizationDensity(file=loc, asym_unit=po.asym_units[asym])
        # Add to ensemble
        e.densities.append(den)
        

    # Replace local links with DOIs
    #repo = ihm.location.Repository(doi="10.5281/zenodo.2598760", root="../..",
    #              top_directory="salilab-imp_deposition_tutorial-1ad5919",
    #              url="https://zenodo.org/record/2598760/files/salilab/"
    #                  "imp_deposition_tutorial-v0.2.zip")
    #po.system.update_locations_in_repositories([repo])

po.flush()



exit()

