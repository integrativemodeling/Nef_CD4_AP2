from modeller import *
from modeller.automodel import *

# Override the 'special_restraints' and 'user_after_single_model' methods:

class MyModel(automodel):
    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms
    
    def special_patches(self, aln):
        self.rename_segments(segment_ids=('A','S','N','B','M','D'), renumber_residues=[9, 1, 33, 1, 1, 405])

    def select_atoms(self):
        # Select  residues 1 and 2 (PDB numbering)
        sele = selection(self.residue_range('1:B', '13:B'),
                         self.residue_range('25:B', '87:B'),
                         self.residue_range('40:N', '48:N'),
                         self.residue_range('73:N', '78:N'),
                         self.residue_range('110:N', '125:N'),
                         self.residue_range('1:M', '5:M'),
                         self.residue_range('21:M', '26:M'),
                         self.residue_range('38:M', '45:M'))
        return sele

env = environ()

env.io.hetatm = True

a = MyModel(env, alnfile='aln_2vgl.pir', knowns=('Fixed_sequence_refined','2vgl_B_nohelix'),sequence='mod_AP2_Nef',assess_methods=(assess.DOPE,assess.GA341, assess.normalized_dope))

a.starting_model = 1
a.ending_model = 10
a.make()                           # do comparative modeling


