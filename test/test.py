import unittest
import os
import shutil
import sys
import subprocess
import ihm.reader

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '..'))


class Tests(unittest.TestCase):

    def run_modeller_script(self, script_name, model_name, resrng):
        """Run a Modeller script and test the output model"""
        os.chdir(os.path.join(TOPDIR, 'data'))
        # Run script
        p = subprocess.check_call(["python", script_name, "--test"])
        # Make sure PDB was produced with the requested residue range
        with open('%s.B99990001.pdb' % model_name) as fh:
            pdb_lines = [x for x in fh.readlines() if x.startswith('ATOM')]
        rng = (int(pdb_lines[0][22:26]), int(pdb_lines[-1][22:26]))
        self.assertEqual(rng, resrng)

    def test_ap2beta2(self):
        """Test generation of AP2 beta2 model with Modeller"""
        self.run_modeller_script('mod_extended_AP2beta2.py',
                                 'mod_AP2_Nef', (9, 418))

    def test_simple(self):
        """Test model building"""
        os.chdir(os.path.join(TOPDIR, 'scripts'))
        p = subprocess.check_call(["python", "mod_Nef_AP2_newhelix.py",
                                   "--test"])
        # todo: assert outputs
        
    def test_mmcif(self):
        """Test generation of mmCIF output"""
        os.chdir(os.path.join(TOPDIR, 'scripts'))
        if os.path.exists("IM_Nef-CD4-AP2.cif"):
            os.unlink("IM_Nef-CD4-AP2.cif")
        # Potentially override methods that need network access
        env = os.environ.copy()
        env['PYTHONPATH'] = os.path.join(TOPDIR, 'test', 'mock') \
                            + ':' + env.get('PYTHONPATH', '')
        p = subprocess.check_call(
                ["python", "mod_Nef_AP2_newhelix.py", "--mmcif", "--test"],
                env=env)
        # Check output file
        self._check_mmcif_file('IM_Nef-CD4-AP2.cif')

    def _check_mmcif_file(self, fname):
        with open(fname) as fh:
            s, = ihm.reader.read(fh)
        self.assertEqual(len(s.citations), 1)
        self.assertEqual(s.citations[0].doi, '10.1038/s41594-020-0463-z')
        self.assertEqual(len(s.software), 3)
        self.assertEqual(len(s.orphan_starting_models), 6)
        # Should be 1 state
        self.assertEqual(len(s.state_groups), 1)
        state1, = s.state_groups[0]
        # Should be 1 model
        self.assertEqual(sum(len(x) for x in state1), 1)
        # Check # of spheres and atoms in each model
        m = state1[0][0]
        self.assertEqual(len(m._spheres), 1698)
        self.assertEqual(len(m._atoms), 0)
        # Should be 1 ensemble
        self.assertEqual([e.num_models for e in s.ensembles], [9999])
        # One set of crosslinks
        xl1, = s.restraints
        self.assertEqual(len(xl1.experimental_cross_links), 90)
        self.assertEqual(len(xl1.cross_links), 90)
        self.assertEqual(xl1.linker.auth_name, 'DSSO')


if __name__ == '__main__':
    unittest.main()
