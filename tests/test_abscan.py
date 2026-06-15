import unittest
import os
import tempfile
import numpy as np
from abscan.utils import get_vina_box, get_ligand_coords, get_obabel_path
from abscan.scanner import AlanineScanner

class TestAbscanUtils(unittest.TestCase):
    def test_get_vina_box(self):
        # Create dummy coordinates
        coords = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [2.0, 2.0, 2.0]
        ])
        center, box_size = get_vina_box(coords, padding=10.0)
        
        # Center should be [1.0, 1.0, 1.0]
        self.assertEqual(center, [1.0, 1.0, 1.0])
        # Box size should be max - min + padding = 2 - 0 + 10 = 12
        self.assertEqual(box_size, [12.0, 12.0, 12.0])

    def test_get_ligand_coords_parsing(self):
        # Create a temporary PDB file with HETATM records
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
            tmp.write("HETATM    1  C1  ZMR A 466      -2.787  16.675  -5.274  1.00 52.80\n")
            tmp.write("HETATM    2  C2  ZMR A 466      -1.924  17.878  -5.440  1.00 51.41\n")
            tmp_path = tmp.name
            
        try:
            coords = get_ligand_coords(tmp_path)
            self.assertEqual(coords.shape, (2, 3))
            self.assertAlmostEqual(coords[0][0], -2.787)
            self.assertAlmostEqual(coords[1][2], -5.440)
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

    def test_obabel_path_exists(self):
        try:
            path = get_obabel_path()
            self.assertTrue(os.path.exists(path))
        except RuntimeError:
            # If obabel is not installed, it raises RuntimeError
            pass


class TestAlanineScannerInit(unittest.TestCase):
    def test_invalid_scoring_function(self):
        with self.assertRaises(ValueError):
            # Should raise ValueError for unsupported scoring function
            AlanineScanner(
                pdb_file="dummy.pdb",
                resno=466,
                scoring_function="unsupported_sf"
            )

if __name__ == '__main__':
    unittest.main()
