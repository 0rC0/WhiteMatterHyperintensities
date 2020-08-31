"""Convert minc to nifti format.
    Credits: https://gist.github.com/ofgulban/46d5c51ea010611cbb53123bb5906ca9
"""

import os
import numpy as np
from nibabel import load, save, Nifti1Image
import sys

def minc2nii(in_file):
    minc = load(in_file)
    basename = minc.get_filename().split(os.extsep, 1)[0]

    affine = np.array([[0, 0, 1, 0],
                       [0, 1, 0, 0],
                       [1, 0, 0, 0],
                       [0, 0, 0, 1]])

    out = Nifti1Image(minc.get_fdata(), affine=affine)
    save(out, basename + '.nii.gz')

if __name__=="__main__":
    minc2nii(sys.argv[1])