"""
The weights included in the pipeline used a version of joblib<0.21. For versions >0.23 sklearn.externals.joblib is
deprecated and removed, so the pkl model can not be loaded anymore.

This script is to update the pkl files to use the standalone joblib.

References: https://stackoverflow.com/questions/61893719/importerror-cannot-import-name-joblib-from-sklearn-externals

Usage:
    conda create -n update_pkl scikit-learn=0.21 joblib
    conda activate update_pkl
    python update_pkl.py [pklfile] [output_pkl]
"""

import sklearn.externals.joblib as extjoblib
import joblib
from sys import argv


def update_pkl(in_pkl, out_pkl):
    old_pkl = extjoblib.load(in_pkl)
    joblib.dump(old_pkl, out_pkl)


if __name__=="__main__":
    if not argv[1] or not argv[2]:
        raise FileNotFoundError("Usage:python update_pkl.py [pklfile] [output_pkl]")
    update_pkl(argv[1], argv[2])