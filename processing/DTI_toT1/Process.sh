#!/bin/bash

# requires IRTK and FSL

# Implements rigid registration of DTI pre-processed data onto \sigma_HF space according to:

#1. Studholme, C., Hill, D.L.G., and Hawkes, D.J. (1999). An overlap invariant entropy measure of 3D
#   medical image alignment. Pattern recognition 32(1), 71-86.

# Data paths
PATHDATA_T1_NII='../../data/nifties/T1/TR700.nii'
PATHDATA_DTI_NII='../DTI_Preprocessing/DTI_denoised_bias_dist.nii'
PATHDATA_DTI_NII_out='../CTI/DTI.nii'
PATHBVEC='../../data/nifties/DTI/DTI.bvec'
PATHBVAL='../../data/nifties/DTI/DTI.bval'
PATHBVEC_out='../CTI/DTI.bvec'
PATHBVAL_out='../CTI/DTI.bval'

rm DTItoT1.dof
rm b0.nii

bet $PATHDATA_T1_NII T1_skullstripped.nii -R
gunzip *

fslroi $PATHDATA_DTI_NII b0.nii 0 1
gunzip *

rreg b0.nii T1_skullstripped.nii -dofout DTItoT1.dof
dofinvert DTItoT1.dof DTItoT1.dof

convert $PATHDATA_DTI_NII DTIdiv.nii -double
fslmaths DTIdiv.nii -div 1 DTIdiv.nii
gunzip *
transformation DTIdiv.nii $PATHDATA_DTI_NII_out -bspline -target $PATHDATA_T1_NII -dofin DTItoT1.dof

cp $PATHBVEC $PATHBVEC_out
cp $PATHBVAL $PATHBVAL_out
gunzip *
