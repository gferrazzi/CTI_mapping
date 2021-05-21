#!/bin/bash

# requires IRTK and FSL

# Implements rigid registration of MP-RAGE data onto \sigma_HF space according to:

#1. Studholme, C., Hill, D.L.G., and Hawkes, D.J. (1999). An overlap invariant entropy measure of 3D
#   medical image alignment. Pattern recognition 32(1), 71-86.

# Data paths
PATHDATA_MPRAGE_NII="$MOTHER_PATH"'../../data/nifties/MPRAGE/MPRAGE.nii.gz'
PATHDATA_MPRAGE_NII_out="$MOTHER_PATH"'../CTI/MPRAGE_T1SPACE.nii'
PATHDATA_T1_NII="$MOTHER_PATH"'../../data/nifties/T1/TR700.nii'
PATHDATA_MICHEL_NII="$MOTHER_PATH"'../Michel_et_al/conductivity.nii'
PATHDATA_MICHEL_NII_out="$MOTHER_PATH"'../CTI/conductivity_hf.nii'
PATHDATA_MPRAGE_BET_NII_out="$MOTHER_PATH"'../CTI/MPRAGE_T1SPACE_BET_MASK.nii'

rm MPRAGEtoT1.dof
rm MPRAGE_skull.nii

bet $PATHDATA_T1_NII T1_skullstripped.nii -R
gunzip *

bet $PATHDATA_MPRAGE_NII MPRAGE_skullstripped.nii -R
gunzip *

rreg MPRAGE_skullstripped.nii T1_skullstripped.nii -dofout MPRAGEtoT1.dof
dofinvert MPRAGEtoT1.dof MPRAGEtoT1.dof

transformation $PATHDATA_MPRAGE_NII $PATHDATA_MPRAGE_NII_out -bspline -target $PATHDATA_T1_NII -dofin MPRAGEtoT1.dof

cp $PATHDATA_MICHEL_NII $PATHDATA_MICHEL_NII_out

bet $PATHDATA_MPRAGE_NII_out MPRAGE_T1SPACE_BET.nii -R -m
gunzip *
rm MPRAGE_T1SPACE_BET.nii
mv MPRAGE_T1SPACE_BET_mask.nii $PATHDATA_MPRAGE_BET_NII_out
