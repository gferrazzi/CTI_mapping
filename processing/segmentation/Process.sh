# requires IRTK and FSL

# Segments CSF, WM and GM according to:

#1. Zhang, Y., Brady, M., and Smith, S. (2001). Segmentation of brain MR images through a hidden
#   Markov random field model and the expectation-maximization algorithm. IEEE transactions
#   on medical imaging 20(1), 45-57.

# creates a mask for the CSF by thresholding conductivity 

PATHDATA_MPRAGE_NII=$MOTHER_PATH'../CTI/MPRAGE_T1SPACE.nii'
PATHDATA_CTI_NII=$MOTHER_PATH'../CTI/1000/CTItensor.nii'
PATHDTI_b0=$MOTHER_PATH'../CTI/DTI.nii'

cp $PATHDATA_MPRAGE_NII .
cp $PATHDATA_CTI_NII .

bet $PATHDATA_MPRAGE_NII MPRAGE_T1SPACE_brain.nii -f 0.3 -g -0.3 -R

region $PATHDATA_CTI_NII Cxx.nii -Rt1 0 -Rt2 1
threshold Cxx.nii csf.nii 20000
threshold MPRAGE_T1SPACE_brain.nii.gz mask.nii 1
fslmaths csf.nii -mul mask.nii csf.nii
rm Cxx.nii
gunzip *

fast MPRAGE_T1SPACE_brain.nii
gunzip *
