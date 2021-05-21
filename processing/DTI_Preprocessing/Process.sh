#!/bin/bash

# requires FSL and MRtrix3, 

# Implements denoising of the DTI data according to:

#1. Veraart, J., Fieremans, E., and Novikov, D.S. (2016a). Diffusion MRI noise mapping using random
#   matrix theory. Magnetic resonance in medicine 76(5), 1582-1593
#2. Veraart, J., Novikov, D.S., Christiaens, D., Ades-Aron, B., Sijbers, J., and Fieremans, E. (2016b).
#   Denoising of diffusion MRI using random matrix theory. Neuroimage 142, 394-406

# Implements bias field corrections according to

#3. Tournier, J.D., Smith, R., Raffelt, D., Tabbara, R., Dhollander, T., Pietsch, M., et al. (2019).
#   MRtrix3: A fast, flexible and open software framework for medical image processing and
#   visualisation. NeuroImage 202, 116137

# Implements distortion correction according to

#4. Andersson, J.L.R., Skare, S., and Ashburner, J. (2003). How to correct susceptibility distortions in
#   spin-echo echo-planar images: application to diffusion tensor imaging. Neuroimage 20(2),
#   870-888

# Data paths
PATHDATA_NII='../../data/nifties/DTI/DTI.nii.gz'
PATHDENOISED_NII='DTI_denoised.nii'
PATHDENOISED_BIASFIELD_NII='DTI_denoised_bias.nii'
PATHDENOISED_BIASFIELD_DIST_NII='DTI_denoised_bias_dist.nii'
PATHDATA_AP_NII="$MOTHER_PATH"'../../data/nifties/DTI/DTI-AP-A.nii.gz'
PATHDATA_PA_NII='../../data/nifties/DTI/DTI-AP-P.nii.gz'
BVAL='../../data/nifties/DTI/DTI.bval'
BVEC='../../data/nifties/DTI/DTI.bvec'

# Denoising + bias field correction
dwidenoise $PATHDATA_NII $PATHDENOISED_NII -force
dwibiascorrect fsl $PATHDENOISED_NII $PATHDENOISED_BIASFIELD_NII -fslgrad $BVEC $BVAL -bias biasfield.nii -force
 
# Distortion correction
mrcat $PATHDATA_AP_NII $PATHDATA_PA_NII b0s.nii -axis 3 -force
mkdir dwi/
cd dwi/
topup --imain=../b0s.nii --datain=../acqparams.txt --config=b02b0.cnf --out=topUpped --fout=my_field.nii --iout=my_unwarped_images.nii -v
applytopup --imain=../DTI_denoised_bias.nii --datain=../acqparams.txt --inindex=2 --topup=topUpped --out=../DTI_denoised_bias_dist.nii --method=jac
cd ../
gunzip *
