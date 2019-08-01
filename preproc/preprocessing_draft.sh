#!/usr/bin/env bash

# Usage: preprocessing_draft.sh [Subject ID]

sub_id=$1
sub_dir='../../20190801'

# Coregister T1 to MNI
fslreorient2std $sub_dir/$sub_id/anat/$sub_id.anat.nii.gz $sub_dir/$sub_id/anat/$sub_id.anat_ro.nii.gz
antsRegistrationSyNQuick.sh -d 3 -n 8 -m $sub_dir/$sub_id/anat/$sub_id.anat_ro.nii.gz -f /usr/share/fsl/data/standard/MNI152_T1_1mm.nii.gz -o $sub_dir/$sub_id/anat/$sub_id.anat_t1_2mni

# two steps coregistration FLAIR to T1 and then to MNI
fslreorient2std $sub_dir/$sub_id/flair/$sub_id.flair.nii.gz $sub_dir/$sub_id/flair/$sub_id.flair_ro.nii.gz
antsRegistrationSyNQuick.sh -d3 -n8 -m $sub_dir/$sub_id/flair/$sub_id.flair_ro.nii.gz -f $sub_dir/$sub_id/anat/$sub_id.anat_ro.nii.gz -o $sub_dir/$sub_id/flair/$sub_id.flair_to_t1
antsApplyTransforms -d 3 -i $sub_dir/$sub_id/flair/$sub_id.flair_to_t1Warped.nii.gz -r /usr/share/fsl/data/standard/MNI152_T1_1mm.nii.gz -t $sub_dir/$sub_id/anat/$sub_id.anat_t1_2mni0GenericAffine.mat -t $sub_dir/$sub_id/anat/$sub_id.anat_t1_2mni1Warp.nii.gz -o $sub_dir/$sub_id/flair/$sub_id.flair_2t1_2mni_coreg.nii.gz

sub_id = "${sub_id}_T2"
# Coregister T2 to MNI
fslreorient2std $sub_dir/$sub_id/anat/$sub_id.anat.nii.gz $sub_dir/$sub_id/anat/$sub_id.anat_ro.nii.gz
antsRegistrationSyNQuick.sh -d 3 -n 8 -m $sub_dir/$sub_id/anat/$sub_id.anat_ro.nii.gz -f /usr/share/fsl/data/standard/MNI152_T1_1mm.nii.gz -o $sub_dir/$sub_id/anat/$sub_id.anat_t2_2mni

# two steps coregistration T2 FLAIR to T2 and then to MNI
fslreorient2std $sub_dir/$sub_id/flair/$sub_id.flair.nii.gz $sub_dir/$sub_id/flair/$sub_id.flair_ro.nii.gz
antsRegistrationSyNQuick.sh -d3 -n8 -m $sub_dir/$sub_id/flair/$sub_id.flair_ro.nii.gz -f $sub_dir/$sub_id/anat/$sub_id.anat_ro.nii.gz -o $sub_dir/$sub_id/flair/$sub_id.flair_to_t2
antsApplyTransforms -d 3 -i $sub_dir/$sub_id/flair/$sub_id.flair_to_t2Warped.nii.gz -r /usr/share/fsl/data/standard/MNI152_T1_1mm.nii.gz -t $sub_dir/$sub_id/anat/$sub_id.anat_t2_2mni0GenericAffine.mat -t $sub_dir/$sub_id/anat/$sub_id.anat_t2_2mni1Warp.nii.gz -o $sub_dir/$sub_id/flair/$sub_id.flair_2t2_2mni_coreg.nii.gz

#N4 Correction for T1 and FLAIRs

N4BiasFieldCorrection -d 3 -i $sub_dir/$sub_id/anat/$sub_id.anat_t1_2mniWarped.nii.gz -o $sub_dir/$sub_id/anat/$sub_id.anat_t1_2mniWarped_n4.nii.gz

N4BiasFieldCorrection -d 3 -i $sub_dir/$sub_id/flair/$sub_id.flair_2t1_2mni_coreg.nii.gz \
                      -o $sub_dir/$sub_id/flair/$sub_id.flair_2t1_2mni_coreg_n4.nii.gz

N4BiasFieldCorrection -d 3 -i $sub_dir/$sub_id/flair/$sub_id.flair_2t2_2mni_coreg.nii.gz \
                      -o $sub_dir/$sub_id/flair/$sub_id.flair_2t2_2mni_coreg_n4.nii.gz

#Brain Extraction
bet $sub_dir/$sub_id/anat/$sub_id.anat_t1_2mniWarped_n4.nii.gz $sub_dir/$sub_id/anat/$sub_id.anat_t1_2mniWarped_n4_brain.nii.gz
fslmaths $sub_dir/$sub_id/anat/$sub_id.anat_t1_2mniWarped_n4_brain.nii.gz -bin $sub_dir/$sub_id/anat/$sub_id.anat_t1_2mniWarped_n4_brain_mask.nii.gz

#Apply Brain Mask to FLAIRs
fslmaths $sub_dir/$sub_id/flair/$sub_id.flair_2t1_2mni_coreg_n4.nii.gz -mul $sub_dir/$sub_id/anat/$sub_id.anat_t1_2mniWarped_n4_brain_mask.nii.gz
fslmaths $sub_dir/$sub_id/flair/$sub_id.flair_2t2_2mni_coreg_n4.nii.gz -mul $sub_dir/$sub_id/anat/$sub_id.anat_t1_2mniWarped_n4_brain_mask.nii.gz
