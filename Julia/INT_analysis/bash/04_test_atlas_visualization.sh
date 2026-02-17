#!/bin/bash
# Test atlas alignment on one subject

master_atlas=/home/jkokino/meditation_project/templates/glasser_nonself_atlas_4mm_spheres.nii.gz
test_bold=/home/jkokino/meditation_project/data/preprocessed/BIDS_nifti/meditators/sub-032/func/sub-032_task-ca2_bold.nii.gz
output_dir=/home/jkokino/meditation_project/data/preprocessed/BIDS_nifti/meditators/sub-032/func

# Resample master atlas to this subject's functional space
3dresample \
    -input $master_atlas \
    -master $test_bold \
    -prefix $output_dir/sub-032_task-ca2_atlas_TEST.nii.gz \
    -overwrite

