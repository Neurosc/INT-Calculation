#!/bin/bash
# Normalize all subjects with deobliquing and 3mm template

new_storage=/BICNAS2/group-northoff/jkokino
base_dir=/home/jkokino/meditation_project/data/preprocessed/BIDS_nifti
mni_template=/home/jkokino/meditation_project/templates/mni_icbm152_3mm.nii.gz

# Loop through groups and subjects
for group in controls meditators; do
    for subject_dir in $base_dir/$group/sub-*/; do
        subject=$(basename $subject_dir)
        
        anat_dir=$subject_dir/anat
        func_dir=$subject_dir/func
        
        
        mkdir -p $new_storage/normalized/$group/$subject/func
        
        cd $anat_dir
        
        transform_matrix=${subject}_T1w_at.nii.Xaff12.1D
        
        # Process each functional run
        for bold_file in $func_dir/*_bold.nii.gz; do
            [[ $bold_file == *"_space-MNI"* ]] && continue
            
            filename=$(basename $bold_file .nii.gz)
            output_file=$new_storage/normalized/$group/$subject/func/${filename}_space-MNI.nii.gz
            
            # Skip if already done
            [ -f $output_file ] && continue
                        
            # Deoblique first
            temp_deoblique=/tmp/${subject}_${filename}_deoblique.nii.gz
            3dWarp -deoblique -prefix $temp_deoblique $bold_file > /dev/null 2>&1
            
            # Normalize
            3dAllineate \
                -1Dmatrix_apply $transform_matrix \
                -input $temp_deoblique \
                -master $mni_template \
                -prefix $output_file \
                -final wsinc5 \
                > /dev/null 2>&1
            
            3drefit -space MNI $output_file
            
            # Clean up
            rm -f $temp_deoblique
        done
    done
done

echo "All subjects normalized!"