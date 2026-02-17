#!/bin/bash
# Create master Glasser atlas using 3dUndump

template_1mm=/home/jkokino/meditation_project/templates/MNI_files/mni_icbm152_nlin_asym_09c.nii
template_3mm=/home/jkokino/meditation_project/templates/MNI_files/mni_icbm152_3mm.nii.gz
coords_file=/home/jkokino/meditation_project/templates/glasser_files/glasser_coordinates_nonself.txt
output_dir=/home/jkokino/meditation_project/templates/glasser/files

cd $output_dir

# Convert coordinates file to 3dUndump format (X Y Z Value)
tail -n +2 $coords_file | awk '{print $3, $4, $5, $1}' > glasser_coords_undump.txt

# Create all spheres at once using 3dUndump
3dUndump \
    -prefix glasser_atlas_1mm \
    -master $template_1mm \
    -srad 4 \
    -xyz glasser_coords_undump.txt

# Convert to NIfTI
3dAFNItoNIFTI -prefix glasser_atlas_1mm.nii.gz glasser_atlas_1mm+tlrc

# Downsample to 3mm
3dresample \
    -input glasser_atlas_1mm.nii.gz \
    -master $template_3mm \
    -prefix glasser_nonself_atlas_4mm.nii.gz \
    -rmode NN

# Clean up intermediate files
rm -f glasser_atlas_1mm+tlrc.* glasser_coords_undump.txt

echo " Master atlas created: glasser_nonself_atlas_4mm.nii.gz"
