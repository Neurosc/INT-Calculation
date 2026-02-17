#!/bin/bash


cd /BICNAS2/group-northoff/jkokino/normalized

find . -name "*_task-conc1_space-MNI.nii.gz" | while read file; do
    echo "Testing: $file"
    3dinfo -n4 "$file" > /dev/null 2>&1
    if [ $? -ne 0 ]; then
        echo "  → CORRUPTED: $file"
    fi
done