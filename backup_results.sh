project_dir=$1

project=$(basename $(dirname ${project_dir}))
mkdir -p tmp/tmp_fastqs_${project}

# Remove project prefix from demuxed untrimmed fqs and change location
for file in ${project_dir}/*/*/01-cutadapt/assigned/*; do
    filename=$(basename ${file})
    newname=${filename#${project}_}
    mv ${file} tmp/tmp_fastqs_${project}/${newname}   
done

# Rename trimmed fqs and change location
for file in ${project_dir}/*/*/01-cutadapt/all-primers-trimmed/*; do
    filename=$(basename ${file})
    newname=${filename#trimmed_trimmed_${project}_}
    assay=${newname#*_}
    assay=${assay%%_*}
    
    mkdir -p ${project_dir}/fastqs/${assay}
    mv ${file} ${project_dir}/fastqs/${assay}/${newname}
done

# Backup project directory
rclone copy ${project_dir} s3:ocom-edna/illumina-meta-analysed-data/${project}
rclone check ${project_dir} s3:ocom-edna/illumina-meta-analysed-data/${project}

# Backup demuxed untrimmed fqs to different location
rclone copy tmp/tmp_fastqs_${project}/* s3:ocom-edna/illumina-meta-sra/${project}
rclone check tmp/tmp_fastqs_${project}/* s3:ocom-edna/illumina-meta-sra/${project}

# Cleanup
rm -r tmp/tmp_fastqs_${project}
