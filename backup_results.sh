project_dir=$1

project=$(basename $(dirname $project_dir))

for file in ${project_dir}/*/*/01-cutadapt/assigned/*; do
    filename=$(basename "$file")
    newname=${filename#${project}_}

    mv $file" "$(dirname "$file")/$newname
done

rclone copy ${project_dir} s3:ocom-edna/illumina-meta-analysed-data/${project}
rclone check ${project_dir} s3:ocom-edna/illumina-meta-analysed-data/${project}

rclone copy ${project_dir}/*/*/01-cutadapt/assigned/* s3:ocom-edna/illumina-meta-sra/${project}
rclone check ${project_dir}/*/*/01-cutadapt/assigned/* s3:ocom-edna/illumina-meta-sra/${project}
