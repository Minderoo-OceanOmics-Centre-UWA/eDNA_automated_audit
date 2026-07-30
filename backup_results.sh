project_dir=$1

project=$(basename "$(dirname "$project_dir")")

mkdir ${project_dir}/fastqs
mv ${project_dir}/*/*/01-cutadapt/assigned/* ${project_dir}/fastqs/
rm -r ${project_dir}/*/*/01-cutadapt/assigned/

rclone copy ${project_dir} s3:ocom-edna/illumina-meta-analysed-data/${project}
rclone check ${project_dir} s3:ocom-edna/illumina-meta-analysed-data/${project}
