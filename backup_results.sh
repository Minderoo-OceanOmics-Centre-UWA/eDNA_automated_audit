project=$1

rclone copy ${project} s3:ocom-edna/illumina-meta-analysed-data/${project}
rclone check ${project} s3:ocom-edna/illumina-meta-analysed-data/${project}

rclone copy ${project}/*/*/01-cutadapt/assigned* s3:ocom-edna/illumina-meta-sra/${project}/
rclone check ${project}/*/*/01-cutadapt/assigned* s3:ocom-edna/illumina-meta-sra/${project}/
