#!/bin/sh
#PBS -l nodes=node11:ppn=16
#PBS -q fat

samples=(P26770_1001 P26770_1002 P26770_1003)

# Base directories
base_dir="/home/users/yliu/projects/scRNA/Kenny_v2/cellranger/count"
cellranger="/home/users/yliu/software/cellranger-7.2.0/cellranger"
transcriptome="/home/users/yliu/projects/scRNA/Kenny/ref/refdata-gex-GRCm39-2024-A"
fastq_base="/home/users/yliu/projects/scRNA/Kenny/P26770"

for sample in "${samples[@]}"; do
    cd "${base_dir}/${sample}"    
    $cellranger count \
        --id="run_count_${sample}" \
        --transcriptome="$transcriptome" \
        --fastqs="${fastq_base}/${sample}/02-FASTQ/220726_A01901_0013_AHC2YGDRX2" \
        --sample="$sample" \
        --localcores=16
done

# change samples
samples=(P32113_1001 P32113_1002 P32113_1003 P32113_1004 P32113_1005 P32113_1006 P32113_1007 P32113_1008)

# change fastq file directory
fastq_base="/home/users/yliu/data/scRNA/Keying/DataDelivery_2024-07-08_23-18-33_ngisthlm00923/files/P32113"

for sample in "${samples[@]}"; do
    cd "${base_dir}/${sample}"
    $cellranger count \
        --id="run_count_${sample}" \
        --transcriptome="$transcriptome" \
        --fastqs="${fastq_base}/${sample}/02-FASTQ/20240705_LH00202_0123_B227W3VLT4" \
        --sample="$sample" \
        --localcores=16
done
