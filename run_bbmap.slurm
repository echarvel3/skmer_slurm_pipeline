#!/bin/bash
#SBATCH --job-name="bbmap_pipeline"
#SBATCH --partition=shared
#SBATCH --account=uot138
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --array=1-24
#SBATCH -o logs/bbmap_%j-%a.out
#SBATCH -e logs/bbmap_%j-%a.err
#SBATCH -t 6:00:00

source ${HOME}/.bashrc
source ${HOME}/programs/miniconda3/etc/profile.d/conda.sh
set -x

echo "========"
RANDOM_SEED=0
NTHREADS=${NTHREADS:-15}
INPUT="/expanse/lustre/projects/uot138/echarvel3/cunner_reads/raw/folder_${SLURM_ARRAY_TASK_ID}/"
INPUT="/expanse/lustre/scratch/echarvel3/temp_project/skmer2_realdata/komodo_data/komodo_download/folder_${SLURM_ARRAY_TASK_ID}/"
READ_1="_1.fastq"
READ_2="_2.fastq"

date
echo "========"
mkdir -p ./bbmap_reads/folder_${SLURM_ARRAY_TASK_ID}/


for file in $(realpath "${INPUT}/*${READ_1}"); do
        read1=$file
        read2=${read1%${READ_1}}${READ_2}
        out_read=${file##*/}

        if [ ! -f "./bbmap_reads/folder_${SLURM_ARRAY_TASK_ID}/${out_read%${READ_1}}.fastq" ]; then
	    /expanse/lustre/projects/uot138/echarvel3/cunner_reads/scripts/genome_size_scripts/scripts/other_scripts/bbmap_pipeline.sh "$read1" "$read2" "./bbmap_reads/folder_${SLURM_ARRAY_TASK_ID}/${out_read%${READ_1}}.fastq"
        fi
done


conda activate skimming_scripts
echo "========"
date
