#!/bin/bash
#SBATCH --job-name="consult-2"
#SBATCH --partition=shared
#SBATCH --account=uot138
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --array=1-24
#SBATCH --mem=20G
#SBATCH -o logs/consult-2_%j_%a.out
#SBATCH -e logs/consult-2_%j_%a.err
#SBATCH -t 48:00:00
echo "========"
RANDOM_SEED=0
NTHREADS=${NTHREADS:-10}
OUTPUT="./nodecontam_skmer/folder_${SLURM_ARRAY_TASK_ID}"

date
source ${HOME}/.bashrc
source ${HOME}/miniconda3/etc/profile.d/conda.sh

conda activate skimming_scripts
#conda activate skmer2

echo "========"
#bash /expanse/lustre/scratch/echarvel3/temp_project/skmer2_realdata/skimming_scripts-echarvel/run_skmer.sh -i ./bbmap_reads/folder_${SLURM_ARRAY_TASK_ID} \
#	-o ${OUTPUT} \
#	-t ${NTHREADS} \
#	-p 2 && echo "The batch ${SLURM_ARRAY_TASK_ID} has been successfully processed."

SCRIPT_DIR="../skimming_scripts-echarvel"
out_dir="./"

mkdir --parent "${out_dir}/consult_reads/folder_${SLURM_ARRAY_TASK_ID}"
mkdir --parent "${out_dir}/consult_reads/folder_class-info_${SLURM_ARRAY_TASK_ID}"

${SCRIPT_DIR}/CONSULT-II/consult_search \
	-i "${SCRIPT_DIR}/CONSULT-II/library-h13_t2_l2_b16/" \
	-o "${out_dir}/consult_reads/folder_${SLURM_ARRAY_TASK_ID}" \
	--query-path "${out_dir}/bbmap_reads/folder_${SLURM_ARRAY_TASK_ID}" \
	--number-of-matches 2 \
	--thread-count ${NTHREADS} \
	--unclassified-out "${out_dir}/consult_reads/folder_${SLURM_ARRAY_TASK_ID}"

echo "========"
date
