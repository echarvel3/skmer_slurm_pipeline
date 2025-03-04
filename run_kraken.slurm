#!/bin/bash
#SBATCH --job-name="kraken-2"
#SBATCH --partition=shared
#SBATCH --account=uot138
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --array=1-24
#SBATCH --mem=20G
#SBATCH -o logs/kraken-2_%j_%a.out
#SBATCH -e logs/kraken-2_%j_%a.err
#SBATCH -t 12:00:00
echo "========"
RANDOM_SEED=0
NTHREADS=${NTHREADS:-10}
OUTPUT="./nodecontam_skmer/folder_${SLURM_ARRAY_TASK_ID}"

date
source ${HOME}/.bashrc
source ${HOME}/miniconda3/etc/profile.d/conda.sh

conda activate skimming_scripts
#conda activate skmer2
conda activate kraken2

echo "========"
#bash /expanse/lustre/scratch/echarvel3/temp_project/skmer2_realdata/skimming_scripts-echarvel/run_skmer.sh -i ./bbmap_reads/folder_${SLURM_ARRAY_TASK_ID} \
#	-o ${OUTPUT} \
#	-t ${NTHREADS} \
#	-p 2 && echo "The batch ${SLURM_ARRAY_TASK_ID} has been successfully processed."

SCRIPT_DIR="./skimming_scripts-echarvel"
out_dir="./"

mkdir --parents ${out_dir}/kraken_reads/folder_${SLURM_ARRAY_TASK_ID}

for file in $(ls "${out_dir}/consult_reads/folder_${SLURM_ARRAY_TASK_ID}/"); do
	    kraken2 \
		--db ../krakenlib/ \
		"${out_dir}/consult_reads/folder_${SLURM_ARRAY_TASK_ID}/${file}" \
		--threads ${NTHREADS} \
		--unclassified-out ${out_dir}/kraken_reads/folder_${SLURM_ARRAY_TASK_ID}/${file%".fastq"}.fastq
done


echo "========"
date
