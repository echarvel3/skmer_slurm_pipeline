#!/bin/bash
#SBATCH --job-name="downsample"
#SBATCH --partition=shared
#SBATCH --account=uot138
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=10G
#SBATCH -o logs/respect_full_%j_%a.out
#SBATCH -e logs/respect_full_%j_%a.err
#SBATCH -t 00:40:00
#SBATCH --array 1-24

module load cpu/0.15.4
module load parallel/20200822

source ${HOME}/.bashrc
source ${HOME}/miniconda3/etc/profile.d/conda.sh

conda activate skimming_scripts

echo "========"
set -x
NTHREADS=${NTHREADS:-20}

FASTQ_DIR="../decontam_temp/folder_${SLURM_ARRAY_TASK_ID}/"
FASTQ_DIR="./kraken_reads/folder_${SLURM_ARRAY_TASK_ID}/"
JELLY_DIR="./decontam_jellyfish/folder_${SLURM_ARRAY_TASK_ID}/"
RESPECT_DIR="./respect-hists/6x_downsampled_output/"
PARAMETER_FILE="./full_cov_respect_output/estimated-parameters.txt"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
SCRIPT_DIR="/expanse/lustre/projects/uot138/echarvel3/cunner_reads/scripts/genome_size_scripts/scripts/"
OUT_DIR="./downsampled_fastq/folder_${SLURM_ARRAY_TASK_ID}/"

date

mkdir --parents "${OUT_DIR}/6x_downsampled"
mkdir --parents "${OUT_DIR}/4x_downsampled"
mkdir --parents "${OUT_DIR}/2x_downsampled"
mkdir --parents "${OUT_DIR}/1x_downsampled"

echo "========"

#Uses respect coverage estimates and read count estimates to downsample.


for x in $(realpath ${FASTQ_DIR}/unclass*); do
	filename=${x##*/}
	#line_count=$(grep "${filename}" "${JELLY_DIR}/read_counts.txt" | cut -f2)
	#coverage=$(grep "${filename%.fastq}" ${PARAMETER_FILE} | cut -f4)
	
	# FOR DOWNSAMPLING WITH RESPECT:
	#line_count=$(grep "${filename%.fastq}" "${RESPECT_DIR}/estimated-parameters.txt" | cut -f9)
	#genome_size=$(grep "${filename%.fastq}" "${RESPECT_DIR}/estimated-parameters.txt" | cut -f5)
	for y in 1 2 4 6; do
		target_cov=$y

		#sample_num=$(printf %.0f $(echo "${line_count}/${coverage}*${target_cov}" | bc -l))

		#sample_num=$((728370775*"${y}"/150)) #CUNNER
		#sample_num=$((225000000*"${y}"/150)) #BEE
		#sample_num=$((143000000*"${y}"/150)) #DROSOPHILA
		sample_num=$((1508000000*"${y}"/150)) #KOMODO
		#sample_num=$((${genome_size}*"${y}"/${line_count})) #RESPECT BASED
	
		${SCRIPT_DIR}/other_scripts/bbmap/reformat.sh in="${x}" out="${OUT_DIR}/${y}x_downsampled/${filename%.fastq}.fastq" samplereadstarget="${sample_num}" sampleseed="${y}" overwrite=true
		#/usr/bin/time -v seqtk sample -s${y} -2 "${x}" "${sample_num}" > "${OUT_DIR}/${filename%.fastq}_rep${y}.fastq"
	done
done

echo "========"
date
