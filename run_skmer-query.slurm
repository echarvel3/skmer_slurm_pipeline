#!/bin/bash
#SBATCH --job-name="skimming_pipeline"
#SBATCH --partition=shared
#SBATCH --account=uot138
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --array=1-24%50
#SBATCH --mem=30G
#SBATCH -o logs/subsample-and-estimate_%j_%a.out
#SBATCH -e logs/subsample-and-estimate_%j_%a.err
#SBATCH -t 10:00:00

module load cpu/0.15.4
module load parallel/20200822

source ${HOME}/.bashrc
source ${HOME}/miniconda3/etc/profile.d/conda.sh

echo "========"
set -x

date
conda activate skimming_scripts

REFERENCE_GENOME="/expanse/lustre/scratch/echarvel3/temp_project/skmer2_realdata/haploid-bee_data/ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.hist"
REFERENCE_GENOME="/expanse/lustre/scratch/echarvel3/temp_project/skmer2_realdata/haploid-drosophila_data/only_chroms_noY.hist"

if [ ! -d "./dummy_lib" ]; then
	mkdir ./dummy_lib
	echo -e "kmer_length\t31\nsketch_size\t100000\nsketching_seed\t42" > ./dummy_lib/CONFIG
fi

species=komodo

for OUTPUT in 1x_downsampled 2x_downsampled 4x_downsampled 6x_downsampled; do
	#REFERENCE_GENOME="../../respect-hists/4x_downsampled_output/skmer-reference.txt"
	REFERENCE_RESPECT=$(realpath "./respect-hists/${OUTPUT}_output/skmer-reference.txt")

	INPUT="./downsampled_fastq/folder_${SLURM_ARRAY_TASK_ID}/${OUTPUT}/"

	for x in $(realpath ./${INPUT}/*); do
		# RUNNING SKMER1
		mkdir --parents "./${species}_skmer/${OUTPUT}/"
		pushd "./${species}_skmer/${OUTPUT}/"
		skmer query $x "../../dummy_lib/" -a -p 2
		popd

		# RUNNING SKMER REPEAT + ASSEMBLY
		#mkdir --parents "./${species}_skmer-repeat-assembly/${OUTPUT}/"
		#pushd "./${species}_skmer-repeat-assembly/${OUTPUT}/"
		#python /expanse/lustre/projects/uot138/echarvel3/cunner_reads/scripts/Skmer-2/skmer/new_dist_newSmoothSizeCorr.py query $x "../../dummy_lib/" -a -p 2 -r ${REFERENCE_GENOME}
		#popd

		# RUNNING SKMER REPEAT + RESPECT
		# mkdir --parents "./${species}_skmer-repeat-respect/${OUTPUT}/"
		# pushd "./${species}_skmer-repeat-respect/${OUTPUT}/"
		# python /expanse/lustre/projects/uot138/echarvel3/cunner_reads/scripts/Skmer-2/skmer/new_dist_newSmoothSizeCorr.py query $x "../../dummy_lib/" -a -p 2 -r ${REFERENCE_RESPECT}
		# popd

		# RUNNING SKMER DIPLOID
		mkdir --parents "./${species}_skmer-diploid/${OUTPUT}/"
		pushd "./${species}_skmer-diploid/${OUTPUT}/"
		python /expanse/lustre/projects/uot138/echarvel3/cunner_reads/scripts/Skmer-2/skmer/__diploid_skmer_n__.py query $x "../../dummy_lib/" -a -p 2
		popd
	done
done

