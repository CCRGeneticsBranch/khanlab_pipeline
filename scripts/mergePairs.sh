#!/bin/bash
#SBATCH --partition=ccr,norm
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=3-00:00:00
#SBATCH --gres=lscratch:200

module load bwa
module load samtools
module load samblaster
module load preseq
module load juicer

#example:
#sbatch mergePairs.sh 'out_yong2/RH4_DMSO_24h_*/*PT.pairs.gz' RH4_DMSO_24h_merged out_yong2
#sbatch mergePairs.sh 'out_yong2/RH4_PFI90_24h_*/*PT.pairs.gz' RH4_PFI90_24h_merged out_yong2
in=$1
prefix=$2
out_dir=$3

source /data/khanlab/projects/hsienchao/conda/etc/profile.d/conda.sh
conda activate dovetail

mkdir -p $out_dir
mkdir -p $out_dir/${prefix}

ref=/data/khanlab/projects/hsienchao/ref/bwa/ucsc.hg19.fasta

pairtools merge -o $out_dir/${prefix}/${prefix}.pairs.gz --tmpdir /lscratch/$SLURM_JOB_ID --memory 120G --nproc $SLURM_CPUS_PER_TASK --nproc-in $SLURM_CPUS_PER_TASK --nproc-out $SLURM_CPUS_PER_TASK $in

./contact_map.sh $out_dir/${prefix}/${prefix}.pairs.gz $ref $out_dir/${prefix}/${prefix} $SLURM_CPUS_PER_TASK


