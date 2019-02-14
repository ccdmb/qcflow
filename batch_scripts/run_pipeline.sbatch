#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --partition=workq
#SBATCH --time=1-00:00:00
#SBATCH --account=y95
#SBATCH --mail-type=ALL
#SBATCH --mail-user=darcy.a.jones@postgrad.curtin.edu.au
#SBATCH --export=NONE

module load nextflow/19.01.0.5050-bin

nextflow run -resume -profile pawsey_zeus ./main.nf \
  --max_cpus 28 \
  --fastq "data/*_{1,2}.fq.gz" \
  --krakendb "databases/krakendb" \
  --contaminants "databases/krakendb/library/{bacteria,archaea,protozoa,viral,UniVec_Core,human}/library.fna" \
  --references "references/*.fasta" \
  --synthetic_contaminants "{data/synthetic_contaminants.fasta,databases/krakendb/library/UniVec_Core/library.fna}"
