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
  --outdir results \
  --fastq "fastq/*_R{1,2}.fastq.gz" \
  --references "references/*.fasta" \
  --adapter1 data/truseq_fwd.fasta \
  --adapter2 data/truseq_rev.fasta \
  --scontaminants data/synth_cont.fasta \
  --filter_phred 5 \
  --trim_phred 2 \
  --min_read_length 50 
 
#  --map
#  --krakendb "databases/krakendb" \
#  --kraken_low_mem
