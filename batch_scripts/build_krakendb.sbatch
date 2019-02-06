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

module load kraken2/2.0.7-beta-gcc
module load blast+/2.7.1

export DBNAME=${DBNAME-krakendb}


SCRIPT_FILE=".download_kraken-$(date +%s).sh"

cat << EOF > ${SCRIPT_FILE}
#!/usr/bin/env bash

set -eu

mkdir -p ${DBNAME}

kraken2-build --download-taxonomy --db ${DBNAME}

kraken2-build --threads 28 --download-library bacteria --db ${DBNAME}
kraken2-build --threads 28 --download-library archaea --db ${DBNAME}
kraken2-build --threads 28 --download-library viral --db ${DBNAME}
kraken2-build --threads 28 --download-library UniVec_Core --db ${DBNAME}
kraken2-build --threads 28 --download-library fungi --db ${DBNAME}

kraken2-build \
    --threads 28 \
    --build \
    --db ${DBNAME} \
    --kmer-len 35 \
    --minimizer-len 31 \
    --minimizer-spaces 6

kraken2-build --clean --db ${DBNAME}
EOF

chmod a+x ${SCRIPT_FILE}

srun --export=all -n 1 -N 1 -c 28 ${SCRIPT_FILE}
#rm -f ${SCRIPT_FILE}
