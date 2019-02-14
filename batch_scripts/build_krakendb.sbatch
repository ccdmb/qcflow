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

export NCPU=28
export DBNAME=${DBNAME-krakendb}
export DBS=${DBS-"bacteria archaea protozoa viral UniVec_Core fungi human"}

# You can add your own genomes to REFS.
# An example of how to add taxids required by kraken using sed is below.
# Just replace 321614 with the actual ncbi taxid
# sed -r '/^>/{s/(^>[^[:space:]]*)/\1|kraken:taxid|321614/g}' \
#   SN15v9_OM_Chr_and_tigs.fasta \
#   > PnodSn15.fasta
export REFS=${REFS-""}

SCRIPT_FILE="download_kraken-$(date +%s).sh"

cat << EOF > ${SCRIPT_FILE}
#!/usr/bin/env bash

set -eu

mkdir -p ${DBNAME}

kraken2-build --download-taxonomy --db ${DBNAME}

for db in ${DBS}; do
  kraken2-build \
    --threads ${NCPU} \
    --download-library \${db} \
    --db ${DBNAME}
done

for fasta in ${REFS}; do
    kraken2-build \
      --threads ${NCPU} \
      --no-masking \
      --add-to-library \${fasta} \
      --db ${DBNAME}
done

kraken2-build \
    --threads ${NCPU} \
    --build \
    --db ${DBNAME} \
    --kmer-len 35 \
    --minimizer-len 31 \
    --minimizer-spaces 6

# kraken2-build --clean --db ${DBNAME}
EOF

chmod a+x ${SCRIPT_FILE}

srun --export=all -n 1 -N 1 -c ${NCPU} ${SCRIPT_FILE}
#rm -f ${SCRIPT_FILE}
