#!/bin/bash
#SBATCH --job-name=gtdbtk_de_novo_wf
#SBATCH -p himem
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cat3@hi.is
#SBATCH -N 1
#SBATCH --ntasks-per-node=1

echo $HOSTNAME

source /users/home/cat3/.bashrc

conda activate gtdbtk

WD=/users/work/cat3/projects/mime/results/de_novo_wf
GENOMES_DIR=$WD/surface-arch-mags/
OUTDIR=$WD/surface-arch-tree/

GTDBTK_DATA_PATH=/users/work/cat3/db/gtdbk

cd $WD

gtdbtk de_novo_wf --genome_dir $GENOMES_DIR --archaea --out_dir $OUTDIR --outgroup_taxon p__Altarchaeota --extension fa --cpus 16
