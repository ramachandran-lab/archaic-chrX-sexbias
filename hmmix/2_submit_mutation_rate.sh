#!/bin/sh

#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=


source ~/hmmix_env/bin/activate
module load python/3.9.0 bcftools/1.9 vcftools/0.1.16

hmmix mutation_rate -outgroup=outgroup.txt -weights=strickmask.bed -window_size=100000 -out=mutationrate.bed
