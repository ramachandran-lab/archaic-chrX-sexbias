#!/bin/sh

#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=

#SBATCH -e arrayjob-%a.err
#SBATCH -o arrayjob-%a.out

#SBATCH -t 4:00:00
#SBATCH --mem 4Gb

#SBATCH --array=1-1835

source ~/hmmix_env/bin/activate
module load python/3.9.0 bcftools/1.9 vcftools/0.1.16

ingroup_hap=`cat ingroup_haps.txt | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

echo "Ingroup_hap: "$ingroup_hap

# # ensure folders exist
# mkdir obs_files
# mkdir trained_files
# mkdir decoded_autosomes
# mkdir decoded_chrXs  # would need to change for diploid mode.  see decode_chrX.py

# create_ingroup
echo "create_ingroup is starting"
ingroup_file="obs.${ingroup_hap}.txt"
hmmix create_ingroup -ind=$ingroup_hap -out=$ingroup_file -vcf='bcf_local/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,bcf_local/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,bcf_local/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,bcf_local/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,bcf_local/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,bcf_local/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,bcf_local/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,bcf_local/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,bcf_local/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,bcf_local/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,bcf_local/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,bcf_local/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,bcf_local/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,bcf_local/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,bcf_local/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,bcf_local/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,bcf_local/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,bcf_local/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,bcf_local/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,bcf_local/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,bcf_local/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,bcf_local/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf,bcf_local/ALL.chrX.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf' -out=obs -ancestral='anc_files/homo_sapiens_ancestor_10.fa,anc_files/homo_sapiens_ancestor_11.fa,anc_files/homo_sapiens_ancestor_12.fa,anc_files/homo_sapiens_ancestor_13.fa,anc_files/homo_sapiens_ancestor_14.fa,anc_files/homo_sapiens_ancestor_15.fa,anc_files/homo_sapiens_ancestor_16.fa,anc_files/homo_sapiens_ancestor_17.fa,anc_files/homo_sapiens_ancestor_18.fa,anc_files/homo_sapiens_ancestor_19.fa,anc_files/homo_sapiens_ancestor_1.fa,anc_files/homo_sapiens_ancestor_20.fa,anc_files/homo_sapiens_ancestor_21.fa,anc_files/homo_sapiens_ancestor_22.fa,anc_files/homo_sapiens_ancestor_2.fa,anc_files/homo_sapiens_ancestor_3.fa,anc_files/homo_sapiens_ancestor_4.fa,anc_files/homo_sapiens_ancestor_5.fa,anc_files/homo_sapiens_ancestor_6.fa,anc_files/homo_sapiens_ancestor_7.fa,anc_files/homo_sapiens_ancestor_8.fa,anc_files/homo_sapiens_ancestor_9.fa,anc_files/homo_sapiens_ancestor_X.fa' -weights=strickmask.bed -outgroup=outgroup.txt

# train
echo "train is starting"
trained_file="trained_files/trained.${ingroup_hap}.haploid.json"
hmmix train -obs=$ingroup_file -weights=strickmask.bed -mutrates=mutationrate.bed -out=$trained_file -haploid

# decode autosomes
echo "decode is starting for autosomes"
decoded_file="decoded_autosomes/decoded.${ingroup_hap}.autosome"
hmmix decode -obs=$ingroup_file -weights=strickmask.bed -mutrates=mutationrate.bed -param=$trained_file -out=$decoded_file -haploid


mv $ingroup_file obs_files/


# decode chrX
echo "decode is starting for chrX using decode_chrX.py"
python decode_chrX.py $ingroup_hap haploid
