#!/bin/bash

#######################################
####

subset="~/RNA-seq/SNP_calling/All_SNPs/sub_set_samples" #folder where stores the sample list. 
chrom="~/RNA-seq/SNP_calling/All_SNPs/chromsome.list.txt" #chromosome list file
Beagle="~/RNA-seq/beagle.21Sep19.ec3.jar" #Beagle software
run7="~/RNA-seq/run7/genotypes" #folder where stores the reference vcf file.
##########################################################
cd ${subset}
sub_set=`ls * | head -n $SLURM_ARRAY_TASK_ID | tail -n 1` #Beagle require a large amount of memory. We divided all samples into subsets and imputed them seperately.

#######################################
for region in 2:1-136231102
do
chr=`echo ${region} | cut -d ":" -f 1`
java -jar -Xmx120g -Djava.io.tmpdir=$TMPDIR ${Beagle} gt=${sub_set}.${region}.vcf.gz ref=${run7}/Chr${chr}-Run7-TAU-Beagle-toDistribute.vcf.gz chrom=${chr} nthreads=18 out=../After_imputation/${sub_set}.${chr}-beagle.vcf.gz
done


