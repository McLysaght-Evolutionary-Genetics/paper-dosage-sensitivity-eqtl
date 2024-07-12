#!/bin/bash

##########################################################################################################################
################################################################cis-eQTL detection pipeline###############################
##########################################################################################################################

#Take Liver tissue as an example
tissue="Liver"

#minor allele frequencies >=0.01; with the minor allele observed in at least 4 samples.
bcftools view -q 0.01:minor -v snps -c 4:minor ./vcf_tissues/${tissue}.vcf | bgzip -cf > ${tissue}.filtered.vcf.gz
tabix -p vcf ${tissue}.filtered.vcf.gz

#Prepare the phenotype matrix (gene expression matrix) and covariance matrix
Rscript ~/prepare_for_eQTL.detection.r ${tissue} --no-save

sort -k1,1 -k2,2n ${tissue}.phenotype.bed | bgzip -cf > ${tissue}.phenotype.bed.gz
tabix -p bed ${tissue}.phenotype.bed.gz
bgzip -cf ${tissue}.covariance.txt > ${tissue}.covariance.txt.gz



##################################
#Run FastQTL: Nominal mode
for j in $(seq 1 60); 
do
fastQTL.static --vcf ${tissue}.filtered.vcf.gz --bed ${tissue}.phenotype.bed.gz --cov ${tissue}.covariance.txt.gz  --normal --out ./Nominal/${tissue}.nominals.chunk${j}.txt.gz --chunk $j 60
done

#Combine chunks into one file
cd ./Nominal
zcat ${tissue}.nominals.chunk*.txt.gz | gzip -c > ${tissue}.nominals.txt.gz
rm ${tissue}.nominals.chunk*.txt.gz

#################################
#Run FastQTL: Permutation mode
for j in $(seq 1 60); 
do
fastQTL.static --vcf ${tissue}.filtered.vcf.gz --bed ${tissue}.phenotype.bed.gz --cov ${tissue}.covariance.txt.gz  --permute 1000 10000 --normal --out ./Permutation/${tissue}.permutations.chunk${j}.txt.gz --chunk $j 60& 
done

#Combine chunks into one file
cd ./Permutation
zcat ${tissue}.permutations.chunk*.txt.gz | gzip -c > ${tissue}.permutations.txt.gz
rm ${tissue}.permutations.chunk*.txt.gz

#################################
#set the criteria for eQTLs in each eGene.
Rscript ~/eQTL.p-value.nominal.correction_basePermutation.r ./Permutation/${tissue}.permutations.txt.gz 0.05 ./Nominal/${tissue}.nominals.txt.gz ${tissue}.nominals.sig.txt
