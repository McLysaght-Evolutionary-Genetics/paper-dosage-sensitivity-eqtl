#!/bin/bash

####################################################################################################
#############################################trans-eQTL detection pipeline##########################
####################################################################################################
#Set the file directory
repeatmasker="~/ARS-UCD1.2/ARS-UCD1.2.repeatmasker.bed"
simplerepeat="~/ARS-UCD1.2/ARS-UCD1.2.simple_repeats.bed"
trans="~/eQTL/trans"
used_genes="~/ARS-UCD1.2/ARS1.2.genebody.protein_coding-lncRNA.bed"
mappability="~/eQTL/trans/Cross_map/Cow_alignability/gene_mappability/gene_mappability.txt"
snp_mappability="~/eQTL/trans/Cross_map/Cow_alignability/bosTau9_mappability_75mer_2mismatch.nochr.map_1.bed"
	
cd /home/shuli.liu/bull_scr/RNA-seq/SNP_calling/Final_vcf_file/filter_MAF_0.05_DR2_0.8/eQTL/vcf_tissues
tissue="Liver"
cd ${trans}
mkdir ${tissue}
#Remove the SNPs overlapped with repeats(repeatmaskers and simple repeats); SNPs with minor allele frequency less than 0.05 and keep the SNPs with at least 10 minor allele samples.
bcftools view -T ^${repeatmasker} ${tissue}.vcf > ${trans}/${tissue}/${tissue}.filtered2.vcf
cd ${trans}/${tissue}
bcftools view -T ^${simplerepeat} ${tissue}.filtered2.vcf > ${tissue}.filtered1.vcf
bcftools view -q 0.05:minor -v snps -c 10:minor ${tissue}.filtered1.vcf > ${tissue}.filtered0.vcf

##Only keep the independent SNPs
plink --indep 50 5 2 --double-id --vcf ${tissue}.filtered0.vcf --chr-set 29 -out ${tissue}.subset.beagle
vcftools --vcf ${tissue}.filtered0.vcf --snps ${tissue}.subset.beagle.prune.in --recode --out ${tissue}.filtered

#Remove intermediate files
rm ${tissue}.filtered0.vcf
rm ${tissue}.filtered1.vcf
rm ${tissue}.filtered2.vcf

#Run the MatrixeQTL to get the trans-eQTLs.
Rscript ~/trans-eQTL.matrixeQTL.r ${tissue}

#####################################################################
#########################Filter trans-eQTLs##########################
#Only focus on protein_coding genes and lincRNA genes
sed '1d; s/"//g' ${tissue}.trans.txt | awk 'FNR==NR{ a[$7]; next}{k=$3; if ( k in a) print}' ${used_genes} - > filter_step1.${tissue}.trans.txt.tmp
#Only keep the mappability >= .8 genes.
awk 'FNR==NR{a[$1]=$2;next}{k=$3; if ( (k in a) && (a[k]!="NA") && (a[k]>0.8)) print }' ${mappability} filter_step1.${tissue}.trans.txt.tmp > filter_step2.${tissue}.trans.txt.tmp
#Filter out the gene-variant pairs where the cross-mapped genes of target gene within the 1Mb of the variant
Rscript ~/filter_out_genes_cross_mappability.r ${tissue}
#Filter out the gene-snp pairs
awk 'FNR==NR{a[$1"\t"$2];next}{k=$2"\t"$3; if (!(k in a)) print}' ${tissue}.filter_out_cross_map.txt.tmp filter_step2.${tissue}.trans.txt.tmp > ${tissue}.trans.sig.filtered.txt

#Only keep SNPs with mappability = 1
awk -F "[_\t]" '{print $2"\t"$3"\t"$3"\t"$0}' ${tissue}.trans.sig.filtered.txt | intersectBed -a - -b ${snp_mappability} > ${tissue}.trans.sig.filtered.final.txt

rm *.tmp

###################################################################################################
##################################Trans-eQTL detection: Adjust the fine-mapped cis-eQTLs###########
###################################################################################################
