#!/bin/bash

#########RNA SNP calling pipeline: based on GATK best practice############
folder_bam="~/bam_files/shuli"  #change as you want
folder_snp="~/SNP_calling/shuli"  # change as you want

cd ${folder_bam}

#Set up the file name, sample names,directories for references and software
bam="SRS2217961-STARAligned.sortedByCoord.out.bam"
sample="SRS2217961"
reference="~/Bovine_genome/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa"
dbSNP="~/SNPS/ARS1.2_BQSR.vcf" #vcf format
interval="~/SNPS/ARS1.2_BQSR.bed"
picard="/software/7/apps/picard/64/2.9.2" 
###########################################################################
##1. add read groups, sort, mark duplicates, and create index

cd ${folder_snp}

mkdir ${sample}
java -jar -Xmx15g -Djava.io.tmpdir=$TMPDIR ${picard}/picard.jar AddOrReplaceReadGroups I=${folder_bam}/${bam} O=./${sample}/rg_added_${bam} \
RGID=4 RGLB=lib1 RGPL=illumina RGPU=run RGSM=20 CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate

cd ${sample}
java -jar -Xmx15g -Djava.io.tmpdir=$TMPDIR ${picard}/picard.jar MarkDuplicates I=rg_added_${bam} O=dedupped_${bam}  \
CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=marked_dup_metrics.txt

##2. Split N trim and reassign mapping qualities
#Reassign OneMapping Quality: reassign all good alignments (MAPQ=225) to the default value of 60. 
gatk SplitNCigarReads \
   --spark-runner LOCAL \
   --TMP_DIR ${TMPDIR} \
      -R ${reference} \
      -I dedupped_${bam} \
      -O SplitNCigarReads_${bam}

 #rm dedupped_${bam} #remove

##3. Base recalibration (BQSR). ( used or not, determine based on the run time; the effect is marginal)
 gatk BaseRecalibrator \
   --spark-runner LOCAL \
   --TMP_DIR ${TMPDIR} \
   -I SplitNCigarReads_${bam} \
   -R ${reference} \
   --known-sites ${dbSNP} \
   -O ${sample}-recal.table

 gatk ApplyBQSR \
  --spark-runner LOCAL \
   -R ${reference} \
   -I SplitNCigarReads_${bam} \
   --bqsr-recal-file ${sample}-recal.table \
   -O BQSR_${bam}
  
 #rm SplitNCigarReads_${bam} #remove
 
##4. Run the haplotypecaller. 

 gatk HaplotypeCaller  \
   --spark-runner LOCAL \
   --TMP_DIR ${TMPDIR} \
   -R ${reference} \
   -I BQSR_${bam} \
   -O ${sample}.vcf.gz \
   -dbsnp ${dbSNP} \
    -L ${interval} \
   --dont-use-soft-clipped-bases \
   --output-mode EMIT_ALL_CONFIDENT_SITES \
   -stand-call-conf 0
   
# rm BQSR_${bam} #remove

 #5. Filter the snps to get the high quality ones
 
 gatk VariantFiltration \
 --TMP_DIR ${TMPDIR} \
 -R ${reference} \
 -V ${sample}.vcf.gz \
 -O ${sample}.filtered.vcf.gz \
 -window 35 -cluster 3 \
 --filter-name one \
 --filter-expression "FS>30.0" \
 --filter-name two \
 --filter-expression "QD<2.0" 

##6. Selected variants
gatk SelectVariants \
-R ${reference} \
-V  ${sample}.filtered.vcf.gz \
--exclude-filtered \
-O selected.${sample}.filtered.vcf.gz
 
##7. Get the ASE matrix.

gatk ASEReadCounter \
-R ${reference} \
-I rg_added_${bam} \
-V selected.${sample}.filtered.vcf.gz \
-min-depth 10 \
-O ${sample}.ASE.table.txt


#rm rg_added_${bam} #remove

 
   

   
   