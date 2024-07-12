#!/bin/bash


######################################################################
######################################################################
##############1. Get Reference metadata
######################################################################
######################################################################


REPO="~/summary-gwas-imputation/src"
DATA="~/RNA-seq/run7_yak_indicus"

module load bcftools

#Run chromosome 1
chr="1"

#convert vcf into genotype matrix.
filter_and_convert ()
{
echo -ne "varID\t" | gzip > $2
bcftools view $1 --types snps --force-samples -Ou |  bcftools query -l | tr '\n' '\t' | sed 's/\t$/\n/' | gzip >> $2

NOW=$(date +%Y-%m-%d/%H:%M:%S)
echo "Starting at $NOW"
bcftools view $1 --types snps --force-samples -Ou | \
bcftools annotate --set-id +'chr%CHROM\_%POS\_%REF\_%FIRST_ALT' -Ou | \
bcftools query -f '%ID[\t%GT]\n' | \
awk '
{
for (i = 1; i <= NF; i++) {
    if (substr($i,0,1) == "c") {
        printf("%s",$i)
    } else if ( substr($i, 0, 1) == ".") {
        printf("\tNA")
    } else if ($i ~ "[0-9]|[0-9]") {
        n = split($i, array, "|")
        printf("\t%d",array[1]+array[2])
    } else {
        #printf("\t%s",$i)
        printf("Unexpected: %s",$i)
        exit 1
    }
}
printf("\n")
}
' | gzip >> $2

NOW=$(date +%Y-%m-%d/%H:%M:%S)
echo "Ending at $NOW"
}

filter_and_convert $DATA/genotypes/Chr${chr}-Run7-TAU-Beagle-toDistribute.vcf.gz $DATA/parquet_run7/Chr${chr}-Run7-TAU-Beagle-toDistribute.txt.gz

#Get the annotation files. chromosome position id allele_0 allele_1 allele_1_frequency rsid
format_annotation ()
{
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\t%CHROM:%POS\n' $1 | awk 'FNR==NR{a[$2]=$1;next}{k=$6; if (k in a)print $1, $2, "chr"$1"_"$2"_"$3"_"$4, $3, $4, $5, a[k]}' OFS="\t" $2 - | \
sed '1i chromosome\tposition\tid\tallele_0\tallele_1\tallele_1_frequency\trsid' | gzip > $3
}

format_annotation $DATA/genotypes/Chr${chr}-Run7-TAU-Beagle-toDistribute.vcf.gz /home/shuli.liu/uvm_mckay/WGBS_others/tissue-specific-MHL/SNPS/ARS1.2_BQSR.snp.rsid.txt $DATA/parquet_run7/Chr${chr}-Run7-TAU-Beagle-toDistribute.annot.txt.gz

#variant selection/Variant metadata
python $REPO/get_reference_metadata.py \
-genotype $DATA/parquet_run7/Chr${chr}-Run7-TAU-Beagle-toDistribute.txt.gz \
-annotation $DATA/parquet_run7/Chr${chr}-Run7-TAU-Beagle-toDistribute.annot.txt.gz \
-filter MAF 0.01 \
-filter TOP_CHR_POS_BY_FREQ \
-rsid_column rsid \
-output $DATA/parquet_run7/Chr${chr}_maf0.01_monoallelic_variants.txt.gz

#Arrange the metadata file
zcat Chr1_maf0.01_monoallelic_variants.txt.gz | sed -n '1p' > variant_metadata.txt
for chr in $(seq 1 29)
do
gzip -cd Chr${chr}_maf0.01_monoallelic_variants.txt.gz | sed '1d' >> variant_metadata.txt
done
gzip -c variant_metadata.txt > variant_metadata.txt.gz

#Arrange the variant file
zcat Chr1-Run7-TAU-Beagle-toDistribute.txt.gz | sed -n '1p' > Run7-TAU-Beagle-toDistribute.txt
for chr in $(seq 1 29)
do
gzip -cd Chr${chr}-Run7-TAU-Beagle-toDistribute.txt.gz | sed '1d' >>  Run7-TAU-Beagle-toDistribute.txt
done
gzip -c Run7-TAU-Beagle-toDistribute.txt > Run7-TAU-Beagle-toDistribute.txt.gz

#Genotype compilation
python $REPO/model_training_genotype_to_parquet.py \
-input_genotype_file $DATA/parquet_run7/Run7-TAU-Beagle-toDistribute.txt.gz \
-snp_annotation_file $DATA/parquet_run7/variant_metadata.txt.gz METADATA \
-parsimony 9 \
--impute_to_mean \
--only_in_key \
--split_by_chromosome \
-rsid_column rsid \
-output_prefix $DATA/parquet_run7/ARS_UCD1.2_maf0.01_monoallelic_variants




######################################################################
######################################################################
##############2. Homonization
######################################################################
######################################################################


##GWAS_TOOL Folder##
GWAS_TOOLS="~/summary-gwas-imputation/src"
##GWAS File Folder##
DATA="~/Cattle_GWAS_ARS1.2"
Reference="~/run7_yak_indicus/parquet_run7"
OUTPUT="~/2_RUN_SprediXcan/Prepare_GWAS"

#get gwas file and sample sizes for each gwas
cd ${OUTPUT}
gwas=`cat GWAS_sample_size.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 1`
n_sample=`cat GWAS_sample_size.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 2`
zcat ${gwas} | awk 'FNR==NR{a[$1"\t"$2]=$5"\t"$4; next}{ k=$1"\t"$2; if ( k in a) print $4,$1,$2,$7, a[k], $5,$6}' OFS="\t" ${Reference}/variant_metadata.txt - | awk '{print $2"_"$3"_"$6"_"$5,$2"_"$3"_"$6"_"$5,$2,$3,$4,$5,$6,$7,$8}' OFS="\t" | sed '1i variant_id\tpanel_variant_id\tchromosome\tposition\tpvalue\teffect_allele\tnon_effect_allele\teffect_size\tstandard_error' | gzip > ${OUTPUT}/harmonized_gwas/HM_${gwas}

##Harmonization 1
python $GWAS_TOOLS/gwas_parsing.py \
-gwas_file $OUTPUT/harmonized_gwas/HM_${gwas} \
-snp_reference_metadata $Reference/variant_metadata.txt.gz METADATA \
-output_column_map variant_id variant_id \
-output_column_map non_effect_allele non_effect_allele \
-output_column_map effect_allele effect_allele \
-output_column_map effect_size effect_size \
-output_column_map pvalue pvalue \
-output_column_map chromosome chromosome \
--chromosome_format \
-output_column_map position position \
-output_column_map standard_error standard_error \
--insert_value sample_size ${n_sample} \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size \
-output $OUTPUT/harmonized_gwas/HM2_${gwas}




######################################################################
######################################################################
##############3. SNP imputation
######################################################################
######################################################################


date

##GWAS_TOOL Folder##
GWAS_TOOLS="~/summary-gwas-imputation/src"

##GWAS File Folder##
Reference="~/run7_yak_indicus/parquet_run7"
OUTPUT="~/2_RUN_SprediXcan/Prepare_GWAS"
LD="~/Bovine_genome"

##Run the imputation. sub batches: 10
trait=$1
for chr in $(seq 1 29)
do
for i in $(seq 0 9)
do
python $GWAS_TOOLS/gwas_summary_imputation.py \
-gwas_file $OUTPUT/harmonized_gwas/HM2_${trait}.* \
-by_region_file ${LD}/genome.1mb.win.bed.gz \
-parquet_genotype ${Reference}/ARS_UCD1.2_monoallelic_variants.chr${chr}.variants.parquet \
-parquet_genotype_metadata ${Reference}/ARS_UCD1.2_monoallelic_variants.variants_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome ${chr} \
-regularization 0.1 \
-sub_batches 10 \
-sub_batch ${i} \
--standardise_dosages \
-output $OUTPUT/summary_imputation/${trait}.chr${chr}_${i}.txt.gz
done
done

python $GWAS_TOOLS/gwas_summary_imputation_postprocess.py \
-gwas_file $OUTPUT/harmonized_gwas/HM2_${trait}.*.gz \
-folder $OUTPUT/summary_imputation \
-pattern ${trait}.* \
-parsimony 7 \
-output $OUTPUT/processed_summary_imputation/imputed_${trait}.txt.gz
date

#Because our prediction model use CHR_POS_ALLELE1_ALLELE2 as the rsid; so change the rsid here into CHR_POS_ALLELE1_ALLELE2 format.
zcat $OUTPUT/processed_summary_imputation/imputed_${trait}.txt.gz | awk '{if(NR==1) print $0; else print $2,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' OFS="\t" | sed 's/^chr//g' | gzip > $OUTPUT/processed_summary_imputation/rsid_imputed_${trait}.txt.gz




######################################################################
######################################################################
##############3. Run S-PrediXcan
######################################################################
######################################################################


##set the directory for SPrediXcan and other folders
SPrediXcan="~/RNA-seq/TWAS/MetaXcan/software/SPrediXcan.py" #SPrediXcan.py
GWAS_dir="~/RNA-seq/TWAS/2_RUN_SprediXcan/Prepare_GWAS/processed_summary_imputation" #Imputated GWAS
model_dir="~/RNA-seq/TWAS/PredictDB_Pipeline_GTEx_v7/model_training" #Prediction model: cross-validated elastic net model
result_dir="~/RNA-seq/TWAS/2_RUN_SprediXcan/Prepare_GWAS/TWAS_results"# folder to store the TWAS results

cd ${GWAS_dir}

gwas=$1
for tissue in `cat ${model_dir}/scripts/gtex_tissues.txt`
do
python ${SPrediXcan} \
--model_db_path ${model_dir}/dbs/gtex_v1_${tissue}_imputed_Bos_signif.db \
--covariance ${model_dir}/covariances/gtex_v1_${tissue}_imputed_Bos_covariances.txt.gz \
--gwas_file rsid_${gwas} \
--snp_column variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--pvalue_column pvalue \
--keep_non_rsid \
--overwrite \
--output_file ${result_dir}/${tissue}/twas_${gwas}
done


















