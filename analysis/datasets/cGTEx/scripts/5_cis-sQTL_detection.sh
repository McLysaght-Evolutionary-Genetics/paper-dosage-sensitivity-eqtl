# Take Liver as an example

# Step 1. Converting bams to juncs
for i in $(seq 1 10)
do
	sh leafcutter/scripts/bam2junc.sh liver${i}.bam liver${i}.junc
done


# Step 2. Intron clustering
ls *.junc > liver_juncfiles.txt
python leafcutter/clustering/leafcutter_cluster.py -j liver_juncfiles.txt -m 50 -o liver -l 500000


# Step 3. Calculate the PCA 
python leafcutter/scripts/prepare_phenotype_table.py liver_perind.counts.gz -p 10 # note: -p 10 specific you want to calculate for sQTL to use as covariates
sh liver_perind.counts.gz_prepare.sh
head -6 liver_perind.counts.gz.PCs > liver_perind.counts.gz.PCs.PC5


# Step 4. sQTL detection
for j in $(seq 1 29)
do
	fastQTL.static --vcf liver.filtered.vcf.gz --bed liver_perind.counts.gz.qqnorm_chr${j}.gz --cov liver_perind.counts.gz.PCs.PC5 --permute 1000 10000 --normal --out liver_perind.permutation.chr${j} --chunk 1 1
done

# go to R
library(qvalue)
d <- read.table("liver_perind.permutation.chrAll.txt", header = F, stringsAsFactors=F)
colnames(d) = c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope","ppval", "bpval")

pdf(paste0(Tissue,"liver_perind.permutation.pdf"),width=6,height=6)
plot(d$ppval, d$bpval, xlab="Direct method", ylab="Beta approximation", main="Check plot")
abline(0, 1, col="red")
dev.off()

d$bh = p.adjust(d$bpval, method="fdr")
write.table(d[which(d$bh <= 0.05), ], "liver_perind.permutation.fdr.txt", quote = F, row.names=F, col.names = T, sep = "\t")
d$st = qvalue(d$bpval,lambda=0.5)$qvalues
write.table(d[which(d$st <= 0.05), ], "liver_perind.permutation.storey.txt", quote = F, row.names=F, col.names = T, sep = "\t")


# Step 5. Ref genes
cat liver_perind.counts.gz.qqnorm_chr* > liver_perind.counts.gz.qqnorm_chrAll
bedtools intersect -a Bos_taurus.ARS-UCD1.2.98_gene.txt -b liver_perind.counts.gz.qqnorm_chrAll -wa -wb > liver_perind.counts.gz.qqnorm_chrAll_ref
