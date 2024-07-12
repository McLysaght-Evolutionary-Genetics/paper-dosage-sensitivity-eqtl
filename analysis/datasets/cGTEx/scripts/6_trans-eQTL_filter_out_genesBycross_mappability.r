library("data.table")
tissue="Adipose"
cross_mappability<-read.table("/home/shuli.liu/bull_scr/RNA-seq/SNP_calling/Final_vcf_file/filter_MAF_0.05_DR2_0.8/eQTL/trans/Cross_map/Cow_alignability/cross_mappability3/cross_mappability.1Mb.genebody.txt")
step2<-read.table(paste0("filter_step2.",tissue,".trans.txt.tmp"))
colnames(cross_mappability)<-c("gene1","gene2","score","gene2_chr","gene2_1mbup","gene2_1mbdown")
colnames(step2)<-c("id","snp","gene","quality","pvalue","fdr","slope")
filter<-merge(step2,cross_mappability,by.x="gene",by.y="gene1")
filter<-as.data.frame(filter)
as.data.frame(setDT(filter)[,c("gene1_chr","gene1_snp"):=tstrsplit(snp,"_",keep=c(1,2))])
filter_genes<-filter[gene1_chr==gene2_chr & gene1_snp>=gene2_1mbup & gene1_snp<=gene2_1mbdown,]
gene_snp_pair<-filter_genes[,c("gene","snp")]
write.table(gene_snp_pair,paste0(tissue,"filter_out_cross_map.txt.tmp"),sep="\t",col.names=F,row.names=F,quote =FALSE)
