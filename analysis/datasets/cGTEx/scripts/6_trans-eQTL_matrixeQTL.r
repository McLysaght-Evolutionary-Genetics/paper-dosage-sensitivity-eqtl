###function#########
moveme <- function (invec, movecommand) {
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]], 
                                 ",|\\s+"), function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x) {
    Where <- x[which(x %in% c("before", "after", "first", 
                              "last")):length(x)]
    ToMove <- setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp) - 1
      }
      else if (A == "after") {
        after <- match(ba, temp)
      }
    }
    else if (A == "first") {
      after <- 0
    }
    else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
}





##########Prepare covariance for fastQTL 



#########genotype#########
##use SNPRelate to calculate the pca; select the first four pcas.
library(SNPRelate)
tissue=commandArgs(T)[1]
k=commandArgs(T)[2]
vcf.fn<-paste0(tissue,".filtered.vcf")
snpgdsVCF2GDS(vcf.fn,"ccm.gds",method="biallelic.only")
genofile<-openfn.gds("ccm.gds")
ccm_pca<-snpgdsPCA(genofile)

if(length(ccm_pca$sample.id) < 150 ) {
pca_genotype<-ccm_pca$eigenvect[,1:3]
colnames(pca_genotype)<-c("pc1","pc2","pc3")
}else if (length(ccm_pca$sample.id)  < 250){
pca_genotype<-ccm_pca$eigenvect[,1:5]
colnames(pca_genotype)<-c("pc1","pc2","pc3","pc4","pc5")
}else{
pca_genotype<-ccm_pca$eigenvect[,1:10]
colnames(pca_genotype)<-c("pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")
}

rownames(pca_genotype)<-ccm_pca$sample.id

#####expression######
#Phenotype matrix: 
#Chromosome ID [string]
#Start genomic position of the phenotype (e.g. TSS) [integer]
#End genomic position of the phenotype (e.g. TTS) [integer]
#Phenotype ID [string] Then additional columns give phenotype quantifications for all samples. 
library(peer)
library(preprocessCore)
library(RNOmni)
load("/home/shuli.liu/uvm_mckay/WGBS_others/RNA-seq/all.LJA.8742samples.RData") 
expr=subset(TPM_tmp0,rownames(TPM_tmp0)%in%ccm_pca$sample.id)
expr_matrix00=as.data.frame(t(expr))
expr_matrix00<-expr_matrix00[,order(factor(colnames(expr_matrix00),levels=ccm_pca$sample.id))]
count_0.1<-rowSums(expr_matrix00>0.1)
#calculate the sample number:
nsamples=length(ccm_pca$sample.id)
#keep the genes with >0.1 tpm in >=20% samples.
expr_matrix00<-expr_matrix00[count_0.1>=(0.2*nsamples),] ##row is gene; column is sample
###expression values (TPM) were firstly quantile normalized between samples.and then inverse normal transformed across samples.
expr_matrix00_qn<-normalize.quantiles(as.matrix(expr_matrix00)) #normalize.quantiles needs column a chip (sample);row is the gene. normalize between samples. 

###expression for each gene were inverse normal distributed across samples.
INT<-function(x){
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
}
expr_matrix00_qn_ind<-t(apply(expr_matrix00_qn,MARGIN=1,FUN=INT))#apply to each row, each row represents one gene, observed values for all the samples. scale across samples.
colnames(expr_matrix00_qn_ind)<-colnames(expr_matrix00)
rownames(expr_matrix00_qn_ind)<-rownames(expr_matrix00)
expr_peer<-expr_matrix00_qn_ind
expr_matrix00_qn_ind<-as.data.frame(-expr_matrix00_qn_ind)
expr_matrix00_qn_ind$id<-row.names(expr_matrix00)


region_annot<-read.table("/home/shuli.liu/bismark_genome/ARS-UCD1.2/ARS1.2.TSS.bed")
colnames(region_annot)<-c("#Chr","start","end","name","score","strand","ID")
expr_matrix0<-merge(region_annot,expr_matrix00_qn_ind,by.x="ID",by.y="id")
expr_matrix<-expr_matrix0[,-which(names(expr_matrix0) %in% c("name","score","strand"))]
expr_matrix<-expr_matrix[moveme(names(expr_matrix),"ID after end")]

gene_expr<-expr_matrix[,-c(1:3)]
write.table(gene_expr,paste0(tissue,".gene.txt"),sep="\t",row.names=F,quote =FALSE)


model=PEER()
PEER_setPhenoMean(model,as.matrix(expr_peer)) #NULL response means no err# #N rows (samples); G columns (Genes)
dim(PEER_getPhenoMean(model))
  PEER_setNk(model,10)
  PEER_getNk(model)
  PEER_update(model)
  factors=PEER_getX(model)
  rownames(factors)<-row.names(as.data.frame(expr)) 
 #estimate how many cofounder factors to estimate #15 factors for < 150 samples; 30 factors for 150-250 samples; 35 factors for > 250 samples.


####prepare covariance file 
######merge above covariance### 
covariance0<-cbind(pca_genotype,factors)
#######Other covariance#######
##example: Species
info<-read.table("/home/shuli.liu/uvm_mckay/WGBS_others/RNA-seq/data_info_tissue_breed.txt",header=T,sep="\t")
tissue_info<-subset(info,info$Sample%in%ccm_pca$sample.id)
tissue_info<-tissue_info[order(factor(tissue_info$Sample, levels=ccm_pca$sample.id)),]
rownames(tissue_info)<-tissue_info$Sample
Species<-as.data.frame(tissue_info$Species)
colnames(Species)<-"Species"
##merge with above covariances
covariance<-t(cbind(Species,covariance0))
covariance<-as.data.frame(covariance)
covariance$id<-row.names(covariance)
covariance$Species<-as.character(covariance$Species)
covariance$Species[covariance$Species=="Bos taurus"]=0
covariance$Species[covariance$Species=="Bos indicus"]=1
covariance$Species[covariance$Species=="Bos indicus x Bos taurus"]=2
covariance$Species[covariance$Species=="Bos grunniens"]=3
covariance$Species[covariance$Species=="Bos grunniens x Bos taurus"]=4
covariance$Species[covariance$Species=="Bos frontalis"]=5
covariance$Species[covariance$Species=="Buffalo"]=6
covariance<-covariance[moveme(names(covariance),"id first")]
write.table(covariance,paste0(tissue,".cvrt.txt"),sep="\t",row.names=F,quote =FALSE)


#############################
###load files
library("MatrixEQTL")


snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(paste0(tissue,".snp.txt"));

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(paste0(tissue,".gene.txt"))

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$LoadFile(paste0(tissue,".cvrt.txt"))


###snpspos
snpspos=NULL
snpspos$snp<-read.gdsn(index.gdsn(genofile, "snp.rs.id"))
snpspos$chr<-read.gdsn(index.gdsn(genofile, "snp.chromosome"))
snpspos$pos<-read.gdsn(index.gdsn(genofile, "snp.position"))

##genepos
genepos<-expr_matrix[c(1:4)]
genepos<-genepos[moveme(names(genepos),"ID first")]
colnames(genepos)<-c("geneid","chr","s1","s2")

###########################################################
###run matrix eQTL

# Output file name
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis =1e-5;
pvOutputThreshold_tra =1e-5;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;
# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR;

me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = NULL,
  pvOutputThreshold.cis = NULL,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = TRUE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);


## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
#show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
#show(me$trans$eqtls)
min.pv.gene<-me$trans$min.pv.gene
write.table(as.data.frame(min.pv.gene),paste(tissue,".trans.min.pv.gene.txt"),sep="\t")
write.table(as.data.frame(me$trans$eqtls[me$trans$eqtls$FDR<=0.05,]),paste(tissue,".trans.txt"),sep="\t")

