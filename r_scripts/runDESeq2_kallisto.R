setwd("..")

############ Functions ############

autoGparFontSizeMatrix<-function(n){ #Calcule automatiquement la taille de police selon le nombre de colonnes ou lignes (empirique)
  n=max(n,50)
  n=min(n,1000)
  return(gpar(fontsize=1/n*600))
}

corrDist<-function(x) return(as.dist((1 - cor(Matrix::t(x)))/2))

rowScale<-function(data,center=TRUE,scaled=FALSE){
  data<-t(data)
  data<-t(scale(data,center=center,scale=scaled))
  return(data)
}

unsupervisedClustering<-function(x,transpose=TRUE,method.dist="pearson",method.hclust="ward.D2",bootstrap=FALSE,nboot=10){
  if(transpose) x<-t(x)
  if(bootstrap){
    require(pvclust)
    resClust<-pvclust(t(x),nboot=nboot,method.hclust = method.hclust,parallel = TRUE,method.dist = method.dist)$hclust
  }else{
    if(method.dist=="pearson"){
      resDist<-corrDist(x)
    }else{
      resDist<-dist(x, method = method.dist)
    }
    resClust<-stats::hclust(resDist,method = method.hclust)
  }
  return(resClust)
}

############ Packages ############

require(DESeq2)
library(fdrtool)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(IHW)
library(tximport)
library(dplyr)
library(tidyr)
library(svglite)

############ Arguments & data ############

args <- commandArgs(trailingOnly = TRUE)

INFILE <- args[1]
OUTDIR <- args[2]
LOGFCTHRESHOLD <- 1
FDRTHRESHOLD <- 0.05

#############################
# 		with kallisto 		#
#############################

palette <- c("#8fd7d7","#bdd373","#ffcd8e","#00b0be","#98c127","#ffb255")
names(palette) <- c('RH30_CTL_14h','RH30_APA_14h','RH30_D11_14h','RH30_CTL_48h','RH30_APA_48h','RH30_D11_48h')

sampleTable <- read.table(INFILE, header=TRUE,sep = '\t')
conditions<-factor( sampleTable[ , 3] )
uniq_conds <- unique(conditions)

files <- sampleTable[,2]
names(files) <- sampleTable[,1]
txi.kallisto <- tximport(files, type="kallisto", txOut=T)
dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~condition)
dds <- DESeq(dds)

matrix <- counts(dds,normalized=TRUE)
# Transforming raw counts
rld <- rlogTransformation(dds, blind=TRUE)
mrld <- assay(rld)

## Component plot
rv = rowVars(assay(rld))
select = order(rv, decreasing=TRUE)[seq_len(min(1000, length(rv)))]
pca = prcomp(t(assay(rld)[select,]))
  
variance = pca$sdev^2 / sum(pca$sdev^2)
variance = round(variance, 3) * 100
data_pca = as.data.frame(pca$x)
data_pca$samplename = rownames(data_pca)
data_pca = merge(data_pca,sampleTable[,c(1,3)])

outPlot = c(paste(OUTDIR,"/PCAplot_with_names.pdf",sep=""))
pdf(outPlot,width=5.9,height=5.1)
ggplot(data_pca, aes(x=PC1,y=PC2, label=samplename, color=condition)) +
	geom_point(size=4) +
	geom_text_repel(color="black") +
	theme_light() +
	xlab(paste(c("PC1 (variance = ",variance[1],"%)"),collapse="")) +
	ylab(paste(c("PC2 (variance = ",variance[2],"%)"),collapse="")) +
	guides(color=guide_legend(title="Condition")) + 
	scale_colour_manual(values=palette)
dev.off()

outPlot = c(paste(OUTDIR,"/PCAplot_without_names.pdf",sep=""))
pdf(outPlot,width=5.9,height=5.1)
ggplot(data_pca, aes(x=PC1,y=PC2, label=samplename, color=condition)) +
	geom_point(size=5) +
	theme_light() +
	xlab(paste(c("PC1 (variance = ",variance[1],"%)"),collapse="")) +
	ylab(paste(c("PC2 (variance = ",variance[2],"%)"),collapse="")) +
	guides(color=guide_legend(title="Condition")) + 
	scale_colour_manual(values=palette)
dev.off()

########## Boxplot genes ########

# correspondance_symbol <- read.table(args[3],sep="\t",h=F)
# fpkm_exp <- fpkm(dds)

# gene_selected <- "SNAI1"
# count_norm <- fpkm_exp
# transcripts_selected <- correspondance_symbol[which(correspondance_symbol$V1 == gene_selected),2]
# data_boxplot <- t(count_norm[which(rownames(count_norm) %in% transcripts_selected),])
# data_boxplot_2 <- data.frame('expression'=c(),'sample'=c(),'condition'=c(),'transcript'=c())

# for (transcript in colnames(data_boxplot)){
#   sub_data <- data.frame('expression'= as.numeric(data_boxplot[,transcript]),'sample'=sampleTable[,1],'condition'=sampleTable[,3],'transcript'=transcript)
#   data_boxplot_2 <- rbind(data_boxplot_2,sub_data)
# }
# data_boxplot <- data_boxplot_2

# outPlot = paste(OUTDIR,"/gene_SNAI1_expression.pdf",sep="")
# pdf(outPlot,width=8,height=7)
# ggplot(data_boxplot, aes(x = condition, y = expression, fill=condition)) +
#   geom_boxplot(alpha=0.8) +
#   theme_minimal() +
#   guides(size="none") +
#   xlab("Conditions") +
#   ylab(expression('Normalized expression of '~italic(SNAI1)~' transcripts (FPKM)')) +
#   facet_wrap(~transcript) +
#   ggtitle(expression(italic(SNAI1)~' Transcripts')) +
#   theme(plot.title=element_text(hjust=0.5), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
#   scale_fill_manual(values=palette) +
#   guides(fill=guide_legend(title="Condition"))
# dev.off()

# gene_selected <- "SNAI2"
# count_norm <- fpkm_exp
# transcripts_selected <- correspondance_symbol[which(correspondance_symbol$V1 == gene_selected),2]
# data_boxplot <- t(count_norm[which(rownames(count_norm) %in% transcripts_selected),])
# data_boxplot_2 <- data.frame('expression'=c(),'sample'=c(),'condition'=c(),'transcript'=c())

# for (transcript in colnames(data_boxplot)){
#   sub_data <- data.frame('expression'= as.numeric(data_boxplot[,transcript]),'sample'=sampleTable[,1],'condition'=sampleTable[,3],'transcript'=transcript)
#   data_boxplot_2 <- rbind(data_boxplot_2,sub_data)
# }
# data_boxplot <- data_boxplot_2

# outPlot = paste(OUTDIR,"/gene_SNAI2_expression.pdf",sep="")
# pdf(outPlot,width=8,height=7)
# ggplot(data_boxplot, aes(x = condition, y = expression, fill=condition)) +
#   geom_boxplot(alpha=0.8) +
#   theme_minimal() +
#   guides(size="none") +
#   xlab("Conditions") +
#   ylab(expression('Normalized expression of '~italic(SNAI2)~' transcripts (FPKM)')) +
#   facet_wrap(~transcript) +
#   ggtitle(expression(italic(SNAI2)~' Transcripts')) +
#   theme(plot.title=element_text(hjust=0.5), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
#   scale_fill_manual(values=palette) +
#   guides(fill=guide_legend(title="Condition"))
# dev.off()

# gene_selected <- "KCNN3"
# count_norm <- fpkm_exp
# transcripts_selected <- correspondance_symbol[which(correspondance_symbol$V1 == gene_selected),2]
# data_boxplot <- t(count_norm[which(rownames(count_norm) %in% transcripts_selected),])
# data_boxplot_2 <- data.frame('expression'=c(),'sample'=c(),'condition'=c(),'transcript'=c())

# for (transcript in colnames(data_boxplot)){
#   sub_data <- data.frame('expression'= as.numeric(data_boxplot[,transcript]),'sample'=sampleTable[,1],'condition'=sampleTable[,3],'transcript'=transcript)
#   data_boxplot_2 <- rbind(data_boxplot_2,sub_data)
# }
# data_boxplot <- data_boxplot_2

# outPlot = paste(OUTDIR,"/gene_KCNN3_expression.pdf",sep="")
# pdf(outPlot,width=8,height=7)
# ggplot(data_boxplot, aes(x = condition, y = expression, fill=condition)) +
#   geom_boxplot(alpha=0.8) +
#   theme_minimal() +
#   guides(size="none") +
#   xlab("Conditions") +
#   ylab(expression('Normalized expression of '~italic(KCNN3)~' transcripts (FPKM)')) +
#   facet_wrap(~transcript) +
#   ggtitle(expression(italic(KCNN3)~' Transcripts')) +
#   theme(plot.title=element_text(hjust=0.5), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
#   scale_fill_manual(values=palette) +
#   guides(fill=guide_legend(title="Condition"))
# dev.off()

# gene_selected <- "PAX3FOXO1"
# count_norm <- fpkm_exp
# #transcripts_selected <- correspondance_symbol[which(correspondance_symbol$V1 == gene_selected),2]
# data_boxplot <- t(count_norm[which(rownames(count_norm) == "PAX3FOXO1"),])
# data_boxplot_2 <- data.frame('expression'=c(),'sample'=c(),'condition'=c(),'transcript'=c())

# for (transcript in colnames(data_boxplot)){
#   sub_data <- data.frame('expression'= as.numeric(data_boxplot[,transcript]),'sample'=sampleTable[,1],'condition'=sampleTable[,3],'transcript'=transcript)
#   data_boxplot_2 <- rbind(data_boxplot_2,sub_data)
# }
# data_boxplot <- data_boxplot_2

# outPlot = paste(OUTDIR,"/gene_PAX3FOXO1_expression.pdf",sep="")
# pdf(outPlot,width=8,height=7)
# ggplot(data_boxplot, aes(x = condition, y = expression, fill=condition)) +
#   geom_boxplot(alpha=0.8) +
#   theme_minimal() +
#   guides(size="none") +
#   xlab("Conditions") +
#   ylab(expression('Normalized expression of '~italic(PAX3FOXO1)~' transcripts (FPKM)')) +
#   facet_wrap(~transcript) +
#   ggtitle(expression(italic(PAX3FOXO1)~' Transcripts')) +
#   theme(plot.title=element_text(hjust=0.5), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
#   scale_fill_manual(values=palette) +
#   guides(fill=guide_legend(title="Condition"))
# dev.off()

######## Differential transcripts ###########

condCol="condition"
logFCthreshold=1
AdjPValthreshold=0.05
GenesInFig=50
bootstrap=FALSE
nboot=30

pw_tests <- list()
pw_tests[[1]] <- c('RH30_APA_14h','RH30_CTL_14h')
pw_tests[[2]] <- c('RH30_D11_14h','RH30_CTL_14h')
pw_tests[[3]] <- c('RH30_APA_48h','RH30_CTL_48h')
pw_tests[[4]] <- c('RH30_D11_48h','RH30_CTL_48h')

for (i in 1:length(pw_tests)){
  cond1 <- pw_tests[[i]][1]
  cond2 <- pw_tests[[i]][2]

  comp <- paste(cond1,cond2,sep="__vs__")

  sampleAnnot <- as.data.frame(sampleTable[,3])
  rownames(sampleAnnot) <- sampleTable[,2]
  colnames(sampleAnnot) <- "condition"
  samples <- rownames(sampleAnnot)

  res <- results(dds,contrast=c(condCol,cond1,cond2),independentFiltering=T)
  res$meanInComp <- rowMeans(mrld[,sampleAnnot[,condCol]%in% c(cond1,cond2)])
  res$padj <- adj_pvalues(ihw(pvalues = res$pvalue,covariates = res$meanInComp,alpha=AdjPValthreshold))

  DE <- data.frame(res)

  DE.sel <- list()
  DE.sel$up <- DE[which(DE$padj < AdjPValthreshold & DE$log2FoldChange > logFCthreshold),]
  DE.sel$down <- DE[which(DE$padj < AdjPValthreshold & DE$log2FoldChange < -logFCthreshold),]

  DE.sel$isDE=rbind(DE.sel$up,DE.sel$down)
  DE.sel$notDE=DE[setdiff(rownames(DE),rownames(DE.sel$isDE)),]

  DE$DE="NONE"
  DE[rownames(DE.sel$up),"DE"]="UP"
  DE[rownames(DE.sel$down),"DE"]="DOWN"
  DE$DE=factor(DE$DE,levels=c("DOWN","NONE","UP"))
  DE$transcript <- rownames(DE)
  correspondance_symbol <- read.table(args[3],sep="\t",h=T)
  DE <- merge(DE,correspondance_symbol)
  write.table(DE, paste(paste(OUTDIR,comp,sep="/"),"DE.tsv",sep="_"),sep="\t",row.names=F,quote=F)


  DE$DE_label <- NA
  list_genes <- c("KCNN3","SNAI1","SNAI2")
  DE$DE_label[which(DE$gene %in% list_genes)] <- DE$transcript[which(DE$gene %in% list_genes)]
  DE_palette = c("firebrick3", "gray60", "palegreen3")
  names(DE_palette) = c("DOWN","NONE","UP")

  pdf(paste(c(OUTDIR,"/",comp,"_Volcano_with_names.pdf"),collapse=""),width=5.9,height=5.9)
    print(ggplot(data=DE, aes(x=log2FoldChange, y=-log10(padj), color=DE, label=DE_label)) + 
    geom_point() +
    theme_minimal() +
    geom_text_repel(color="black", min.segment.length=0, nudge_x=2) +
    scale_color_manual(values=DE_palette))
  dev.off()

}


