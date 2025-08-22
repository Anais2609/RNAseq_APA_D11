setwd("..")

############ Functions ############

autoGparFontSizeMatrix<-function(n){ #Calcule automatiquement la taille de police selon le nombre de colonnes ou lignes (empirique)
  n=max(n,50)
  n=min(n,1000)
  return(gpar(fontsize=1/n*500))
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

suppressMessages(library(DESeq2))
suppressMessages(library(fdrtool))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(IHW))
suppressMessages(library(svglite))
suppressMessages(library(circlize))

############ Arguments & data ############

args <- commandArgs(trailingOnly = TRUE)

INFILE <- args[1]
COUNTDIRECTORY <- args[2]
OUTDIR <- args[3]
LOGFCTHRESHOLD <- 1
FDRTHRESHOLD <- 0.05

# Data
sampleTable <- read.table(INFILE, header=TRUE,sep = '\t')
directory <- c(COUNTDIRECTORY)


############ Data preparation ############ 

# List of all the combination possible of paired conditions
conditions<-factor( sampleTable[ , 3] )
uniq_conds <- unique(conditions)

palette <- c("#8fd7d7","#bdd373","#ffcd8e","#00b0be","#98c127","#ffb255")
names(palette) <- c('RH30_CTL_14h','RH30_APA_14h','RH30_D11_14h','RH30_CTL_48h','RH30_APA_48h','RH30_D11_48h')

############ Data manipulation ############ 

# Counts recovery
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design= ~condition)

# deseq2 analysis
ddsHTSeq <- DESeq(ddsHTSeq)

# Normalized counts matrix
matrix <- counts(ddsHTSeq,normalized=TRUE)

# Transforming raw counts
rld <- rlogTransformation(ddsHTSeq, blind=TRUE)
mrld <- assay(rld)
write.table(cbind(rep(NA,nrow(matrix)),matrix),paste(OUTDIR,"/../../RNAseq_APA_D11_expression_v0.gct",sep=""),sep="\t",quote=F)

############ Plots ############

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
  scale_color_manual(values=palette)
dev.off()

outPlot = c(paste(OUTDIR,"/PCAplot_without_names.pdf",sep=""))
pdf(outPlot,width=5.9,height=5.1)
ggplot(data_pca, aes(x=PC1,y=PC2, label=samplename, color=condition)) +
	geom_point(size=5) +
	theme_light() +
	xlab(paste(c("PC1 (variance = ",variance[1],"%)"),collapse="")) +
	ylab(paste(c("PC2 (variance = ",variance[2],"%)"),collapse="")) +
	guides(color=guide_legend(title="Condition"))+
  scale_color_manual(values=palette)
dev.off()

####### Genes expression ######

# fpm_exp <- fpm(ddsHTSeq)

# gene_selected <- "SNAI1"
# count_norm <- fpm_exp
# data_boxplot <- t(count_norm[which(rownames(count_norm)==gene_selected),])
# data_boxplot <- as.data.frame(as.numeric(data_boxplot))

# colnames(data_boxplot) <- gene_selected
# data_boxplot$sample <- sampleTable[,1]
# data_boxplot$condition <- sampleTable[,3]

# outPlot = c(paste(OUTDIR,"/gene_SNAI1_expression.pdf",sep=""))
# pdf(outPlot,width=5.9,height=5.1)
# ggplot(data_boxplot, aes(x = condition, y = get(gene_selected), fill=condition)) +
#   geom_boxplot(alpha=0.5) +
#   theme_minimal() +
#   guides(size="none") +
#   xlab("Cells") +
#   ylab(expression("Normalized expression of "~italic(SNAI1)~" (FPM)")) +
#   guides(fill="none") +
#   scale_fill_manual(values=palette)
# dev.off()

# gene_selected <- "SNAI2"
# count_norm <- fpm_exp
# data_boxplot <- t(count_norm[which(rownames(count_norm)==gene_selected),])
# data_boxplot <- as.data.frame(as.numeric(data_boxplot))

# colnames(data_boxplot) <- gene_selected
# data_boxplot$sample <- sampleTable[,1]
# data_boxplot$condition <- sampleTable[,3]

# outPlot = c(paste(OUTDIR,"/gene_SNAI2_expression.pdf",sep=""))
# pdf(outPlot,width=5.9,height=5.1)
# ggplot(data_boxplot, aes(x = condition, y = get(gene_selected), fill=condition)) +
#   geom_boxplot(alpha=0.5) +
#   theme_minimal() +
#   guides(size="none") +
#   xlab("Cells") +
#   ylab(expression("Normalized expression of "~italic(SNAI2)~" (FPM)")) +
#   guides(fill="none") +
#   scale_fill_manual(values=palette)
# dev.off()

####### Analysis ###########

pw_tests <- list()
pw_tests[[1]] <- c('RH30_APA_14h','RH30_CTL_14h')
pw_tests[[2]] <- c('RH30_D11_14h','RH30_CTL_14h')
pw_tests[[3]] <- c('RH30_APA_48h','RH30_CTL_48h')
pw_tests[[4]] <- c('RH30_D11_48h','RH30_CTL_48h')

# Variables
condCol="condition"
logFCthreshold=LOGFCTHRESHOLD
AdjPValthreshold=FDRTHRESHOLD
GenesInFig=50
bootstrap=FALSE
nboot=30

for (i in 1:length(pw_tests)){
  cond1 <- pw_tests[[i]][1]
  cond2 <- pw_tests[[i]][2]
  comp <- paste(cond1,cond2,sep="__vs__")

  sampleAnnot <- as.data.frame(sampleTable[,3])
  rownames(sampleAnnot) <- sampleTable[,2]
  colnames(sampleAnnot) <- "condition"
  samples <- rownames(sampleAnnot)

  res <- results(ddsHTSeq,contrast=c(condCol,cond1,cond2),independentFiltering=T)

  ## Expression matrices : normalised and transformed
  exprDat <- matrix
  exprDatT <- mrld

  res$meanInComp <- rowMeans(exprDat[,sampleAnnot[,condCol]%in% c(cond1,cond2)])
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
  DE$symbol <- rownames(DE)
  write.table(DE[which(DE$DE!="NONE"),], file=paste(c(OUTDIR,"/",comp,"_","DE_genes.tsv"),collapse=""), sep="\t", row.names=F, quote=F)
  write.table(DE, file=paste(c(OUTDIR,"/",comp,"_","DE_all_genes.tsv"),collapse=""), sep="\t", row.names=F, quote=F)


  DE$DE_label <- NA
  DE$DE_label[which(DE$DE != "NONE")] <- DE[which(DE$DE != "NONE"),"symbol"]
  #list_genes <- c("KCNN1","KCNN3")
  #DE$DE_label[which(rownames(DE) %in% list_genes)] <- rownames(DE[which(rownames(DE) %in% list_genes),])
  DE_palette = c("firebrick3", "gray60", "palegreen3")
  names(DE_palette) = c("DOWN","NONE","UP")

  #Volcano-plot
  pdf(paste(c(OUTDIR,"/",comp,"_Volcano_with_names.pdf"),collapse=""),width=5.9,height=5.9)
    print(ggplot(data=DE, aes(x=log2FoldChange, y=-log10(padj), col=DE, label=DE_label)) + 
    geom_point() +
    theme_minimal() +
    geom_text_repel(color="black", min.segment.length=0, nudge_x=0.5) +
    scale_color_manual(values=DE_palette))
  dev.off()


  pdf(paste(c(OUTDIR,"/",comp,"_Volcano_without_names.pdf"),collapse=""),width=5.9,height=5.9)
    print(ggplot(data=DE, aes(x=log2FoldChange, y=-log10(padj), col=DE, label=DE_label)) + 
    geom_point() +
    theme_minimal() +
    scale_color_manual(values=DE_palette))
  dev.off()

}

# exprDatT <- mrld
# go_gmt <- read.table(paste(c(OUTDIR,"/../../GO_human_2025_genes_only.gmt"),collapse=""),sep=";",h=F,quote="")
# get_gene_and_ontology <- function(line){
#     line <- unlist(strsplit(line,"\t"))
#     this_data <- data.frame("gene"=line[3:length(line)],"vector"=rep(line[1],length(line)-2))
#     return(this_data)
# }
# data_ontologies <- do.call("rbind",apply(go_gmt, 1, function(x) get_gene_and_ontology(x)))
# data_ontologies <- data_ontologies[which(data_ontologies$gene!=""),]

# these_genes <- unique(data_ontologies[which(data_ontologies[,2] == "cell migration"),1])
# these_genes <- these_genes[which(these_genes %in% rownames(exprDatT))]
# sampleAnnot <- as.data.frame(sampleTable[,3])
# rownames(sampleAnnot) <- sampleTable[,2]
# colnames(sampleAnnot) <- "Condition"
# samples <- rownames(sampleAnnot)
# condCol="Condition"
# exprDat <- matrix
# exprDatT <- mrld
# annot=sampleAnnot[condCol]
# annot[,condCol] <- as.factor(annot[,condCol])
# colTopAnnot <- vector("list", ncol(annot))
# names(colTopAnnot) <- "Condition"
# i<-1
# for(col in colnames(annot)){
#   colTopAnnot[[col]]<-palette
#   i<-i+1
# }
# ha<-HeatmapAnnotation(df=annot,col = colTopAnnot,name = condCol)
# DEgenes.names=these_genes 
# sampleHt<-colnames(exprDatT)
# haByComp<-ha
# exprDE=exprDatT[DEgenes.names,sampleHt,drop=FALSE]
# exprDE.scaled=na.omit(rowScale(exprDE,center=T,scaled=T))
# hclustGeneDE<-unsupervisedClustering(exprDE.scaled,transpose = F,nboot=nboot,bootstrap = F)
# #hclustSampleDE<-unsupervisedClustering(exprDE.scaled,transpose = T,nboot=nboot,bootstrap = bootTemp,method.dist = "euclidean")
# kmeans_samples <- c(1,1,1,2,2,2,3,3,3,4,4,4)
# names(kmeans_samples) <- sampleTable[,1]
# split <-factor(kmeans_samples, levels=c(kmeans_samples[1],kmeans_samples[7],kmeans_samples[4],kmeans_samples[10]))
# rowScaledExpr<-rowScale(exprDatT,center=TRUE,scale=TRUE)
# quantile.expr<-quantile(unlist(rowScaledExpr),seq(0,1,.01),na.rm = T)
# colHA=colorRamp2(c(quantile.expr[2],0,quantile.expr[100]),c("#189fc8","white","red"))
# row_font_size <- autoGparFontSizeMatrix(nrow(exprDE.scaled))
# pdf(paste(OUTDIR,"Heatmap_GO_cell_migration.pdf",sep="/"),width=7,height=9)
# draw(Heatmap(exprDE.scaled,top_annotation = haByComp, row_names_gp = row_font_size,
#               cluster_rows = hclustGeneDE,col = colHA,name="Z-score",cluster_columns = F,
#               column_names_gp = autoGparFontSizeMatrix(ncol(exprDE.scaled)),
#               column_split=split,cluster_column_slices = F,column_title = NULL),
#               merge_legend = TRUE)
# dev.off()


