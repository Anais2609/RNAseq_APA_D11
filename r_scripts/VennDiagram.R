library(ggplot2)
library(VennDiagram)

upregulated_RH4 <- read.table("../DESEQ2/results_all/RH4_siKCNN3__vs__RH4_siCT_DE_all_genes.tsv",sep="\t",h=T)
upregulated_RH4 <- upregulated_RH4[which(upregulated_RH4$DE=="UP"),9]
upregulated_RH30 <- read.table("../DESEQ2/results_all/RH30_siKCNN3__vs__RH30_siCT_DE_all_genes.tsv",sep="\t",h=T)
upregulated_RH30 <- upregulated_RH30[which(upregulated_RH30$DE=="UP"),9]

downregulated_RH4 <- read.table("../DESEQ2/results_all/RH4_siKCNN3__vs__RH4_siCT_DE_all_genes.tsv",sep="\t",h=T)
downregulated_RH4 <- downregulated_RH4[which(downregulated_RH4$DE=="DOWN"),9]
downregulated_RH30 <- read.table("../DESEQ2/results_all/RH30_siKCNN3__vs__RH30_siCT_DE_all_genes.tsv",sep="\t",h=T)
downregulated_RH30 <- downregulated_RH30[which(downregulated_RH30$DE=="DOWN"),9]

my_plot <- venn.diagram(
  x = list(upregulated_RH4,upregulated_RH30),
  category.names = c("RH4","RH30"),
  filename = NULL,
  output = TRUE ,
          imagetype="svg" ,
          height = 480 , 
          width = 480 , 
          resolution = 300,
          compression = "lzw",
          lwd = 1,
          cex = 2,
          fontfamily = "sans",
          col=c("#c70039","#f28824"),
          fill = c(alpha("#c70039",0.3), alpha("#f28824",0.3)),
          cat.cex = 2.2,
          cat.fontface = "bold",
          cat.fontfamily = "sans",
          cat.col = c("#c70039","#f28824"),
          cat.default.pos = "outer",
          cat.pos = c(-35,35),
          cat.dist = c(0.05,0.05)
        )
ggsave(my_plot, file='../DESEQ2/results_all/venn_upregulated.svg', device = "svg")
write.table(upregulated_RH4[which((upregulated_RH4 %in% upregulated_RH30))],'../DESEQ2/results_all/common_upregulated.txt',sep="\t",quote=F,row.names=F)

my_plot <- venn.diagram(
  x = list(downregulated_RH4,downregulated_RH30),
  category.names = c("RH4","RH30"),
  filename = NULL,
  output = TRUE ,
          imagetype="svg" ,
          height = 480 , 
          width = 480 , 
          resolution = 300,
          compression = "lzw",
          lwd = 1,
          cex = 2,
          fontfamily = "sans",
          col=c("#c70039","#f28824"),
          fill = c(alpha("#c70039",0.3), alpha("#f28824",0.3)),
          cat.cex = 2.2,
          cat.fontface = "bold",
          cat.fontfamily = "sans",
          cat.col = c("#c70039","#f28824"),
          cat.default.pos = "outer",
          cat.pos = c(-35,35),
          cat.dist = c(0.05,0.05)
        )
ggsave(my_plot, file='../DESEQ2/results_all/venn_downregulated.svg', device = "svg")
write.table(downregulated_RH4[which((downregulated_RH4 %in% downregulated_RH30))],'../DESEQ2/results_all/common_downregulated.txt',sep="\t",quote=F,row.names=F)

