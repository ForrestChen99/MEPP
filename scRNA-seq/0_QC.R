getwd()
options(stringsAsFactors = F)
set.seed(42)

suppressMessages(library(Seurat))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(gdata))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(harmony))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(qs))
suppressPackageStartupMessages(library(RColorBrewer))
suppressMessages(library(patchwork))


###########Data import###########

setwd("/work/cwt/MOSE/Data")
datalist = list.dirs(recursive = F);datalist


wholedata_list = lapply(datalist,function(x){
  xdata = Read10X(data.dir = paste(x,"/",sep = ""),gene.column = 2) #gene.column define gene.column in features.tsv ##华大样本column = 1，BD和10x gene.column=2？
  xf = CreateSeuratObject(counts = xdata,project=x)
  double_score = read.delim(paste("/work/cwt/MOSE/scrublet_result",gsub("[.]","",x),"doublet.txt",sep = ""),sep = ",")
  xf$scrublet_doublet_score = double_score$doublet_scores
  xf$scrublet_doublet = double_score$predicted_doublets
  print(x)
  return(xf)
}) 

wholedata = wholedata_list[[1]]
wholedata = merge(x = wholedata,y = wholedata_list[2:length(wholedata_list)],
                  merge.data = T,project = "MOSE")
table(wholedata$orig.ident)

wholedata$orig.ident = gsub("[./]","",wholedata$orig.ident)

setwd('/work/cwt/MOSE/')
qsave(wholedata,'0_QC/wholedata_raw.qs')



wholedata = JoinLayers(wholedata)


mito.genes <- grep(pattern = "^mt-", x = rownames(wholedata@assays[["RNA"]]), value = TRUE)
wholedata$percent.mito <- Matrix::colSums(wholedata@assays[["RNA"]]$counts[mito.genes, ])/Matrix::colSums(wholedata@assays[["RNA"]]$counts)
ribo.genes <- grep(pattern = "^Rps|^Rpl", x = rownames(wholedata@assays[["RNA"]]), value = TRUE)
wholedata$percent.ribo <- Matrix::colSums(wholedata@assays[["RNA"]]$counts[ribo.genes, ])/Matrix::colSums(wholedata@assays[["RNA"]]$counts)

{
  
  qcp1 = VlnPlot(object = wholedata, features = c("nFeature_RNA"),group.by="orig.ident",pt.size = 0)+
    geom_hline(yintercept = seq(2000,9000,by=1000),linetype="dashed",color="red")
  
  qcp2 = VlnPlot(object = wholedata, features = c("nCount_RNA"),group.by="orig.ident",pt.size = 0)+
    geom_hline(yintercept = seq(2000,9000,by=1000),linetype="dashed",color="red")
  
  qcp3 = VlnPlot(object =wholedata, features = c("percent.mito"),group.by="orig.ident",pt.size = 0)+
    geom_hline(yintercept = seq(0.05,0.4,by=0.05),linetype="dashed")
  
  qcp4 = VlnPlot(object =wholedata, features = c("percent.ribo"),group.by="orig.ident",pt.size = 0)+
    geom_hline(yintercept = seq(0.05,0.4,by=0.05),linetype="dashed")
  
  qcp5 = FeatureScatter(object = wholedata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by="orig.ident")
  
  qcp = qcp1+qcp2+qcp3+qcp4+qcp5
  qcp
}
ggsave(qcp,path = "/work/cwt/MOSE/0_QC/",filename = "QCfigure.pdf",width = 4.66,height = 13.5)

#Filter
wholedata_cutoff=subset(wholedata,subset=scrublet_doublet == "False")###remove doublets
dim(wholedata_cutoff)
dim(wholedata)

wholedata_cutoff <- subset(x = wholedata_cutoff, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mito >  -Inf & percent.mito < 0.3 & percent.ribo < 0.3)
dim(wholedata_cutoff)

#####filter genes 保留在所有细胞中count>=10的基因
counts <- GetAssayData(object =wholedata_cutoff, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
wholedata_cutoff <- CreateSeuratObject(filtered_counts, meta.data = wholedata_cutoff@meta.data)
#####
num_cell=cbind(table(wholedata$orig.ident),table(wholedata_cutoff$orig.ident))
colnames(num_cell)=c("raw","filtered")
write.csv(num_cell,file="0_QC/cell_num_QC_statics.csv")

qsave(wholedata_cutoff,'0_QC/wholedata_cutoff.qs')

wholedata = wholedata_cutoff
keep(wholedata,sure = T)








