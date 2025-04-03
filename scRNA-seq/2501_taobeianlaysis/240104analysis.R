#250104analysis
#comparison of MOSE & 13D TAOBEI


#QC
suppressMessages(library(Seurat))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(gdata))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(harmony))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(qs))
suppressPackageStartupMessages(library(RColorBrewer))
suppressMessages(library(patchwork))

select = dplyr::select
filter = dplyr::filter

#read in data
getwd()
list.files('/home/cwt/Data/Private/TAOBEI/')
xdata = Read10X(data.dir = '/home/cwt/Data/Private/TAOBEI/',gene.column = 2)

wholedata = CreateSeuratObject(counts = xdata,project='13D')
table(wholedata$orig.ident)

wholedata$orig.ident = 'TAOBEI'

mito.genes <- grep(pattern = "^MT-", x = rownames(wholedata@assays[["RNA"]]), value = TRUE)
wholedata$percent.mito <- Matrix::colSums(wholedata@assays[["RNA"]]$counts[mito.genes, ])/Matrix::colSums(wholedata@assays[["RNA"]]$counts)
ribo.genes <- grep(pattern = "^RPS|^RPL", x = rownames(wholedata@assays[["RNA"]]), value = TRUE)
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
ggsave(qcp,path = "/home/cwt/Project/MOSE/0104analysis/",filename = "QCfigure.pdf",width = 4.66,height = 13.5)

#Filter
#wholedata_cutoff=subset(wholedata,subset=scrublet_doublet == "False")###remove doublets
#dim(wholedata_cutoff)
#dim(wholedata)

wholedata_cutoff <- subset(x = wholedata, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA > 500 & nCount_RNA < 50000 & percent.mito >  -Inf & percent.mito < 0.3 & percent.ribo < 0.3)
dim(wholedata_cutoff)

#####filter genes 保留在所有细胞中count>=10的基因
counts <- GetAssayData(object =wholedata_cutoff, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
wholedata_cutoff <- CreateSeuratObject(filtered_counts, meta.data = wholedata_cutoff@meta.data)

qsave(wholedata_cutoff,'0104analysis/TAOBEI_cutoff.qs')


##
wholedata = wholedata_cutoff
keep(wholedata,sure = T)

wholedata = wholedata %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000)
VGENES = VariableFeatures(wholedata)
VGENES=setdiff(VGENES,VGENES[grep("^MT|^RPL|^RPS",VGENES)])
wholedata = ScaleData(wholedata,features = VGENES,vars.to.regress = c("percent.mito","nCount_RNA","nFeature_RNA","percent.ribo")) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>% FindNeighbors(dims = 1:30)
wholedata = wholedata %>% RunUMAP(dims = 1:30)

table(wholedata$orig.ident)

for(i in seq(0.1,0.5,by = 0.1)){
  wholedata = FindClusters(wholedata,resolution = i,algorithm = 1)
}
DimPlot(wholedata,label = T)
FeaturePlot(wholedata,features = c('PAX8','EPCAM','CD24'))

DotPlot(wholedata,features = c('PAX8','EPCAM','CD24','PTPRC','DCN','COL1A1','MSLN'))

wholedata.marker = wt_downsamplefindmarker_human(wholedata)
writexl::write_xlsx(wholedata.marker,'taobei.marker.xlsx')

wholedata = subset(wholedata,idents = c('6','8'),invert = T)

wholedata = wholedata %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000)
VGENES = VariableFeatures(wholedata)
VGENES=setdiff(VGENES,VGENES[grep("^MT|^RPL|^RPS",VGENES)])
wholedata = ScaleData(wholedata,features = VGENES,vars.to.regress = c("percent.mito","nCount_RNA","nFeature_RNA","percent.ribo")) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>% FindNeighbors(dims = 1:30)
wholedata = wholedata %>% RunUMAP(dims = 1:30)

for(i in seq(0.1,0.5,by = 0.1)){
  wholedata = FindClusters(wholedata,resolution = i,algorithm = 1)
}
DimPlot(wholedata,group.by = 'RNA_snn_res.0.3')
wholedata@active.ident = wholedata$RNA_snn_res.0.3

whole.marker = wt_downsamplefindmarker_human(wholedata)
writexl::write_xlsx(whole.marker,'taobei.marker.xlsx')

wholedata = RenameIdents(wholedata,`0` = 'tb_epi0',`1` = 'tb_epi1',`2` = 'tb_epi2',
                         `3` = 'tb_epi3',`4` = 'tb_epi4'
                           )

FeaturePlot(wholedata,features = c('PAX8','EPCAM','CD24','WT1','CK7'))
DotPlot(wholedata,features = c('PAX8','EPCAM','CD24','WT1','CK7'))
DimPlot(wholedata)

getwd()
setwd( "/home/cwt/Project/MOSE/taobei_analysis")

qsave(wholedata,file= 'wholedata_ident.qs')

wholedata@active.ident = factor(wholedata@active.ident,levels = c(
  'tb_epi0','tb_epi1','tb_epi2','tb_epi3','tb_epi4'
))
wholedata$dcelltype = wholedata@active.ident

epi.color = c('#F8B739','#E44985','#9DE5FF','#ACA8FF','#73777B')
names(epi.color) = levels(wholedata@active.ident)
DimPlot(wholedata,cols = epi.color,pt.size = 0.001)
ggsave('tb_dimplot.pdf',width = 5.1,height = 3.4)

#marker gene 展示##
library(scibet)

expr_matrix <- GetAssayData(wholedata, assay = "RNA", slot = "data")
expr_df <- as.data.frame(as.matrix(t(expr_matrix)))
expr_df$label <- wholedata$dcelltype

etest_gene <- SelectGene(expr_df, k = 50)
etest_gene

remove_genes <- c("CCNA2", "KIF2C", "TPX2", "CEP55", "RRM2", "ASPM", "TYMS", "HMGB2", "HMGN2", "HMGB3", "CENPF")
filtered_genes <- setdiff(etest_gene, remove_genes)

filtered_genes = c(filtered_genes,'TP63','SOX21','FABP5','FZD7','LEF1')

Marker_heatmap(expr_df, filtered_genes)
ggsave('tb.marker.pdf',width = 12,height = 3.4)


aa

``````````````````````GO annotation of the clusters````````````````````````````````````````
#Go annotation
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

DE = dplyr::filter(whole.marker[[2]],avg_log2FC > 0.585 & p_val_adj <= 0.05)

TC_DE = lapply(c('tb_epi0','tb_epi1','tb_epi2','tb_epi3','tb_epi4'),function(celltype){
  dt = dplyr::filter(DE,cluster == celltype)
  return(dt)
})

wt_enrich = function(genelist){
  enrichGO(gene = rownames(genelist),keyType = 'SYMBOL', #'ENTREZID','ENSEMBL'
           OrgDb = org.Hs.eg.db,ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 0.01,qvalueCutoff = 0.05,readable = FALSE)
}
ego1 = wt_enrich(TC_DE[[1]])
ego2 = wt_enrich(TC_DE[[2]])
ego3 = wt_enrich(TC_DE[[3]])
ego4 = wt_enrich(TC_DE[[4]])
ego5 = wt_enrich(TC_DE[[5]])

dotplot(ego5)

extractego = function(ego,name){
  sub = ego@result
  sub$Cluster = name
  return(sub)
}

ego1_df = extractego(ego1,name = 'tb_epi0')
ego2_df = extractego(ego2,name = 'tb_epi1')
ego3_df = extractego(ego3,name = 'tb_epi2')
ego4_df = extractego(ego4,name = 'tb_epi3')
ego5_df = extractego(ego5,name = 'tb_epi4')

ego1_df$P1 = -log10(ego1_df$p.adjust)
ego2_df$P1 = -log10(ego2_df$p.adjust)
ego3_df$P1 = -log10(ego3_df$p.adjust)
ego4_df$P1 = -log10(ego4_df$p.adjust)
ego5_df$P1 = -log10(ego5_df$p.adjust)

ego_result = dplyr::bind_rows(ego1_df,ego2_df,ego3_df,ego4_df,ego5_df)

writexl::write_xlsx(ego_result,'ego_bp_epi_result.xlsx')
write.csv(ego_result,'ego_bp_epi_result.csv')

##select pathway and make figure##
{
  
  epi.color
  p_epi0 = ggplot(subset(ego1_df,ID %in% c('GO:0071347','GO:0045766','GO:0045765','GO:0050673','GO:0010631')), aes(x = reorder(Description, -P1), y = -log10(p.adjust))) +
    geom_bar(stat = "identity", fill = "#F8B739") +
    coord_flip() +
    theme_pubr() +
    labs(title = 'epi0',
         y = "-log10(p.adjust)",
         x = "Description") +
    geom_text(aes(label = Description), hjust = 1, color = "black", size = 4, position = position_nudge(x = 0))+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0))
  
  p_epi1 = ggplot(subset(ego2_df,ID %in% c('GO:0001666','GO:0061621','GO:0046939','GO:0006096','GO:0070265')), aes(x = reorder(Description, -P1), y = -log10(p.adjust))) +
    geom_bar(stat = "identity", fill = "#E44985") +
    coord_flip() +
    theme_pubr() +
    labs(title = 'epi1',
         y = "-log10(p.adjust)",
         x = "Description") +
    geom_text(aes(label = Description), hjust = 1, color = "black", size = 4, position = position_nudge(x = 0))+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0))
  
  p_epi2 = ggplot(subset(ego3_df,ID %in% c('GO:0007160','GO:0034329','GO:0038128','GO:0051403','GO:0022604')), aes(x = reorder(Description, -P1), y = -log10(p.adjust))) +
    geom_bar(stat = "identity", fill = "#9DE5FF") +
    coord_flip() +
    theme_pubr() +
    labs(title = 'epi2',
         y = "-log10(p.adjust)",
         x = "Description") +
    geom_text(aes(label = Description), hjust = 1, color = "black", size = 4, position = position_nudge(x = 0))+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0))
  
  p_epi3 = ggplot(subset(ego4_df,ID %in% c('GO:0061564','GO:0048762','GO:0016055','GO:0198738','GO:0048863')), aes(x = reorder(Description, -P1), y = -log10(p.adjust))) +
    geom_bar(stat = "identity", fill = "#ACA8FF") +
    coord_flip() +
    theme_pubr() +
    labs(title = 'epi3',
         y = "-log10(p.adjust)",
         x = "Description") +
    geom_text(aes(label = Description), hjust = 1, color = "black", size = 4, position = position_nudge(x = 0))+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0))
  
  
  p_epi4 = ggplot(subset(ego5_df,ID %in% c('GO:0007059','GO:0000280','GO:0140014','GO:0000819','GO:0006260')), aes(x = reorder(Description, -P1), y = -log10(p.adjust))) +
    geom_bar(stat = "identity", fill = "#73777B") +
    coord_flip() +
    theme_pubr() +
    labs(title = 'epi4',
         y = "-log10(p.adjust)",
         x = "Description") +
    geom_text(aes(label = Description), hjust = 1, color = "black", size = 4, position = position_nudge(x = 0))+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0))
  
  cowplot::plot_grid(p_epi0,p_epi1,p_epi2,p_epi3,p_epi4,align = "h",ncol = 5)
  ggsave('epi_go_bp.pdf',width= 14.3,height = 2.8)
}

``````````````````````Cluster similarity with MEPP clusters````````````````````````````````````````
library(homologene)
library(ggplot2)
library(GeneOverlap)

getwd()
mepp = qread("/home/cwt/Project/MOSE/1_clustering/epi_1112.qs")
DimPlot(mepp)

MEPP.marker = wt_downsamplefindmarker_mouse(mepp)
genelist = MEPP.marker[[2]]$gene
genetrans = homologene(genelist,inTax = 10090,outTax = 9606)
mepp.marker = merge(MEPP.marker[[2]],genetrans,by.x = 'gene',by.y='10090')
writexl::write_xlsx(mepp.marker,"/home/cwt/Project/MOSE/1_clustering/mepp_epi_marker.xlsx")

#comparison
taobei.list = filter_data(whole.marker[[2]],clustername = unique(whole.marker[[2]]$cluster))
mepp.list = filter_data_m(mepp.marker,clustername = unique(mepp.marker$cluster))
mepp.list1 = list(mepp.list$epi0,mepp.list$epi1,mepp.list$epi2,mepp.list$epi3)
names(mepp.list1) = c('epi0','epi1','epi2','epi3')
p = wt_cluster_comparison(taobei.list,mepp.list1)
p
pdf('/home/cwt/Project/MOSE/taobei_analysis/mepp_comparison_new.pdf',width = 2.5,height = 2.5)
p
dev.off()

##violin plot of mepp and taobei
#the feature of mepp in taobei
DimPlot(wholedata)
{
  epi1.marker = mepp.marker %>% dplyr::filter(cluster == 'epi1') %>%
    top_n(50,wt = avg_log2FC) %>% pull('9606')
  epi1.list = list(epi1.marker)
  names(epi1.list) = 'epi1.feature'
  
  wholedata = AddModuleScore(wholedata,features = epi1.list,name= 'epi1')
  VlnPlot(wholedata,features = c('epi11'))
  colnames(wholedata@meta.data)
  
  #epi0
  epi0.marker = mepp.marker %>% dplyr::filter(cluster == 'epi0') %>%
    top_n(50,wt = avg_log2FC) %>% pull('9606')
  epi0.list = list(epi0.marker)
  wholedata = AddModuleScore(wholedata,features = epi0.list,name= 'epi0')

  VlnPlot(wholedata,features = c('epi01'))
  
  #epi2 and #epi3
  epi2.marker = mepp.marker %>% dplyr::filter(cluster == 'epi2') %>%
    top_n(50,wt = avg_log2FC) %>% pull('9606')
  epi2.list = list(epi2.marker)
  wholedata = AddModuleScore(wholedata,features = epi2.list,name= 'epi2')
  
  epi3.marker = mepp.marker %>% dplyr::filter(cluster == 'epi3') %>%
    top_n(50,wt = avg_log2FC) %>% pull('9606')
  epi3.list = list(epi3.marker)
  wholedata = AddModuleScore(wholedata,features = epi3.list,name= 'epi3')
  
  zheng.color = c('#65799B','#E23E57','#83CC61','#2FC5CC','#F58B54')
  VlnPlot(wholedata,features= c('epi01','epi11','epi21','epi31'),pt.size = 0, col = zheng.color, ncol = 2)
  ggsave(filename = 'module_vln.pdf',width= 3.3,height= 4.3)
}


#comparison with others
#spectrum comparison
spectrum_tumor = readxl::read_xlsx('/home/cwt/Genelist/2022_Nature_cancer_marker_.xlsx')
colnames(spectrum_tumor)[2] = 'cluster'

clustername_spectrum = unique(spectrum_tumor$cluster)
spectrum_tumor_list = filter_data(spectrum_tumor,clustername  = clustername_spectrum)

p1 = wt_cluster_comparison(spectrum_tumor_list,taobei.list)
p1

#chai comparison
chai = read.csv('/home/cwt/Genelist/EPI_each_cluster_allDEGs_toChen.csv')

cluster.chai = c('0','1','2','3','4','5','6','7','8','9','10','11')

chai.list = list()
for (col_name in cluster.chai){
  column = chai %>% dplyr::filter(cluster == col_name) %>% select('gene') %>% pull()
  chai.list[[col_name]] = column
}
chai.list
names(chai.list) = c('Epi0','Epi1','Epi2','Epi3','Epi4','Epi5','Epi6','Epi7','Epi8','Epi9','Epi10','Epi11')

p.chai = wt_cluster_comparison(chai.list,taobei.list)
p.chai

cowplot::plot_grid(plotlist = list(p1,p.chai))
ggsave('tabei_other_compare.pdf',width = 8.3,height = 3.6)




``````````````````````MEPP comparison with ID8````````````````````````````````````````
ppk = qread('/home/cwt/Project/MOSE/250207/wholedata_ident.qs')
mepp
DimPlot(mepp)

colnames(ppk@meta.data)
table(ppk$group)
id8wt = subset(ppk,subset = group == 'ID8-wt')
DimPlot(id8wt)
colnames(id8wt@meta.data)
table(id8wt$bigclass)

id8wt.tumor.seu = subset(id8wt, subset = bigclass == 'Tumor')

id8wt.tumor.seu$Group = 'ID8'
mepp$Group = 'MEPP'

combine.seu = merge(id8wt.tumor.seu,mepp)
table(combine.seu$Group)
combine.seu = JoinLayers(combine.seu)

combine.seu = wt_firststep(combine.seu)

##原位灶之间比较
table(combine.seu$group)
deg = FindMarkers(combine.seu,ident.1 = 'MOSE-OV',ident.2 = 'ID8-wt',group.by = 'group',
                  min.pct=0.1,logfc.threshold = 0,pseudocount.use = F,only.pos = F)
deg_plot = deg %>% mutate(significant = case_when(
  avg_log2FC >= 0.585 & p_val_adj < 0.05 ~ "up",
  avg_log2FC <= -0.585 & p_val_adj < 0.05 ~ "down",
  TRUE ~ "stable"
))

data_clean <- deg_plot %>%
  mutate(
    avg_log2FC = ifelse(is.infinite(avg_log2FC), NA, avg_log2FC),  # Remove infinite values
    p_val_adj = ifelse(p_val_adj == 0, 1e-300, p_val_adj),  # Replace zero p-values
    log_pval = -log10(p_val_adj)  # Calculate -log10(p-value)
  ) %>%
  drop_na(avg_log2FC)


wt_volcano(data_clean)
ggsave('meppvsid8.pdf',width = 4.5,height = 4.0)


#GSEA analysis
library(msigdbr)
library(fgsea)
list.files('/home/cwt/Genelist')
list.files('/home/cwt/Genelist/msigdb_mm')
HM = gmtPathways('/home/cwt/Genelist/mh.all.v2024.1.Mm.symbols.gmt')

up.mepp =  data_clean %>% rownames_to_column(var = "gene") %>% 
  dplyr::filter(avg_log2FC >0) %>% arrange(desc(avg_log2FC)) %>% dplyr::select(gene,avg_log2FC) %>% deframe()

up.mepp.fgsea = fgsea(HM,stats = up.mepp)


mm.reactome = gmtPathways('/home/cwt/Genelist/msigdb_mm/m2.cp.reactome.v2024.1.Mm.symbols.gmt')
up.mepp.fgsea.reactome = fgsea(mm.reactome,stats = up.mepp,scoreType= 'pos')

mm.gobp = gmtPathways('/home/cwt/Genelist/msigdb_mm/m5.go.bp.v2024.1.Mm.symbols.gmt')
up.mepp.fgsea.gobp = fgsea(mm.gobp,stats = up.mepp,scoreType= 'pos')

mm.wiki = gmtPathways('/home/cwt/Genelist/msigdb_mm/m2.cp.wikipathways.v2024.1.Mm.symbols.gmt')
up.mepp.wiki = fgsea(mm.wiki,stats = up.mepp,scoreType= 'pos')

mm.biocarta = gmtPathways('/home/cwt/Genelist/msigdb_mm/m2.cp.biocarta.v2024.1.Mm.symbols.gmt')
up.mepp.bio = fgsea(mm.biocarta,stats = up.mepp,scoreType= 'pos')



p1 = plotEnrichment(mm.gobp[["GOBP_POSITIVE_REGULATION_OF_MAPK_CASCADE"]],
               up.mepp) + labs(title="GOBP:Positive regulation of MAPK")
p2 = plotEnrichment(HM[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]],
                    up.mepp) + labs(title="HM:EMT")
p3 = plotEnrichment(HM[["HALLMARK_GLYCOLYSIS"]],
                    up.mepp) + labs(title="HM:Glycolysis")
p4 = plotEnrichment(mm.wiki[["WP_OXIDATIVE_STRESS_AND_REDOX_PATHWAY"]],
                    up.mepp) + labs(title="WP_OXIDATIVE_STRESS_AND_REDOX_PATHWAY");p4

plot_grid(plotlist = list(p1,p2,p3,p4),ncol =4)
ggsave('meppvsid8_gsea.pdf',width = 13,height = 2.3)

p_bcell = plotEnrichment(mm.gobp[["GOBP_TERPENOID_METABOLIC_PROCESS"]],
               up.mepp) + labs(title="")
p_bcell1 = plotEnrichment(mm.gobp[["GOBP_ISOPRENOID_METABOLIC_PROCESS"]],
               up.mepp) + labs(title="")

plot_grid(plotlist = list(p_bcell,p_bcell1),ncol=1)
ggsave('fake.pdf',width = 2.8,height = 5.4)


##ID8图
id8wt.tumor.seu
DimPlot(id8wt.tumor.seu)
id8wt.tumor.seu = JoinLayers(id8wt.tumor.seu)

id8wt.tumor.seu = id8wt.tumor.seu %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData(vars.to.regress =c("percent.mito","nCount_RNA","nFeature_RNA","percent.ribo")) %>% 
  RunPCA(npcs = 20, verbose = FALSE) %>% RunUMAP(reduction = "pca", dims = 1:20) %>% 
  FindNeighbors(dims = 1:20)

for(i in seq(0.1,0.3,by = 0.1)){
  id8wt.tumor.seu = FindClusters(id8wt.tumor.seu,resolution = i,algorithm = 1)
}
DimPlot(id8wt.tumor.seu,group.by = 'RNA_snn_res.0.1')
id8wt.tumor.seu@active.ident = id8wt.tumor.seu$RNA_snn_res.0.1
id8.marker = wt_downsamplefindmarker_mouse(id8wt.tumor.seu)
#id8wt.tumor.seu = subset(id8wt.tumor.seu,idents = c('4','5','6'),invert = T)
id8wt.tumor.seu = subset(id8wt.tumor.seu,idents = c('4'),invert = T)
id8wt.tumor.seu = RenameIdents(id8wt.tumor.seu,`0`='id8_epi0',`1`='id8_epi1',`2`='id8_epi2',`3`='id8_epi3')
p_id8 = DimPlot(id8wt.tumor.seu,raster = F,pt.size = 0.001)

##ID8内映射MEPP的list
{
  epi1.marker = mepp.marker %>% dplyr::filter(cluster == 'epi1') %>%
    top_n(50,wt = avg_log2FC) %>% pull('gene')
  epi1.list = list(epi1.marker)
  names(epi1.list) = 'epi1.feature'
  id8wt.tumor.seu = AddModuleScore(id8wt.tumor.seu,features = epi1.list,name= 'epi1')
  VlnPlot(id8wt.tumor.seu,features = c('epi11'))

  
  #epi0
  epi0.marker = mepp.marker %>% dplyr::filter(cluster == 'epi0') %>%
    top_n(50,wt = avg_log2FC) %>% pull('gene')
  epi0.list = list(epi0.marker)
  id8wt.tumor.seu = AddModuleScore(id8wt.tumor.seu,features = epi0.list,name= 'epi0')

  #epi2 and #epi3
  epi2.marker = mepp.marker %>% dplyr::filter(cluster == 'epi2') %>%
    top_n(50,wt = avg_log2FC) %>% pull('gene')
  epi2.list = list(epi2.marker)
  id8wt.tumor.seu = AddModuleScore(id8wt.tumor.seu,features = epi2.list,name= 'epi2')
  
  epi3.marker = mepp.marker %>% dplyr::filter(cluster == 'epi3') %>%
    top_n(50,wt = avg_log2FC) %>% pull('gene')
  epi3.list = list(epi3.marker)
  id8wt.tumor.seu = AddModuleScore(id8wt.tumor.seu,features = epi3.list,name= 'epi3')
  
  zheng.color = c('#65799B','#E23E57','#83CC61','#2FC5CC','#F58B54')
  p_vln = VlnPlot(id8wt.tumor.seu,features= c('epi01','epi11','epi21','epi31'),pt.size = 0, col = zheng.color, ncol = 4)
  
  cowplot::plot_grid(plotlist = list(p_id8,p_vln),ncol = 2,rel_widths = c(1,3))
  ggsave(filename = 'id8_mepp.pdf',width= 13.3,height= 2.5)
}

getwd()
qsave(wholedata,'taobei_250212.qs')
qsave(id8wt,'id8wt_whole.qs')
qsave(id8wt.tumor.seu,'id8wt_tumor.qs')
qsave(combine.seu,'meppplusid8.seu.qs')


##试试用combine.seu。但是ID8-wt和MEPP很难有好的分层
combine.seu ##meppy原位+网膜+ID8
DimPlot(combine.seu)
colnames(combine.seu@meta.data)
table(combine.seu$orig.ident)

DimPlot(combine.seu,group.by = 'orig.ident')
DimPlot(combine.seu,group.by = 'RNA_snn_res.0.1')

```````````````````MEPP-OV OM进一步补充`````````````````````````````````````
DimPlot(mepp)
colnames(mepp@meta.data)
table(mepp$group)
FeaturePlot(object = mepp, features = c("Wt1", "Cdh2", "Akt1", "Mapk3"))
VlnPlot(mepp,features = c("Wt1", "Cdh2", "Akt1", "Mapk3"),group.by = "group")







###############FAKE PLOT 250207######
library(msigdbr)
library(fgsea)
list.files('/home/cwt/Genelist/')
HM = gmtPathways('/home/cwt/Genelist/h.all.v2023.1.Hs.symbols.gmt') 

deg = FindMarkers(wholedata,ident.1 = 'tb_epi0',ident.2 = 'tb_epi1',group.by = 'dcelltype',
                  min.pct=0.1,logfc.threshold = 0,pseudocount.use = F,only.pos = F)

deg_plot = deg %>% mutate(significant = case_when(
  avg_log2FC >= 0.585 & p_val_adj < 0.05 ~ "up",
  avg_log2FC <= -0.585 & p_val_adj < 0.05 ~ "down",
  TRUE ~ "stable"
))

wt_volcano(deg_plot)
ggsave('fake_volca.pdf',width = 4.4,height = 3.3)

genelist = deg %>% rownames_to_column(var = "gene") %>% 
  dplyr::filter(avg_log2FC >0) %>% arrange(desc(avg_log2FC)) %>% dplyr::select(gene,avg_log2FC) %>% deframe()

res = fgsea(HM,stats = genelist)

p1 = plotEnrichment(HM[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
               genelist) + labs(title="")
p2 = plotEnrichment(HM[["HALLMARK_IL6_JAK_STAT3_SIGNALING"]],
                    genelist) + labs(title="");p2
p3 = plotEnrichment(HM[["HALLMARK_INFLAMMATORY_RESPONSE"]],
                    genelist) + labs(title="");p3
p4 = plotEnrichment(HM[["HALLMARK_KRAS_SIGNALING_DN"]],
                    genelist) + labs(title="");p4
library(cowplot)
plot_grid(plotlist = list(p1,p2,p3,p4),ncol =2)
ggsave('fake_gsea.pdf',width = 4.7,height = 3.7)


###############FAKE PLOT 250207######
fakedata = subset(wholedata,downsample = 400)
DimPlot(fakedata)

getwd()
dt = qread('/home/cwt/Project/MOSE/taobei_analysis/id8wt_tumor.qs')
FeaturePlot(dt,features = c('Cd276'),order =T)
VlnPlot(dt,features = c('Cd276'))

library(Seurat)
DimPlot(mepp)
FeaturePlot(mepp,features = c('Msln'))
