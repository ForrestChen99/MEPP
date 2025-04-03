{
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
  suppressMessages(library(ggpubr))
  suppressMessages(library(tidyr))
  suppressMessages(library(pheatmap))
  suppressMessages(library(ComplexHeatmap))
  library(circlize)
  library(Seurat)
  library(presto)
  library(data.table)
  library(tidyverse)
}
packageVersion('Seurat')
remotes::install_version("Matrix", version = "1.6-4")
install.packages('SeuratObject')
install.packages('Seurat')

devtools::install_github('immunogenomics/presto')


getwd()
wholedata = qs::qread('/home/cwt/Project/MOSE/0_QC/wholedata_cutoff.qs')
wholedata
table(wholedata$orig.ident)
wholedata$group = wholedata$orig.ident
wholedata$group = factor(wholedata$group,levels = c('MOSE-OV','MOSE-OM'))
table(wholedata$group)

#
wholedata = wholedata %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000)
VGENES = VariableFeatures(wholedata)
VGENES=setdiff(VGENES,VGENES[grep("^mt|^Rpl|^Rps",VGENES)])

wholedata = ScaleData(wholedata,features = VGENES,vars.to.regress = c("percent.mito","nCount_RNA","nFeature_RNA","percent.ribo")) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>% FindNeighbors(dims = 1:30)

for(i in seq(0.1,0.5,by = 0.1)){
  wholedata = FindClusters(wholedata,resolution = i,algorithm = 1)
}
wholedata = wholedata %>% RunUMAP(dims = 1:30)

wholedata@active.ident = wholedata$RNA_snn_res.0.1
DotPlot(wholedata,features = c('Wt1','Hspd1','Ptprc','Cd14','Cd3e','Postn','Dcn','Prlr','Epcam','Hba-a1','Acta1','Mki67'),
        cols = c('grey','red'))
DimPlot(wholedata,group.by = 'orig.ident')
DimPlot(wholedata,label = T)
#

{
  wholedata_harmony = RunHarmony(wholedata,group.by.vars = 'orig.ident',assay.use = 'RNA',reduction = 'pca',project.dim = F)
  harmony_embeddings = Embeddings(wholedata_harmony,'harmony')
  wholedata_harmony = wholedata_harmony %>% RunUMAP(reduction = 'harmony',dims = 1:30) %>%
    FindNeighbors(reduction = 'harmony',dims = 1:30)
  
  for(i in seq(0.1,0.5,by = 0.1)){
    wholedata_harmony = FindClusters(wholedata_harmony,resolution = i,algorithm = 1)
  }
  
  DimPlot(wholedata_harmony,group.by = 'RNA_snn_res.0.1',label = T)
  DotPlot(wholedata_harmony,features = c('Wt1','Hspd1','Ptprc','Cd14','Cd3e','Postn','Dcn','Prlr','Epcam','Hba-a1','Acta1','Mki67'),
          cols = c('grey','red'))
  DotPlot(wholedata_harmony,features = c('Mki67'),
          cols = c('grey','red'))
  DimPlot(wholedata_harmony,group.by = 'group')
  wholedata_harmony@active.ident = wholedata_harmony$RNA_snn_res.0.1
  
  FeaturePlot(wholedata_harmony,features = c('nFeature_RNA','scrublet_doublet_score'))
  
  whole.marker = wt_downsamplefindmarker_mouse(wholedata_harmony)
  
  
  
  
  wholedata_ident = RenameIdents(wholedata_harmony,`0` = 'Epithelial',`1` = 'Macrophage',`2` = 'TNK',`3` = 'B',`4` = 'RBC',`5` = 'Fibroblast',`6` = 'Granulosa',
                                 `7` = 'Neutrophil',`8` = 'cDC1',`9` = 'Mast')
  wholedata_ident = subset(wholedata_ident,idents = '10',invert = T)
  DimPlot(wholedata_ident)
   
  
  wholedata_ident@active.ident = factor(wholedata_ident@active.ident,levels = c(
    'Epithelial','Granulosa','Macrophage','cDC1','Neutrophil','Mast','TNK','B','Fibroblast','RBC'
  ))
  wholedata_ident$bigclass = wholedata_ident@active.ident
  qsave(wholedata_harmony,file = '1_clustering/wholedata_ident.qs')
}

bigclass.color = c('#3AA6B9','#FFC15E','#758694','#EB5E60','#78BBE6','#FF895D','#FFB4C2','#5AA897','#9A7E6F','#E7CCCC')
names(bigclass.color) = unique(levels(wholedata_ident$bigclass))

#figures
if(F){
  p1 = DimPlot(wholedata_ident,cols = bigclass.color,label=T);p1
  ggsave('whole.dimplot.pdf',path = '1_clustering/',width = 5,height = 4)
  
  whole.marker = wt_downsamplefindmarker_mouse(wholedata_ident)
  writexl::write_xlsx(whole.marker,'wholemarker.xlsx')
  DotPlot(wholedata_ident,features = c(
    'Wt1','Msln','Krt19',
    'Inha','Amh','Ihh',
    'C1qc','Cd68',
    'Clec9a','Xcr1',
    'S100a8','S100a9','Cxcr2',
    'Tpsab1','Ms4a2','Cd3g','Cd3d','Nkg7',
    'Cd79a','Cd19',
    'Igfbp7','Rhoj','Cdh11',
    'Rhag','Alas2'))+theme(axis.text.x = element_text(angle = 90,hjust=1,vjust = 0.5),plot.title = element_text(hjust = 0.5,size = 8, face = "bold"))+scale_color_distiller(palette = "RdBu")
  ggsave('whole.dotplot.pdf',path = '1_clustering/',width = 8.4,height = 3.6)
  
  DimPlot(wholedata_ident,split.by = 'group',cols = bigclass.color)
  ggsave('whole.split.dim.pdf',width = 7.61,height = 3.6,path = '1_clustering/')
}

#clusterprofiling
if(F){
  whole.meta = wholedata_ident@meta.data
  whole.meta$bigclass = factor(whole.meta$bigclass,levels = c(
    'Epithelial','Granulosa','Macrophage','cDC1','Neutrophil','Mast','TNK','B','Fibroblast','RBC'
  ))
  ClusterFreq <- whole.meta[,c("bigclass","group")] %>% table %>%
    data.frame() %>% set_colnames(c("seurat_clusters","stim","Per"))
  ClusterPer <- ClusterFreq %>% tidyr::spread(seurat_clusters,Per)
  ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
  ClusterPer1 <- ClusterPer %>% tidyr::gather(key=seurat_clusters,value=Per,-stim)
  
  ClusterPer1$seurat_clusters = factor(ClusterPer1$seurat_clusters,levels = c(
    'Epithelial','Granulosa','Macrophage','cDC1','Neutrophil','Mast','TNK','B','Fibroblast','RBC'
  ))
  
  ggplot(ClusterPer1,  aes(x = seurat_clusters, y = Per, color = stim)) +
    geom_point(position = position_dodge(width = 0.5)) +
    theme_pubr() +
    labs(y = 'Proportion among all cells (%)', x = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = 'top') +
    scale_color_manual(values = c("MOSE-OV" = "#92CCE1", "MOSE-OM" = "#F68686"))
  ggsave('whole.dot.proportion.pdf',width = 3.8,height = 4,path = '1_clustering/')
  
  
  plotdata = ClusterPer %>% column_to_rownames(var = 'stim')
  plotdata_scale = t(apply(plotdata,1,scale))
  colnames(plotdata_scale) = colnames(plotdata)
  
  
  anno.data = whole.meta %>% dplyr::select(orig.ident,group) %>% arrange(group)
  anno.data <- anno.data %>% distinct(orig.ident, .keep_all = TRUE)
  rownames(anno.data) = NULL
  anno.data$group = factor(anno.data$group,levels = c("MOSE-OV",'MOSE-OM'))
  
  
  group.color = c('#006989','#E88D67')
  names(group.color) = unique(levels(anno.data$group))
  
  top_annotation = HeatmapAnnotation(
    Group = anno.data$group,
    col = list(
      Group = group.color
    )
  )
  
  ht = Heatmap(t(plotdata_scale),cluster_rows = F,cluster_columns = F,show_row_dend = F,show_column_names = F,
               col = colorRamp2(c(-1, 0, 3), c("#1679AB", "white", "#C80036")),top_annotation = top_annotation)
  pdf('1_clustering/whole.heatmap.proportion.pdf',width = 3.8,height = 4)
  draw(ht)
  dev.off()
  
}

#HALLMARK evaluation of all clusters
{
  library(msigdbr)
  library(fgsea)
  library(GSVA)

  list.files('/home/cwt/Genelist/')
  HM = gmtPathways('/home/cwt/Genelist/mh.all.v2024.1.Mm.symbols.gmt')  
  
  dim(wholedata)
  bigclass = unique(wholedata@meta.data$bigclass)
  
  source('/home/cwt/Genelist/scores.r') ##serua5V5 change data acquistion way
  
  dt = wholedata@assays[['RNA']]$data
  
  
  HMscore <- score_cells(seur=wholedata, names=HM, combine_genes='mean', 
                         groups=NULL, group_stat='mean', cells.use=NULL)
  HMmatrix = as.matrix(HMscore) %>% t
  HMmatrix = apply(HMmatrix, 2, function(x)signif(x,digits = 3))
  colnames(HMmatrix) = rownames(wholedata@meta.data)
  HMseuobj = CreateSeuratObject(counts = HMmatrix,data = HMmatrix,meta.data = wholedata@meta.data)
  
  wholedata@assays$HM = HMseuobj@assays$RNA
  DefaultAssay(wholedata) = 'HM'
  
  table(wholedata$group)
  Diff_HM = lapply(bigclass,function(i){
    sub = subset(wholedata,subset = bigclass == i)
    DEfeatures =FindMarkers(sub,ident.1 = 'MOSE-OM',ident.2 = 'MOSE-OV',group.by = 'group',
                            min.pct=0.1,logfc.threshold = 0,pseudocount.use = F,only.pos = F)
    DEfeatures$terms = rownames(DEfeatures)
    DEfeatures$celltype = i
    return(DEfeatures)
    print(i)
  })
  
  process_diff_combined <- function(Diff_HM, dcelltype, selected_hallmark) {
    # Set names for the list
    names(Diff_HM) <- dcelltype
    
    # Combine the list into a single dataframe
    Diff_combined <- dplyr::bind_rows(Diff_HM)
    rownames(Diff_combined) <- NULL
    
    # Modify the Pathway column
    Diff_combined <- Diff_combined %>%
      mutate(Pathway = gsub("HALLMARK-", "", terms))
    
    # Modify the Pathway names
    Diff_combined$Pathway <- gsub('-', '_', Diff_combined$Pathway)
    
    # Filter and factor the Diff_combined dataframe
    Diff_combined_selected <- dplyr::filter(Diff_combined, Pathway %in% selected_hallmark$Name)
    Diff_combined_selected$Pathway <- factor(Diff_combined_selected$Pathway, levels = selected_hallmark$Name)
    
    return(Diff_combined_selected)
  }
  
  selected_hallmark = readxl::read_xlsx('/home/cwt/Genelist/HALLMARK_CATEGORY.xlsx')

  plot.data = process_diff_combined(Diff_HM,bigclass,selected_hallmark)
  
  ggplot(plot.data,aes(y = celltype,x = Pathway))+
    geom_point(aes(color=avg_log2FC,size=-log10(p_val_adj+10^-100)))+
    scale_color_gradient2(low = "#769FCD", mid = "white", high = "#FC5185", midpoint = 0, 
                          limits = c(-0.5, 1))+
    labs(title = 'HALLMARK:OM vs OV')+
    scale_size_continuous(range=c(1,6),breaks=seq(0,30,10),name="FDR")+
    theme(panel.background=element_rect(color="black",fill=NA),legend.key =  element_blank(),axis.ticks=element_blank(),
          panel.grid=element_blank(),axis.text.x=element_text(size=8,angle=90,vjust = 0,hjust = 1,color="black"),
          axis.text.y=element_text(size=8,color="black"),axis.title=element_blank(),legend.position = "bottom")+
    guides(size = guide_legend(nrow = 2))
  
  ggsave('hallmark_whole_omvsov.pdf',width = 11,height= 5.6)
  
  
}

devtools::install_github("junjunlab/GseaVis")



``````````````````````````Subcluster``````````````````````````````````````````````
wholedata = wholedata_ident
qsave(wholedata,'1_clustering/wholedata_1031.qs')


epi = subset(wholedata,idents = 'Epithelial')
epi = wt_firststep(epi)

DimPlot(epi,group.by = 'RNA_snn_res.0.2')
epi@active.ident = epi$RNA_snn_res.0.2
epi.marker = wt_downsamplefindmarker_mouse(epi)

epi1 = subset(epi,idents = c('4','5','6'),invert = T) #remove other cell types


epi1 <- epi1 %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20)
DimPlot(epi1)

epi.marker = wt_downsamplefindmarker_mouse(epi1)

epi = epi1 %>% RenameIdents(`0` = 'epi0',`1` = 'epi1',`2`= 'epi2',`3` = 'epi3')
epi@active.ident = factor(epi@active.ident,levels = c('epi0','epi1','epi2','epi3'))
epi$dcelltype = epi@active.ident

epi.color = c('#CD8D7A','#A0C49D','#96B6C5','#607274')
names(epi.color) = unique(levels(epi$dcelltype))

DimPlot(epi,label= F,cols = epi.color)
DimPlot(epi,label= F,cols = epi.color,split.by = 'group')

#qsave(epi,file = '1_clustering/epi.qs')
#figures
if(F){
  p1 = DimPlot(epi,label= F,cols = epi.color);p1
  ggsave('epi.dimplot.pdf',path = '1_clustering/',width = 5,height = 4)
  
  epi.marker = wt_downsamplefindmarker_mouse(epi)
  writexl::write_xlsx(epi.marker,'epimarker.xlsx')

  
  DimPlot(epi,label= F,cols = epi.color,split.by = 'group')
  ggsave('epi.split.dim.pdf',width = 7.61,height = 3.6,path = '1_clustering/')
}

#clusterprofiling
if(F){
  epi.meta = epi@meta.data
  epi.meta$dcelltype = factor(epi.meta$dcelltype,levels = c(
    'epi0','epi1','epi2','epi3'
  ))
  ClusterFreq <- epi.meta[,c("dcelltype","group")] %>% table %>%
    data.frame() %>% set_colnames(c("seurat_clusters","stim","Per"))
  ClusterPer <- ClusterFreq %>% tidyr::spread(seurat_clusters,Per)
  ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
  ClusterPer1 <- ClusterPer %>% tidyr::gather(key=seurat_clusters,value=Per,-stim)
  
  ClusterPer1$seurat_clusters = factor(ClusterPer1$seurat_clusters,levels = c(
    'epi0','epi1','epi2','epi3'
  ))
  
  ggplot(ClusterPer1,  aes(x = seurat_clusters, y = Per, color = stim)) +
    geom_point(position = position_dodge(width = 0.5)) +
    theme_pubr() +
    labs(y = 'Proportion among all cells (%)', x = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = 'top') +
    scale_color_manual(values = c("MOSE-OV" = "#92CCE1", "MOSE-OM" = "#F68686"))
  ggsave('epi.dot.proportion.pdf',width = 3.1,height = 3.25,path = '1_clustering/')
}


#ROIE of epi clusters
{
  epi.meta = epi@meta.data
  summary = table(epi.meta[,c('dcelltype','group')])
  roe = as.data.frame(ROIE(summary))
  roe$cluster = rownames(roe)
  rownames(roe) = NULL
  res = data.frame()
  res = rbind(res,roe)
  summary(roe)
  
  heatmapdata = gather(roe,group,value,-c(cluster))
  str(heatmapdata)
  
  heatmapdata$group = factor(heatmapdata$group,levels = c('MOSE-OV','MOSE-OM'))
  heatmapdata$cluster = factor(heatmapdata$cluster,levels = c('epi0','epi1','epi2','epi3'))
  ggplot(heatmapdata,aes(x = group,y=cluster,fill = value))+
    geom_tile()+scale_fill_gradient(high = "#FF8787",low = "#FFF7EC")+
    geom_tile()+
    geom_text(aes(label = round(value,2)))+theme_transparent()+
    theme(axis.text.x = element_text(angle = 45,vjust = 0.5,face = "bold"))+ylab("")+xlab("")+
    theme(axis.text.y = element_text(face = "bold"))
  ggsave(filename = "mose_epi_roie.pdf",width = 3.3,height = 3.9)
}

#functional annotation of epi subclusters
epi.marker2 = epi.marker[[2]]

##GO enrichment
{
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  
  DE = dplyr::filter(epi.marker2,avg_log2FC > 0.585 & p_val_adj <= 0.05)
  
  TC_DE = lapply(c('epi0','epi1','epi2','epi3'),function(celltype){
    dt = dplyr::filter(DE,cluster == celltype)
    return(dt)
  })
  
  ego1 = enrichGO(gene = rownames(TC_DE[[1]]),keyType = 'SYMBOL', #'ENTREZID','ENSEMBL'
                  OrgDb = org.Mm.eg.db,ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 0.01,qvalueCutoff = 0.05,readable = FALSE
  )
  ego2 = enrichGO(gene = rownames(TC_DE[[2]]),keyType = 'SYMBOL', #'ENTREZID','ENSEMBL'
                  OrgDb = org.Mm.eg.db,ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 0.01,qvalueCutoff = 0.05,readable = FALSE
  )
  ego3 = enrichGO(gene = rownames(TC_DE[[3]]),keyType = 'SYMBOL', #'ENTREZID','ENSEMBL'
                  OrgDb = org.Mm.eg.db,ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 0.01,qvalueCutoff = 0.05,readable = FALSE
  )
  ego4 = enrichGO(gene = rownames(TC_DE[[4]]),keyType = 'SYMBOL', #'ENTREZID','ENSEMBL'
                  OrgDb = org.Mm.eg.db,ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 0.01,qvalueCutoff = 0.05,readable = FALSE
  )
  
  dotplot(ego1)+ggtitle('epi0')
  p2 = dotplot(ego2) +ggtitle('Epi1');p2
  p3 = dotplot(ego3)+ggtitle('Epi2');p3
  
  
  extractego = function(ego,name){
    sub = ego@result
    sub$Cluster = name
    return(sub)
  }
  
  ego1_df = extractego(ego1,name = 'epi0')
  ego2_df = extractego(ego2,name = 'epi1')
  ego3_df = extractego(ego3,name = 'epi2')
  ego4_df = extractego(ego4,name = 'epi3')
  
  ego1_df$P1 = -log10(ego1_df$p.adjust)
  ego2_df$P1 = -log10(ego2_df$p.adjust)
  ego3_df$P1 = -log10(ego3_df$p.adjust)
  ego4_df$P1 = -log10(ego4_df$p.adjust)
  
  ego_result = dplyr::bind_rows(ego1_df,ego2_df,ego3_df,ego4_df)
  writexl::write_xlsx(ego_result,'ego_bp_epi_result.xlsx')
  write.csv(ego_result,'ego_bp_epi_result.csv')
  
  ##select pathway and make figure
  {
    
    epi.color
    p_epi0 = ggplot(subset(ego1_df,ID %in% c('GO:0007015','GO:0006631','GO:0007266','GO:0061572','GO:0019216')), aes(x = reorder(Description, -P1), y = -log10(p.adjust))) +
      geom_bar(stat = "identity", fill = "#CD8D7A") +
      coord_flip() +
      theme_pubr() +
      labs(title = 'epi0',
           y = "-log10(p.adjust)",
           x = "Description") +
      geom_text(aes(label = Description), hjust = 1, color = "black", size = 4, position = position_nudge(x = 0))+
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.title = element_text(hjust = 0))
    
    p_epi1 = ggplot(subset(ego2_df,ID %in% c('GO:0050678','GO:0033002','GO:0060562','GO:0010631','GO:0090132')), aes(x = reorder(Description, -P1), y = -log10(p.adjust))) +
      geom_bar(stat = "identity", fill = "#A0C49D") +
      coord_flip() +
      theme_pubr() +
      labs(title = 'epi1',
           y = "-log10(p.adjust)",
           x = "Description") +
      geom_text(aes(label = Description), hjust = 1, color = "black", size = 4, position = position_nudge(x = 0))+
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.title = element_text(hjust = 0))
    
    p_epi2 = ggplot(subset(ego3_df,ID %in% c('GO:0007059','GO:0098813','GO:0000280','GO:0006260','GO:0140014')), aes(x = reorder(Description, -P1), y = -log10(p.adjust))) +
      geom_bar(stat = "identity", fill = "#96B6C5") +
      coord_flip() +
      theme_pubr() +
      labs(title = 'epi2',
           y = "-log10(p.adjust)",
           x = "Description") +
      geom_text(aes(label = Description), hjust = 1, color = "black", size = 4, position = position_nudge(x = 0))+
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.title = element_text(hjust = 0))
    
    p_epi3 = ggplot(subset(ego4_df,ID %in% c('GO:0001666','GO:0036293','GO:0071456','GO:0036294','GO:0097193')), aes(x = reorder(Description, -P1), y = -log10(p.adjust))) +
      geom_bar(stat = "identity", fill = "#607274") +
      coord_flip() +
      theme_pubr() +
      labs(title = 'epi3',
           y = "-log10(p.adjust)",
           x = "Description") +
      geom_text(aes(label = Description), hjust = 1, color = "black", size = 4, position = position_nudge(x = 0))+
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.title = element_text(hjust = 0))
    
    cowplot::plot_grid(p_epi0,p_epi1,p_epi2,p_epi3,align = "h",ncol = 4)
    ggsave('epi_go_bp.pdf',width= 11.4,height = 2.5)
  }
  
  
  
}


##progeny
{
  library(progeny)
  CellsClusters <- data.frame(Cell = names(Idents(epi)), 
                              CellType = as.character(Idents(epi)),
                              stringsAsFactors = FALSE)
  
  progeny_obj = progeny(epi,scale = F,organism = 'Mouse',top = 500,return_assay = TRUE)
  
  progeny_obj <- Seurat::ScaleData(progeny_obj, assay = "progeny") 
  
  progeny_scores_df <- 
    as.data.frame(t(GetAssayData(progeny_obj, slot = "scale.data", 
                                 assay = "progeny"))) %>%
    rownames_to_column("Cell") %>%
    gather(Pathway, Activity, -Cell) 
  
  
  progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters,by = 'Cell')
  
  summarized_progeny_scores <- progeny_scores_df %>% 
    group_by(Pathway, CellType) %>%
    summarise(avg = mean(Activity), std = sd(Activity))
  
  summarized_progeny_scores_df <- summarized_progeny_scores %>%
    dplyr::select(-std) %>%   
    spread(Pathway, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 
  
  
  
  htmap = ComplexHeatmap::Heatmap(summarized_progeny_scores_df,cluster_rows = F,cluster_columns = F,
                                  col = colorRamp2(c(-1.5, 0, 1.5), c("#1679AB", "white", "#C80036")))
  
  pdf('progeny_epi.pdf',width = 4.6,height = 2.1)
  draw(htmap)
  dev.off()
}


#score of YHDU
{
  library(msigdbr)
  library(fgsea)
  library(GSVA)
  
  list.files('/home/cwt/Genelist/')
  HM = gmtPathways('/home/cwt/Genelist/mh.all.v2024.1.Mm.symbols.gmt')  
  
  dim(epi)
  dcelltype = unique(epi@meta.data$dcelltype)
  
  source('/home/cwt/Genelist/scores.r') ##serua5V5 change data acquistion way

  HMscore <- score_cells(seur=epi, names=HM, combine_genes='mean', 
                         groups=NULL, group_stat='mean', cells.use=NULL)
  HMmatrix = as.matrix(HMscore) %>% t
  HMmatrix = apply(HMmatrix, 2, function(x)signif(x,digits = 3))
  colnames(HMmatrix) = rownames(epi@meta.data)
  HMseuobj = CreateSeuratObject(counts = HMmatrix,data = HMmatrix,meta.data = epi@meta.data)
  
  epi@assays$HM = HMseuobj@assays$RNA
  DefaultAssay(epi) = 'HM'
  
  table(epi$group)
  Diff_HM = lapply(dcelltype,function(i){
    sub = subset(epi,subset = dcelltype == i)
    DEfeatures =FindMarkers(sub,ident.1 = 'MOSE-OM',ident.2 = 'MOSE-OV',group.by = 'group',
                            min.pct=0.1,logfc.threshold = 0,pseudocount.use = F,only.pos = F)
    DEfeatures$terms = rownames(DEfeatures)
    DEfeatures$celltype = i
    return(DEfeatures)
    print(i)
  })
  
  process_diff_combined <- function(Diff_HM, dcelltype, selected_hallmark) {
    # Set names for the list
    names(Diff_HM) <- dcelltype
    
    # Combine the list into a single dataframe
    Diff_combined <- dplyr::bind_rows(Diff_HM)
    rownames(Diff_combined) <- NULL
    
    # Modify the Pathway column
    Diff_combined <- Diff_combined %>%
      mutate(Pathway = gsub("HALLMARK-", "", terms))
    
    # Modify the Pathway names
    Diff_combined$Pathway <- gsub('-', '_', Diff_combined$Pathway)
    
    # Filter and factor the Diff_combined dataframe
    Diff_combined_selected <- dplyr::filter(Diff_combined, Pathway %in% selected_hallmark$Name)
    Diff_combined_selected$Pathway <- factor(Diff_combined_selected$Pathway, levels = selected_hallmark$Name)
    
    return(Diff_combined_selected)
  }
  #selected_hallmark = readxl::read_xlsx('/home/cwt/Genelist/HALLMARK_CATEGORY.xlsx')
  
  plot.data = process_diff_combined(Diff_HM,dcelltype,selected_hallmark)
  
  ggplot(plot.data,aes(y = celltype,x = Pathway))+
    geom_point(aes(color=avg_log2FC,size=-log10(p_val_adj+10^-100)))+
    scale_color_gradient2(low = "#769FCD", mid = "white", high = "#FC5185", midpoint = 0, 
                          limits = c(-0.5, 0.5))+
    labs(title = 'HALLMARK:OM vs OV')+
    scale_size_continuous(range=c(1,6),breaks=seq(0,30,10),name="FDR")+
    theme(panel.background=element_rect(color="black",fill=NA),legend.key =  element_blank(),axis.ticks=element_blank(),
          panel.grid=element_blank(),axis.text.x=element_text(size=8,angle=90,vjust = 0,hjust = 1,color="black"),
          axis.text.y=element_text(size=8,color="black"),axis.title=element_blank(),legend.position = "bottom")+
    guides(size = guide_legend(nrow = 2))
  
  ggsave('hallmark_epi_omvsov.pdf',width = 9.3,height= 4.7)
  
}


#GSEA analysis
library(GseaVis)

#DEG analysis of OM and OV of epithelial cells
##分别分析epi0和epi1原发灶和转移灶之间的差异基因
if(F){
  deg_epi0 = FindMarkers(subset(epi,idents = 'epi0'),ident.1 = 'MOSE-OM',ident.2 = 'MOSE-OV',group.by = 'group',
                         min.pct=0.1,logfc.threshold = 0,pseudocount.use = F,only.pos = F)
  deg_epi1 = FindMarkers(subset(epi,idents = 'epi1'),ident.1 = 'MOSE-OM',ident.2 = 'MOSE-OV',group.by = 'group',
                         min.pct=0.1,logfc.threshold = 0,pseudocount.use = F,only.pos = F)
  
  deg_epi0_plot = deg_epi0 %>% 
    mutate(significant = case_when(
      avg_log2FC >= 0.585 & p_val_adj < 0.05 ~ "up",
      avg_log2FC <= -0.585 & p_val_adj < 0.05 ~ "down",
      TRUE ~ "stable"
    ))
  deg_epi1_plot = deg_epi1 %>% 
    mutate(significant = case_when(
      avg_log2FC >= 0.585 & p_val_adj < 0.05 ~ "up",
      avg_log2FC <= -0.585 & p_val_adj < 0.05 ~ "down",
      TRUE ~ "stable"
    ))
  
  wt_volcano = function(data,avg_log2FC,p_val_adj){
    ggplot(data,aes(avg_log2FC,-1*log10(p_val_adj)))+
      geom_point(aes(color=significant),size=0.8)+theme_classic()+
      scale_color_manual(values = c("#4CBBD5","grey","#f39b7f"))+
      geom_hline(yintercept = 1.3,linetype=4,size=0.3)+
      geom_vline(xintercept = c(-0.5,0.5),linetype=4,size=0.3)+
      theme(title=element_text(size = 18),text = element_text(size=18),legend.position = 'none')+
      labs(x="log2(Fold Change)",y="-log10(Adjusted Pvalue)")
  }
  p1 = wt_volcano(deg_epi0_plot)+ggtitle('epi0')
  p2 = wt_volcano(deg_epi1_plot)+ggtitle('epi1')
  cowplot::plot_grid(plotlist = list(p1,p2),ncol = 2)
  ggsave('epi_volcano_omvsov.pdf',width = 5.2,height = 2.7)
}

#GSEA分析EPI0和EPI1各自的通路富集 --hallmark通路
##epi0 在原发灶中富集Glycolysis  epi1在转移灶中富集emt
{
  library(fgsea)
  library(cowplot)
  
  HM
  
  epi0_genelist = deg_epi0 %>% rownames_to_column(var = "gene") %>% 
    dplyr::filter(avg_log2FC <0) %>% arrange(desc(avg_log2FC)) %>% dplyr::select(gene,avg_log2FC) %>% deframe()
  
  epi0_fgsea = fgsea(HM,stats = epi0_genelist)
  p1 = plotEnrichment(HM[["HALLMARK_GLYCOLYSIS"]],
                      epi0_genelist) + labs(title="GLYCOLYSIS")
  
  epi1_genelist = deg_epi1 %>% rownames_to_column(var = "gene") %>% 
    dplyr::filter(avg_log2FC >0) %>% 
    mutate(avg_log2FC = ifelse(is.infinite(avg_log2FC), 6, avg_log2FC)) %>% 
    arrange(desc(avg_log2FC)) %>% dplyr::select(gene,avg_log2FC) %>% deframe()
  
  p2 = plotEnrichment(HM[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]],
                      epi1_genelist) + labs(title="EMT")
  
  plot_grid(plotlist = list(p1,p2),ncol =2)
  ggsave('epi_gsea_omvsov.pdf',width = 4.7,height = 2.3)
  
}

##cluster similarity comparison of MOSE-EPI
if(F){
  #install.packages('homologene')
  library(homologene)
  library(ggplot2)
  library(GeneOverlap)
  
  wt_cluster_comparison = function(genelist1,genelist2){
    require(GeneOverlap)
    require(RColorBrewer)
    gom.obj = newGOM(genelist1,genelist2)
    gom.obj
    print(gom.obj)
    pvaltable = getMatrix(gom.obj,name = "pval")
    intersectiontable = getMatrix(gom.obj,name = "intersection")
    union = getMatrix(gom.obj,name = "union")
    odds = getMatrix(gom.obj,name = 'odds.ratio')
    
    pvaltable = as.data.frame(pvaltable)
    pvaltable$genelist1 = rownames(pvaltable)
    p.spread = pvaltable %>% gather(genelist2,pvalue,-genelist1)
    
    intersectiontable = as.data.frame(intersectiontable)
    intersectiontable$genelist1 = rownames(intersectiontable)
    inter.spread = intersectiontable %>% gather(genelist2,overlap,-genelist1)
    
    union = as.data.frame(union)
    union$genelist1 = rownames(union)
    union.spread = union %>% gather(genelist2,union,-genelist1)
    
    odds = as.data.frame(odds)
    odds$genelist1 = rownames(odds)
    odds.spread = odds %>% gather(genelist2,odds,-genelist1)
    
    
    T.comparison = left_join(p.spread,inter.spread)
    T.comparison = left_join(T.comparison,union.spread)
    T.comparison = left_join(T.comparison,odds.spread)
    T.comparison$pct = T.comparison$overlap/(T.comparison$union+T.comparison$overlap)
    
    palette = brewer.pal(11,'RdBu')
    T.comparison$genelist2 = factor(T.comparison$genelist2,levels = paste0(names(genelist2)))
    T.comparison$genelist1 = factor(T.comparison$genelist1,levels = paste0(names(genelist1)))
    p = ggplot(T.comparison,aes(x = genelist1,y = genelist2))+
      geom_point(aes(color = odds,size = -log10(pvalue+10^-100)))+
      scale_color_gradientn(limit=c(0,10),colors= colorRampPalette(brewer.pal(8,'OrRd'))(10))+
      theme(panel.background=element_rect(color="black",fill=NA),legend.key =  element_blank(),axis.ticks=element_blank(),
            panel.grid=element_blank(),axis.text.x=element_text(size=8,angle=90,vjust = 0.5,hjust = 0.5,color="black"),
            axis.text.y=element_text(size=8,color="black"),axis.title=element_blank(),legend.position = "bottom")+
      geom_vline(aes(xintercept = genelist1),linetype = 'dashed',alpha = 0.5,color = 'gray')+
      geom_hline(aes(yintercept = genelist2),linetype = 'dashed',alpha = 0.5,color = 'gray')
    return(p)
  }
  filter_data <- function(data, clustername) {
    result_list <- list()
    
    for (col_name in clustername) {
      filtered_column = data %>% filter(cluster == col_name) %>% select(gene) %>% pull()
      
      result_list[[col_name]] <- filtered_column
    }
    
    return(result_list)
  }
  filter_data_m <- function(data, clustername) {
    result_list <- list()
    
    for (col_name in clustername) {
      filtered_column = data %>% filter(cluster == col_name) %>% select('9606') %>% pull()
      
      result_list[[col_name]] <- filtered_column
    }
    
    return(result_list)
  }
  
  DE_RNA = epi.marker2
  genelist = DE_RNA$gene
  genetrans = homologene(genelist,inTax = 10090,outTax = 9606)
  DE_RNA_merge = merge(DE_RNA,genetrans,by.x = 'gene',by.y='10090')
  
  
  #spectrum_comparison
  {
    getwd()
    spectrum_tumor = readxl::read_xlsx('/home/cwt/Genelist/2022_Nature_cancer_marker_.xlsx')
    colnames(spectrum_tumor)[2] = 'cluster'
    
    clustername_spectrum = unique(spectrum_tumor$cluster)
    spectrum_tumor_list = filter_data(spectrum_tumor,clustername  = clustername_spectrum)
    
    MEPP_tumor_list = filter_data_m(DE_RNA_merge, clustername = unique(DE_RNA$cluster))
    
    p = wt_cluster_comparison(spectrum_tumor_list,MEPP_tumor_list)
    p
    pdf('/work/cwt/ID8/2_comparison/spectrum_comparison.pdf',width = 8,height = 5)
    p
    dev.off()
  }
  
  
  #chai comparison
  chai = read.csv('/home/cwt/Genelist/EPI_each_cluster_allDEGs_toChen.csv')
  chai = readxl::read_xlsx('/home/cwt/Genelist/2024_Chai_epi.xlsx')
  {
    chai_gene = chai$gene
    genetrans = homologene(chai_gene,inTax = 9606,outTax = 10090)
    chai_merge = merge(chai,genetrans,by.x = 'gene',by.y='9606')
    colnames(chai_merge)[c(1,9)] = c('hgene','mgene')
    
    unique(chai_merge$hgene)
    
    
    cluster.chai = c('0','1','2','3','4','5','6','7','8','9','10','11')
    
    chai.list = list()
    for (col_name in cluster.chai){
      column = chai_merge %>% dplyr::filter(cluster == col_name) %>% select('hgene') %>% pull()
      chai.list[[col_name]] = column
    }
    chai.list
    
    names(chai.list) = c('Epi0','Epi1','Epi2','Epi3','Epi4','Epi5','Epi6','Epi7','Epi8','Epi9','Epi10','Epi11')
    
    p.chai = wt_cluster_comparison(chai.list,MEPP_tumor_list)
    p.chai
    
    cowplot::plot_grid(plotlist = list(p,p.chai))
    ggsave('epi_cluster_comparison.pdf',width = 11,height = 4)
  }
}


DefaultAssay(epi) = 'RNA'
p0 = FeaturePlot(epi,features = 'Runx2')
p1 = DotPlot(epi,features = c('Runx2'))
p2 = VlnPlot(epi,features = c('Runx2'),group.by = 'group')+
  stat_compare_means()
cowplot::plot_grid(plotlist = list(p0,p1,p2),ncol = 3)

#####TCGA survival analysis
{
  library(survminer)
  library(survival)
  library(limma)
  list.files('/home/cwt/Project/TCGA_OV/')
  
  tcga_ov_expr = read_rds("/home/cwt/Project/TCGA_OV/OV_TCGA_Expression.rds.gz")
  tcga_ov_survival = read_rds("/home/cwt/Project/TCGA_OV/tcga_survival.rds.gz")
  expr1 = normalizeBetweenArrays(tcga_ov_expr)
  
  expr1 = as.data.frame(t(expr1))
  survival_expr = merge.data.frame(expr1,tcga_ov_survival,by.x = 0,by.y = 'sample')
  
  selected.gene = c('PTEN','TP53')
  dt1 = surv_cutpoint(survival_expr,time = 'OS.time',event = 'OS',variables = selected.gene)
  dt2 = surv_categorize(dt1)
  
  {
    formulae <- list()
    genes <- colnames(dt2)[3:ncol(dt2)]
    formulae <- lapply(genes, function(x) as.formula(paste0("Surv(OS.time,OS) ~ ", x)))
    fits <- surv_fit(formulae, data = dt2)
    
    p <- ggsurvplot_list(fits,
                         data = dt2,
                         risk.table = FALSE,
                         pval = TRUE,
                         palette = c('#C96868','#B5C0D0'),
                         break.time.by = 1000,
                         #ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE)
    
    getwd()
    pdf('surv_pten_tp53.pdf',width = 6.2,height = 3.52)
    p_list = arrange_ggsurvplots(p[1:length(p)], print = TRUE, ncol = 2,nrow = 1)
    dev.off()
  }
  
  dt2 = as_tibble(dt2)
  dt2 = dt2 %>% mutate(test = case_when(
    PTEN == 'low' & TP53 == 'low' ~ 'PTEN-LOW/TP53-LOW',
    PTEN == 'low' & TP53 == 'high' ~ 'PTEN-LOW/TP53-HIGH',
    PTEN == 'high' & TP53 == 'low' ~ 'PTEN-HIGH/TP53-LOW',
    PTEN == 'high' & TP53 == 'high' ~ 'PTEN-HIGH/TP53-HIGH'
  ))
  ggsurvplot(survfit(Surv(OS.time,OS)~ test,data = dt2),pval = T,
             risk.table= F,palette = c('#C96868','#B5C0D0','#FADFA1','#EED3D9'))
  
  
}


#qsave(epi,file = '1_clustering/epi_1112.qs') #1112

#Immune clusters

getwd()
dim(wholedata)
DefaultAssay(wholedata) = 'RNA'

table(wholedata$bigclass)

``````````````````Macrophage``````````````````````````````````````````````````````
#data processing
if(F){
  mac = subset(wholedata,idents = 'Macrophage')
  mac = wt_firststep(mac)
  
  DimPlot(mac)
  mac@active.ident = mac$RNA_snn_res.0.2
  mac.marker = wt_downsamplefindmarker_mouse(mac)
  
  mac1 = subset(mac,idents = c('6','7'),invert = T)
  mac1 = wt_firststep(mac1)
  mac1@active.ident = mac1$RNA_snn_res.0.2
  mac.marker = wt_downsamplefindmarker_mouse(mac1)
  
  DimPlot(mac1)
  FeaturePlot(mac1,features = c('Cxcl9','Spp1'))
  
  DotPlot(mac1,features = c('Wt1','Dcn','Cd3e'))
  
  mac1 = subset(mac1,idents = '6',invert= T)
  DimPlot(mac1)
  
  mac.marker = wt_downsamplefindmarker_mouse(mac1)
  writexl::write_xlsx(mac.marker,'macmarker.xlsx')
  
  
  mac1 = RenameIdents(mac1, `0` = 'mac1_Cx3cr1',`1` = 'mac2_Lyve1',`2`= 'mac3_Il7r',`3` = 'mac4_Cxcl9',`4` = 'mac5_Il1b',`5` = 'mac6_Serpinb2',`7` ='mac2_Lyve1')
  
  mac1@active.ident = factor(mac1@active.ident, levels = c('mac1_Cx3cr1','mac2_Lyve1','mac3_Il7r','mac4_Cxcl9','mac5_Il1b','mac6_Serpinb2'))
  mac1$dcelltype = mac1@active.ident
  
  mac.color = c('#B7E0FF','#E78F81','#77E4C8','#478CCF','#1F4E5F','#FF4545')
  names(mac.color) = unique(mac1$dcelltype)
  qsave(mac1,'mac1118.qs')
  
  
}

#figures
if(F){
  p1 = DimPlot(mac1,cols = mac.color);p1
  ggsave('mac.dimplot.pdf',path = '1_clustering/',width = 4.6,height = 3.5)
  
  mac.marker = wt_downsamplefindmarker_mouse(mac1)
  writexl::write_xlsx(mac.marker,'macmarker.xlsx')
  
  
  p2 = DimPlot(mac1,label= F,cols = mac.color,split.by = 'group',pt.size = 0.001)
  ggsave('mac.split.dim.pdf',width = 7.61,height = 3.6,path = '1_clustering/')
  
  p3 = DotPlot(mac,features = c(
  'Cx3cr1','Cxcl14',
  'Lyve1','Cd163','Ccl8',
  'Il7r','Procr',
  'Cxcl9','Cxcl10',
  'Il1b','Plac8',
  'Serpinb2','Ptgis'
    ))+theme(axis.text.x = element_text(angle = 90,hjust=1,vjust = 0.5),plot.title = element_text(hjust = 0.5,size = 8, face = "bold"))+scale_color_distiller(palette = "RdBu")
  ggsave('mac.dotplot.pdf',path = '1_clustering/',width = 7,height = 2.5)
  
  
}

mac = mac1
#clusterprofiling
if(F){
  mac.meta = mac@meta.data
  mac.meta$dcelltype = factor(mac.meta$dcelltype,levels = c(
    'mac1_Cx3cr1','mac2_Lyve1','mac3_Il7r','mac4_Cxcl9','mac5_Il1b','mac6_Serpinb2'
  ))
  ClusterFreq <- mac.meta[,c("dcelltype","group")] %>% table %>%
    data.frame() %>% set_colnames(c("seurat_clusters","stim","Per"))
  ClusterPer <- ClusterFreq %>% tidyr::spread(seurat_clusters,Per)
  ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
  ClusterPer1 <- ClusterPer %>% tidyr::gather(key=seurat_clusters,value=Per,-stim)
  
  ClusterPer1$seurat_clusters = factor(ClusterPer1$seurat_clusters,levels = c(
    'mac1_Cx3cr1','mac2_Lyve1','mac3_Il7r','mac4_Cxcl9','mac5_Il1b','mac6_Serpinb2'
  ))
  
  p4 = ggplot(ClusterPer1,  aes(x = seurat_clusters, y = Per, color = stim)) +
    geom_point(position = position_dodge(width = 0.5)) +
    theme_pubr() +
    labs(y = 'Proportion among all cells (%)', x = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = 'top') +
    scale_color_manual(values = c("MOSE-OV" = "#92CCE1", "MOSE-OM" = "#F68686"))
  ggsave('mac.dot.proportion.pdf',width = 3.7,height = 3.5,path = '1_clustering/')
}
p_mac = cowplot::plot_grid(plotlist = list(p2,p3,p4),ncol=3,rel_widths = c(2,2,1));p_mac

#m1 and m2 of macrophages?

``````````````````````````````````````````````````TNK``````````````````````````````````````````````````
#data processing
if(F){
  table(wholedata$bigclass)
  
  tnk = subset(wholedata,idents = 'TNK')
  
  tnk = wt_firststep(tnk)
  DimPlot(tnk,group.by = 'RNA_snn_res.0.1')
  tnk@active.ident = tnk$RNA_snn_res.0.1
  tnk.marker = wt_downsamplefindmarker_mouse(tnk)
  
  tnk1 = subset(tnk1,idents = c('8'),invert= T)
  
  tnk1 = wt_firststep(tnk1)
  DimPlot(tnk1,group.by = 'RNA_snn_res.0.1')
  tnk1@active.ident = tnk1$RNA_snn_res.0.1
  tnk.marker = wt_downsamplefindmarker_mouse(tnk1)
  writexl::write_xlsx(tnk.marker,'tnkmarker.xlsx')
  
  
  DotPlot(tnk1,features = c('Cd3e','Cd4','Cd8a'))
  DotPlot(tnk1,features = c('Cd4','Cd8a','Fcer1g','Klrc1','Foxp3','Il7r','Lef1','Tcf7','Cxcl13','Pdcd1','Fkbp5','Itga6'))
  
  tnk = RenameIdents(tnk1,`0` = 't1_Treg',`1` = 't2_CD4_naive',`2` = 't5_NK',`5` = 't3_CD8T',`6` = 't4_Cxcl10')
  tnk@active.ident = factor(tnk@active.ident,levels = c('t1_Treg','t2_CD4_naive','t3_CD8T','t4_Cxcl10','t5_NK'))
  tnk$dcelltype = tnk@active.ident
  
  qsave(tnk,file = 'tnk1118.qs')
}

#Figures
if(F){
  tnk.color = c('#8D448B','#01AAC1','#F9A828','#87A922','#EA728C')
  names(tnk.color) = unique(tnk$dcelltype)
  
  p1 = DimPlot(tnk,cols = tnk.color,pt.size = 0.001);p1
  ggsave('tnk.dimplot.pdf',path = '1_clustering/',width = 4.6,height = 3.5)
  
  tnk.marker = wt_downsamplefindmarker_mouse(tnk)
  writexl::write_xlsx(tnk.marker,'tnkmarker.xlsx')
  
  
  p2 = DimPlot(tnk,label= F,cols = tnk.color,split.by = 'group',pt.size = 0.001)
  ggsave('tnk.split.dim.pdf',width = 7.61,height = 3.6,path = '1_clustering/')
  
  p3 = DotPlot(tnk,features = c(
'Cd3g','Ctla4','Foxp3','Cd4','Il7r','Cd8a','Ccr7','Fkbp2','Cxcl10','Slfn5','Klra4','Gzma'
  ))+theme(axis.text.x = element_text(angle = 90,hjust=1,vjust = 0.5),plot.title = element_text(hjust = 0.5,size = 8, face = "bold"))+scale_color_distiller(palette = "RdBu")
  ggsave('tnk.dotplot.pdf',path = '1_clustering/',width = 7,height = 3)
}

#clusterprofiling
if(F){
  tnk.meta = tnk@meta.data
  tnk.meta$dcelltype = factor(tnk.meta$dcelltype,levels = c(
    't1_Treg','t2_CD4_naive','t3_CD8T','t4_Cxcl10','t5_NK'
  ))
  ClusterFreq <- tnk.meta[,c("dcelltype","group")] %>% table %>%
    data.frame() %>% set_colnames(c("seurat_clusters","stim","Per"))
  ClusterPer <- ClusterFreq %>% tidyr::spread(seurat_clusters,Per)
  ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
  ClusterPer1 <- ClusterPer %>% tidyr::gather(key=seurat_clusters,value=Per,-stim)
  
  ClusterPer1$seurat_clusters = factor(ClusterPer1$seurat_clusters,levels = c(
    't1_Treg','t2_CD4_naive','t3_CD8T','t4_Cxcl10','t5_NK'
  ))
  
  p4 = ggplot(ClusterPer1,  aes(x = seurat_clusters, y = Per, color = stim)) +
    geom_point(position = position_dodge(width = 0.5)) +
    theme_pubr() +
    labs(y = 'Proportion among all cells (%)', x = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = 'top') +
    scale_color_manual(values = c("MOSE-OV" = "#92CCE1", "MOSE-OM" = "#F68686"))
  ggsave('tnk.dot.proportion.pdf',width = 3.7,height = 3.5,path = '1_clustering/')
}

p_tnk = cowplot::plot_grid(plotlist = list(p2,p3,p4),ncol=3,rel_widths = c(2,2,1))
p_mac/p_tnk
ggsave(filename = 'plot_tnk_mac.pdf',width = 16.5,height = 6.45)
qsave(wholedata,'wholedata1118.qs')

``````````````````````````````````````````````````Neutrophil``````````````````````````````````````````````````
neu = subset(wholedata,idents = c('Neutrophil'))
neu = wt_firststep(neu)

DimPlot(neu,group.by = 'RNA_snn_res.0.1')
neu@active.ident = neu$RNA_snn_res.0.1
neu.marker = wt_downsamplefindmarker_mouse(neu)
writexl::write_xlsx(neu.marker,'neumarker.xlsx')

neu_ori = neu
neu = RenameIdents(neu_ori,`0` = 'n1_Rsad2',`1` = 'n2_Ngp',`2` = 'n3_Atp5h',`3` = 'n4_Alpk1')


neu@active.ident = factor(neu@active.ident,levels = c('n1_Rsad2','n2_Ngp','n3_Atp5h','n4_Alpk1'))
neu$dcelltype = neu@active.ident

if(F){
  neu.color = c('#605EA1','#8EA3A6','#F6D6D6','#7BD3EA')
  names(neu.color) = unique(neu$dcelltype)
  
  p1 = DimPlot(neu,cols = neu.color,pt.size = 0.001);p1
  ggsave('neu.dimplot.pdf',path = '1_clustering/',width = 4.6,height = 3.5)
  
  neu.marker = wt_downsamplefindmarker_mouse(neu)
  writexl::write_xlsx(neu.marker,'neumarker.xlsx')
  
  
  p2 = DimPlot(neu,label= F,cols = neu.color,split.by = 'group',pt.size = 0.001)
  ggsave('neu.split.dim.pdf',width = 7.61,height = 3.6,path = '1_clustering/')

}

if(F){
  neu.meta = neu@meta.data
  neu.meta$dcelltype = factor(neu.meta$dcelltype,levels = c(
    'n1_Rsad2','n2_Ngp','n3_Atp5h','n4_Alpk1'
  ))
  ClusterFreq <- neu.meta[,c("dcelltype","group")] %>% table %>%
    data.frame() %>% set_colnames(c("seurat_clusters","stim","Per"))
  ClusterPer <- ClusterFreq %>% tidyr::spread(seurat_clusters,Per)
  ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
  ClusterPer1 <- ClusterPer %>% tidyr::gather(key=seurat_clusters,value=Per,-stim)
  
  ClusterPer1$seurat_clusters = factor(ClusterPer1$seurat_clusters,levels = c(
    'n1_Rsad2','n2_Ngp','n3_Atp5h','n4_Alpk1'
  ))
  
  p4 = ggplot(ClusterPer1,  aes(x = seurat_clusters, y = Per, color = stim)) +
    geom_point(position = position_dodge(width = 0.5)) +
    theme_pubr() +
    labs(y = 'Proportion among all cells (%)', x = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = 'top') +
    scale_color_manual(values = c("MOSE-OV" = "#92CCE1", "MOSE-OM" = "#F68686"));p4
  ggsave('neu.dot.proportion.pdf',width = 3.7,height = 3.5,path = '1_clustering/')
}

``````````````````````````````````````````````````Bcell``````````````````````````````````````````````
wholedata
table(wholedata$bigclass)
Bcell = subset(wholedata,idents = 'B')
Bcell = wt_firststep(Bcell)

DimPlot(Bcell)
DimPlot(Bcell,group.by = 'RNA_snn_res.0.3')
Bcell@active.ident = Bcell$RNA_snn_res.0.1
B.marker = FindAllMarkers(Bcell)

Bcell = subset(Bcell,idents = c('0'))
Bcell = wt_firststep(Bcell)

FeaturePlot(Bcell,features = c('Cd79a','Cd19','Cd20','Ighm','Ighg','Cd38','Cd14'),order = T)

FeaturePlot(Bcell,features = c('Cd138','Cd27','Cd38'))
FeaturePlot(Bcell,features = c('Il6','Il12','Tgfb'))
Bcell
p1 = DimPlot(Bcell,raster = F)
p2 = FeaturePlot(Bcell,features = c('Cd79a','Cd19','Ighm','Cd38','Cd14'),order = T,raster = F,ncol = 5);p2
cowplot::plot_grid(plotlist = list(p1,p2),rel_widths = c(1,5))
ggsave('Bplot.pdf',width = 16,height = 2)

qsave(Bcell,file = 'Bcell1126.qs')

````````````````````````cDC1````````````````````````````````
cdc = subset(wholedata,idents = 'cDC1')
cdc = wt_firststep(cdc)
DimPlot(cdc,group.by = 'RNA_snn_res.0.5')
cdc.marker = wt_downsamplefindmarker_mouse(cdc)
DotPlot(cdc,features = c('Clec9a','Xcr1','C1qc','S100a8','Cd3g'))

cdc@active.ident = cdc$RNA_snn_res.0.5
cdc = subset(cdc,idents = c('0','4'))
DimPlot(cdc)

p1 = DimPlot(cdc)
p2 = FeaturePlot(cdc,features = c('Clec9a','Xcr1','Irf8'),ncol = 3,order = T)
cowplot::plot_grid(plotlist = list(p1,p2),rel_widths = c(1,3))
ggsave('cdc_plot.pdf',width = 10,height = 2.3)
