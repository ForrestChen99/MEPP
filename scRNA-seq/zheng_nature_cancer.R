```zheng nature cancer`
hov = read_rds('/home/cwt/Project/Zheng_NC_2023/raw_object.rds')
hov
hov = UpdateSeuratObject(hov)
table(hov$Groups)

hov@version

meta.hov = hov@meta.data
counts.hov = hov@assays$RNA@counts

hov.seu = CreateSeuratObject(counts = counts.hov,meta.data = meta.hov)
rm(counts.hov)


hov.seu
hov.seu = hov.seu %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 3000)
VGENES = VariableFeatures(hov.seu)
VGENES=setdiff(VGENES,VGENES[grep("^MT|^RPL|^RPS",VGENES)])


hov.seu = ScaleData(hov.seu,features = VGENES,vars.to.regress = c("percent.mt","nCount_RNA","nFeature_RNA","percent.ribo")) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>% FindNeighbors(dims = 1:30)

hov.seu = hov.seu %>% RunUMAP(dims = 1:30) %>% FindClusters(resolution=0.1,algorithm =1)
table(hov.seu$Groups)

table(hov.seu$maintypes_2)

qsave(hov.seu,file= '/home/cwt/Project/Zheng_NC_2023/zheng_processed.qs')

##网膜和原发灶的差异
##提取所有tumor cell
hov.subset.tumor = subset(hov.seu,subset = maintypes_2 == 'Epithelial cells')
hov.tumor = subset(hov.subset.tumor,subset = Groups %in% c('Primary Tumor','Metastatic Tumor'))
hov.tumor
hov.tumor = hov.tumor %>% NormalizeData() %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData(vars.to.regress = c("percent.mt","nCount_RNA","nFeature_RNA","percent.ribo")) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>% FindNeighbors(dims = 1:30) %>% RunUMAP(dims = 1:30) %>% FindClusters(resolution=0.1,algorithm =1)
DimPlot(hov.tumor,raster =T,group.by = 'Groups')


gdata::keep(hov.tumor,epi,sure = T)


#harmony
table(hov.tumor1$orig.ident)
hov.tumor1 = subset(hov.tumor,subset = orig.ident %in% c('OCCC1','UOC1'),invert = T)
hov.tumor1$orig.ident = droplevels(hov.tumor1$orig.ident)

hov.harmony = RunHarmony(hov.tumor1,group.by.vars = 'orig.ident',assay.use = 'RNA',reduction = 'pca',project.dim = F)
harmony_embeddings = Embeddings(hov.harmony,'harmony')
hov.harmony = hov.harmony %>% RunUMAP(reduction = 'harmony',dims = 1:30) %>%
  FindNeighbors(reduction = 'harmony',dims = 1:30)
for(i in seq(0.1,0.3,by = 0.1)){
  hov.harmony = FindClusters(hov.harmony,resolution = i,algorithm = 1)
}

DimPlot(hov.harmony,raster =T,group.by = 'orig.ident')
DimPlot(hov.harmony,raster =T,group.by = 'RNA_snn_res.0.1')


hov.harmony@active.ident = hov.harmony$RNA_snn_res.0.1
hov.harmony = subset(hov.harmony,idents = c('3','5','6'),invert = T)
hov.marker = wt_downsamplefindmarker_human(hov.harmony)

DimPlot(hov.harmony,label = T)
DotPlot(hov.harmony,features= c('MKI67'))

hov.harmony = RenameIdents(hov.harmony,`0` = 'zheng_epi0',`1` = 'zheng_epi1',`2` = 'zheng_epi2',`4` = 'zheng_epi3')
hov.harmony$dcelltype = hov.harmony@active.ident

#figures for zheng et al.
if(F){
  p1 = DimPlot(hov.harmony,raster =T,cols = zheng.color,label = T)+NoLegend()
  p2 = DimPlot(hov.harmony,raster =T,cols = zheng.color,split.by = 'Groups')+NoLegend()
  p3 #ROIE
  
  cowplot::plot_grid(plotlist = list(p1,p2,p3),ncol = 3, rel_widths = c(1,2,1))
  getwd()
  ggsave(filename = 'zheng.dimplot.pdf',width = 12.5,height = 3.6)
}
##cluster similarity between
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
  epi.marker = wt_downsamplefindmarker_mouse(epi)
  DE_RNA = epi.marker[[2]]
  genelist = DE_RNA$gene
  genetrans = homologene(genelist,inTax = 10090,outTax = 9606)
  DE_RNA_merge = merge(DE_RNA,genetrans,by.x = 'gene',by.y='10090')
  
  DE_RNA_merge = DE_RNA_merge %>% dplyr::filter(avg_log2FC>0.585 & p_val_adj <0.05)
  #zheng_nature_cancer_comprison
  hov.marker = wt_downsamplefindmarker_human(hov.harmony)
  zheng = hov.marker[[2]]
  zheng = zheng %>% dplyr::filter(avg_log2FC>0.585 & p_val_adj <0.05)
  cluster.zheng = unique(zheng$cluster)
  zheng_list = filter_data(zheng,clustername = cluster.zheng)
  MEPP_tumor_list = filter_data_m(DE_RNA_merge, clustername = unique(DE_RNA$cluster))
  p = wt_cluster_comparison(zheng_list,MEPP_tumor_list);p
  p
  ggsave(filename = 'zheng_comparison.pdf',width = 5,height = 3.4)
}


#proportion change of epi clusters
if(F){
  table(hov.harmony$Groups)
  hov.harmony$Groups = droplevels(hov.harmony$Groups)
  
  zheng.meta = hov.harmony@meta.data
  anno = zheng.meta %>% dplyr::select(orig.ident,Groups)
  anno <- anno %>% distinct(orig.ident, .keep_all = TRUE)
  rownames(anno) = NULL
  
  colnames(zheng.meta)
  zheng.meta$dcelltype = factor(zheng.meta$dcelltype,levels = c(
    'zheng_epi0','zheng_epi1','zheng_epi2','zheng_epi3'
  ))
  ClusterFreq <- zheng.meta[,c("dcelltype","orig.ident")] %>% table %>%
    data.frame() %>% set_colnames(c("seurat_clusters","stim","Per"))
  ClusterPer <- ClusterFreq %>% tidyr::spread(seurat_clusters,Per)
  ClusterPer[,2:ncol(ClusterPer)] <- t(apply(ClusterPer[,2:ncol(ClusterPer)],1,function(x){x/sum(x)*100}))
  ClusterPer1 <- ClusterPer %>% tidyr::gather(key=seurat_clusters,value=Per,-stim)
  
  ClusterPer1$seurat_clusters = factor(ClusterPer1$seurat_clusters,levels = c(
    'zheng_epi0','zheng_epi1','zheng_epi2','zheng_epi3'
  ))

  ClusterPer2 = merge.data.frame(x = ClusterPer1, y = anno, by.x = 'stim',by.y = 'orig.ident')
  
  ClusterPer2 <- ClusterPer2 %>%
    mutate(Metastatic = if_else(Groups == "Metastatic Tumor ", "Primary", "Metastatic"))
  
  ClusterPer2$Metastatic = factor(ClusterPer2$Metastatic,levels = c('Primary','Metastatic'))
  
  
  ggplot(ClusterPer2,  aes(x = seurat_clusters, y = Per, color = Groups)) +
    geom_point(position = position_dodge(width = 0.5)) +
    theme_pubr() +
    labs(y = 'Proportion among all cells (%)', x = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = 'top') #+
    #scale_color_manual(values = c("MOSE-OV" = "#92CCE1", "MOSE-OM" = "#F68686"))
  
  ggsave('epi.dot.proportion.pdf',width = 3.1,height = 3.25,path = '1_clustering/')
}


#ROIE of epi clusters
{
  zheng.meta = hov.harmony@meta.data
  summary = table(zheng.meta[,c('dcelltype','Groups')])
  roe = as.data.frame(ROIE(summary))
  roe$cluster = rownames(roe)
  rownames(roe) = NULL
  res = data.frame()
  res = rbind(res,roe)
  summary(roe)
  
  heatmapdata = gather(roe,group,value,-c(cluster))
  str(heatmapdata)

  heatmapdata$group = factor(heatmapdata$group,levels = c('Primary Tumor','Metastatic Tumor'))
  heatmapdata$cluster = factor(heatmapdata$cluster,levels = c('zheng_epi0','zheng_epi1','zheng_epi2','zheng_epi3'))
  p3 = ggplot(heatmapdata,aes(x = group,y=cluster,fill = value))+
    geom_tile()+scale_fill_gradient(high = "#FF8787",low = "#FFF7EC")+
    geom_tile()+
    geom_text(aes(label = round(value,2)))+theme_transparent()+
    theme(axis.text.x = element_text(angle = 45,vjust = 0.5,face = "bold"))+ylab("")+xlab("")+
    theme(axis.text.y = element_text(face = "bold"))
  ggsave(filename = "zheng_epi_roie.pdf",width = 3.3,height = 3.9)
}


#原位灶和转移灶之间的差距 #gsea
list.files('/home/cwt/Genelist')
HM_human = gmtPathways('/home/cwt/Genelist/h.all.v2023.1.Hs.symbols.gmt')


zheng_deg = FindMarkers(hov.harmony,ident.1 = 'Metastatic Tumor',ident.2 = 'Primary Tumor',group.by = 'Groups',
                        min.pct=0.1,logfc.threshold = 0,pseudocount.use = F,only.pos = F)
zheng_up = zheng_deg %>% rownames_to_column(var = "gene") %>% 
  mutate(avg_log2FC = ifelse(is.infinite(avg_log2FC), 6, avg_log2FC)) %>% 
  arrange(desc(avg_log2FC)) %>% dplyr::select(gene,avg_log2FC) %>% deframe()

plotEnrichment(HM_human[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]],
               zheng_up) + labs(title="EMT")
  
zheng.up.gsea = fgsea(HM_human,stats = zheng_up)

##MOSE
HM_mouse = gmtPathways('/home/cwt/Genelist/mh.all.v2024.1.Mm.symbols.gmt')  


mose_deg =  FindMarkers(epi,ident.1 = 'MOSE-OM',ident.2 = 'MOSE-OV',group.by = 'group',
                        min.pct=0.1,logfc.threshold = 0,pseudocount.use = F,only.pos = F)
mose_up = mose_deg %>% rownames_to_column(var = "gene") %>% 
  mutate(avg_log2FC = ifelse(is.infinite(avg_log2FC), 6, avg_log2FC)) %>% 
  arrange(desc(avg_log2FC)) %>% dplyr::select(gene,avg_log2FC) %>% deframe()

plotEnrichment(HM_mouse[["HALLMARK_G2M_CHECKPOINT"]],
               mose_up)# + labs(title="EMT")

mose.up.gsea = fgsea(HM_mouse,stats = mose_up)

####Combine and plot
{
  zheng.dt = zheng.up.gsea %>% data.frame() %>% dplyr::mutate(nes_human = NES)
  mose.dt = mose.up.gsea  %>% data.frame() %>% dplyr::mutate(nes_mose = NES)
  combined_data <- merge(zheng.dt, mose.dt, by = "pathway")
  
  
  combined_data <- combined_data %>%
    mutate(
      pathway = str_remove(pathway, "HALLMARK_"),  # Remove "HALLMARK_"
      avg_nes = (nes_human + nes_mose) / 2,       # Calculate the average NES
      label = case_when(
        rank(-avg_nes) <= 5 ~ pathway,  # Top 5 by average NES
        rank(avg_nes) <= 5 ~ pathway,   # Bottom 5 by average NES
        TRUE ~ ""  # No label for other pathways
      )
    )
  
  combined_data = combined_data %>% mutate(
    color = case_when(
      nes_human > 0 & nes_mose > 0 ~ "#BA135D",
      nes_human < 0 & nes_mose < 0 ~ "#80C4E9",
      TRUE ~ "grey"
    ),
    size = -log10((padj.x+padj.y)/2 + 1e-100) # Scale the size based on -log10(pvalue)
  )
  fwrite(combined_data,'gsea_analysis_zheng.csv')
  
  require(ggrepel)
  ggplot(combined_data, aes(x = nes_mose, y = nes_human, color = color, size = size)) +
    geom_point(alpha = 0.7) +  # Adjust transparency for better visibility
    scale_color_identity() +   # Use the defined colors directly
    geom_text_repel(
      aes(label = label),
      size = 3,                # Text size
      max.overlaps = Inf,      # Avoid overlap
      box.padding = 0.5,       # Extra padding around labels
      point.padding = 0.5,color = 'Black'      # Space between points and labels
    ) +
    scale_size_continuous(range = c(1,4),name = 'Mean Padj') +  # Adjust the range of point sizes
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +  # Add diagonal reference line
    theme_minimal() +
    labs(
      title = "Zheng_nature_cancer",
      x = "Mouse metastatic versus primary (NES)",
      y = "Human metastatic versus primary (NES)",
      size = "-log10(mean(padj))",
      color = "Point Type"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    )
  ggsave(filename = 'gsea_nes_zheng.pdf',width = 4.21,height = 3.5)
  
  
}



#epi1's feature in Zheng
FeaturePlot(epi,features= c('Wt1','Pax8'))
epi1.marker = epi.marker[[2]] %>% dplyr::filter(cluster == 'epi1') %>%
  top_n(30,wt = avg_log2FC) %>% rownames()

genetrans = homologene(epi1.marker,inTax = 10090,outTax = 9606)
epi1.marker.trans = genetrans$`9606`
epi1.list = list(epi1.marker.trans)
names(epi1.list) = 'epi1.feature'

hov.harmony1 = AddModuleScore(hov.harmony,features = epi1.list,name= 'epi1')
VlnPlot(hov.harmony,features = c('epi11'))

#epi0
epi0.marker = epi.marker[[2]] %>% dplyr::filter(cluster == 'epi0') %>%
  top_n(30,wt = avg_log2FC) %>% rownames()

genetrans = homologene(epi0.marker,inTax = 10090,outTax = 9606)
epi0.marker.trans = genetrans$`9606`
epi0.list = list(epi0.marker.trans)

hov.harmony1 = AddModuleScore(hov.harmony1,features = epi0.list,name= 'epi0')
table(hov.harmony1$epi01)
DotPlot(hov.harmony1,features = c('epi01'))
VlnPlot(hov.harmony1,features = c('epi01'))

#epi2 and #epi3
epi2.marker = epi.marker[[2]] %>% dplyr::filter(cluster == 'epi2') %>%
  top_n(30,wt = avg_log2FC) %>% rownames()

genetrans = homologene(epi2.marker,inTax = 10090,outTax = 9606)
epi2.marker.trans = genetrans$`9606`
epi2.list = list(epi2.marker.trans)
hov.harmony1 = AddModuleScore(hov.harmony1,features = epi2.list,name= 'epi2')

epi3.marker = epi.marker[[2]] %>% dplyr::filter(cluster == 'epi3') %>%
  top_n(30,wt = avg_log2FC) %>% rownames()

genetrans = homologene(epi3.marker,inTax = 10090,outTax = 9606)
epi3.marker.trans = genetrans$`9606`
epi3.list = list(epi3.marker.trans)
hov.harmony1 = AddModuleScore(hov.harmony1,features = epi3.list,name= 'epi3')

zheng.color = c('#65799B','#E23E57','#83CC61','#2FC5CC')
VlnPlot(hov.harmony1,features= c('epi01','epi11','epi21','epi31'),pt.size = 0, col = zheng.color, ncol = 2)
ggsave(filename = 'module_vln.pdf',width= 3.3,height= 4.3)


````````````````````````````Heatmap to plot immune infiltration````````````````````````````````````````````````````````````


hov.seu = qread(file= '/home/cwt/Project/Zheng_NC_2023/zheng_processed.qs')
colnames(hov.seu@meta.data)

table(hov.seu$maintypes_2)
rm(hov.seu)

hov.mac = subset(hov.seu,subset = maintypes_2 == 'Macrophage')
qsave(hov.mac,file= '/home/cwt/Project/Zheng_NC_2023/zheng_mac.qs')

hov.tnk = subset(hov.seu,idents == c('CD4+ T','CD8+ T','NK'))
qsave(hov.tnk,file= '/home/cwt/Project/Zheng_NC_2023/zheng_tnk.qs')


##cluster simiilarity of tnk, and macrophage
mac = qread('1_clustering/mac1118.qs')
tnk = qread('1_clustering/tnk1118.qs')

mac.marker= wt_downsamplefindmarker_human(mac)
tnk.marker = wt_downsamplefindmarker_human(tnk)

writexl::write_xlsx(mac.marker,'macmarker.xlsx')
writexl::write_xlsx(tnk.marker,'tnkmarker.xlsx')

table(hov.mac$Annotation)
table(hov.tnk$Annotation)

hov.mac@active.ident = factor(hov.mac$Annotation)
hov.tnk@active.ident = factor(hov.tnk$Annotation)

hov.mac.marker = wt_downsamplefindmarker_human(hov.mac)
hov.tnk.marker = wt_downsamplefindmarker_human(hov.tnk)

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
  
  
  DE_RNA = mac.marker[[2]]
  genelist = DE_RNA$gene
  genetrans = homologene(genelist,inTax = 10090,outTax = 9606)
  DE_RNA_merge = merge(DE_RNA,genetrans,by.x = 'gene',by.y='10090')
  DE_RNA_merge = DE_RNA_merge %>% dplyr::filter(avg_log2FC>0.585 & p_val_adj <0.05)
  
  #zheng_nature_cancer_comprison

  zheng = hov.mac.marker[[2]]
  zheng = zheng %>% dplyr::filter(avg_log2FC>0.585 & p_val_adj <0.05)
  cluster.zheng = unique(zheng$cluster)
  zheng_list = filter_data(zheng,clustername = cluster.zheng)
  MEPP_tumor_list = filter_data_m(DE_RNA_merge, clustername = unique(DE_RNA$cluster))
  
  p1 = wt_cluster_comparison(zheng_list,MEPP_tumor_list);p1##mac
  p2 = wt_cluster_comparison(zheng_list,MEPP_tumor_list);p2##tnk
  
  
  cowplot::plot_grid(plotlist = list(p1,p2),ncol = 2)
  ggsave(filename = 'zheng_comparison_immune.pdf',width = 11.85,height = 3.3)
}

##ROIE of clusters of tnk and macrophages !!!only in primary and metastatic tumors

{
  zheng.meta = hov.tnk@meta.data %>% dplyr::filter(Groups %in% c('Primary Tumor','Metastatic Tumor'))
  zheng.meta$Groups = droplevels(zheng.meta$Groups)
  summary = table(zheng.meta[,c('Annotation','Groups')])
  roe = as.data.frame(ROIE(summary))
  roe$cluster = rownames(roe)
  rownames(roe) = NULL
  res = data.frame()
  res = rbind(res,roe)
  summary(roe)
  
  heatmapdata = gather(roe,group,value,-c(cluster))
  str(heatmapdata)
  
  heatmapdata$group = factor(heatmapdata$group,levels = c('Primary Tumor','Metastatic Tumor'))
  heatmapdata$cluster = factor(heatmapdata$cluster,levels = unique(hov.tnk$Annotation))
  
  
  p_mac = ggplot(heatmapdata,aes(x = group,y=cluster,fill = value))+
    geom_tile()+scale_fill_gradient(high = "#FF8787",low = "#FFF7EC")+
    geom_tile()+
    geom_text(aes(label = round(value,2)))+theme_transparent()+
    theme(axis.text.x = element_text(angle = 45,vjust = 0.5))+ylab("")+xlab("")+
    theme(axis.text.y = element_text(face = "bold"))
  
  p_tnk = ggplot(heatmapdata,aes(x = group,y=cluster,fill = value))+
    geom_tile()+scale_fill_gradient(high = "#FF8787",low = "#FFF7EC")+
    geom_tile()+
    geom_text(aes(label = round(value,2)))+theme_transparent()+
    theme(axis.text.x = element_text(angle = 45,vjust = 0.5))+ylab("")+xlab("")+
    theme(axis.text.y = element_text(face = "bold"))
  
  cowplot::plot_grid(plotlist = list(p_mac,p_tnk),ncol = 2)
  ggsave(filename = "zheng_tnk_mac_roie.pdf",width = 8,height = 4)
}




``````1126 analysis for DL````````
hov.seu
hov.harmony

FeaturePlot(hov.harmony,features = c('METTL5'))
DotPlot(hov.harmony,features = c('METTL5'))
ggsave('mettl5.pdf',width = 4,height = 3.5)
VlnPlot(hov.harmony,features = c('METTL5'),group.by = 'Groups',pt.size = 0)+
  stat_compare_means()







