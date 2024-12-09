wt_downsamplefindmarker_mouse = function(data){
  downsample = subset(data,downsample=500)
  markers = FindAllMarkers(downsample,logfc.threshold = 0.25,min.pct = 0.1,only.pos = T)
  markers = markers[!grepl("^Rp[sl]",markers$gene,ignore.case = F),]
  markers = markers[!grepl("^mt-",markers$gene,ignore.case = F),]
  top_10 = markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
  list.marker = list(top_10,markers)
  return(list.marker)
}


wt_downsamplefindmarker_human = function(data){
  downsample = subset(data,downsample=500)
  markers = FindAllMarkers(downsample,logfc.threshold = 0.25,min.pct = 0.1,only.pos = T)
  markers = markers[!grepl("^RP[SL]",markers$gene,ignore.case = F),]
  markers = markers[!grepl("^MT-",markers$gene,ignore.case = F),]
  top_10 = markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
  list.marker = list(top_10,markers)
  return(list.marker)
}


wt_firststep = function(sub1){
  DefaultAssay(sub1)="RNA"
  sub1  <- NormalizeData(sub1) %>% FindVariableFeatures(nfeatures = 2000)
  VGENES=VariableFeatures(sub1)
  VGENES=setdiff(VGENES,VGENES[grep("^mt|^Rpl|^Rps",VGENES)])
  sub1 <- ScaleData(sub1,features =VGENES,vars.to.regress =c("percent.mito","nCount_RNA","nFeature_RNA","percent.ribo")) %>% RunPCA(npcs = 20, verbose = FALSE)
  sub1 <- sub1 %>% 
    RunHarmony("orig.ident", plot_convergence = TRUE)
  harmony_embeddings <- Embeddings(sub1, 'harmony')
  sub1 <- sub1 %>% 
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20)
  for(i in seq(0.1,0.5,by = 0.1)){
    sub1 = FindClusters(sub1,resolution = i,algorithm = 1)
  }
  return(sub1)
}
