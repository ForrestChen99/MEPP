##generate supplementary figures?

```````````````````````wholedata`````````````````````````````````````

wholedata = qread('1_clustering/wholedata1118.qs')
wholedata

key.features = c('Wt1','Inha','C1qc','Clec9a','S100a8','Tpsab1','Cd3g','Cd79a','Igfbp7','Rhag')
FeaturePlot(wholedata,features = key.features,raster = F,ncol = 5,order= T)
ggsave('wholedata_feature.pdf',width = 12.3,height = 4.5)

table(wholedata$group)


```````````````````````epi`````````````````````````````````````
#Heatmap of the epi
{
  selected.epi.marker = epi.marker[[2]] %>% group_by(cluster) %>% top_n(10,wt=avg_log2FC)
  {
    DefaultAssay(epi)
    gene = selected.epi.marker %>% arrange(cluster) %>% pull(gene)
    
    gene_cell_exp = AggregateExpression(epi,features = gene,group.by = 'dcelltype')
    df <- data.frame(colnames(as.data.frame(gene_cell_exp$RNA)))
    colnames(df) <- 'class'
    color.use = epi.color
    celltype = unique(epi$dcelltype)
    celltype = as.factor(celltype)
    
    marker_exp <- t(scale(t(as.matrix(gene_cell_exp$RNA)),scale = T))
    marker_exp <- marker_exp[match(gene, rownames(marker_exp)), ]
    
    require(circlize)
    col_fun = colorRamp2(c(-1,0,1.5),c("#1679AB", "white", "#C80036"))
    
    pdf("epi_heatmap_long.pdf",width = 3.9,height = 6.8)
    cell_size <- unit(5, "mm") 
    Heatmap(marker_exp,
            name = "Z Score",
            cluster_rows = F,
            cluster_columns = F,
            show_row_names = T,
            show_column_names = T,
            col =col_fun,
            row_names_side = "left",
            column_names_side = "top")
            #row_gap = unit(1, "mm"),
            #column_gap = unit(1, "mm"),width = cell_size * ncol(marker_exp),
            #height = cell_size * nrow(marker_exp))
    dev.off()
  }
}



#featureplot of the epi
key.features = c('Fxyd3','Angptl7','Mgp','Fmo1','Nuf2','Mki67','Ero1l','Nos2')
FeaturePlot(epi,features = key.features,raster = F,ncol = 4,order= T)
ggsave('epi_featureplot.pdf',width = 11.4,height = 4.5)


#