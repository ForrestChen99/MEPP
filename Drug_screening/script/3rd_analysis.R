
````````````````````````````3rd analysis``````````````````````````````````````````````````````````````
setwd("/Users/forrest/Desktop/Rtest/Drug_screening/")
library(tidyverse)

select = dplyr::select
filter = dplyr::filter

library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ComplexHeatmap)


````````````````````````````Data import``````````````````````````````````````````````````````````````
list.files('3rd/processed/')
filepath = list.files('3rd/processed',full.names = T)
file = list.files('3rd/processed')

datalist = list()
data_import <- function(filepath) {
  data <- read.csv(filepath)
  vector <- as.vector(t(data))
  return(vector)
}

datalist = lapply(filepath,data_import)

##import map data
map_2nd = read.csv('2nd/Drugmap_2nd.csv')

map_ID = map_2nd[seq(1,nrow(map_2nd),2),]
ID = as.vector(t(map_ID))
print(ID)

map_drug = map_2nd[seq(2,nrow(map_2nd),2),]
drug = as.vector(t(map_drug))
print(drug)

data_3rd = data.frame(ID = ID,
                      Drug = drug)
                  
column_names <- gsub("\\.csv$", "", basename(file))
for (i in seq_along(datalist)) {
  data_3rd[[column_names[i]]] <- datalist[[i]]
}

data_3rd = data_3rd %>% filter_all(all_vars(. != 'EMPTY'))



```````````````Analysis: mouse cell line`````````````````````````````````````````````

data_mouse = data_3rd %>% select(ID,Drug,ID8,MOSE)

# Define thresholds
low_threshold <- 0.75
fold_change_threshold <- 1.5

# Function to determine specificity
is_specific <- function(low, other, low_threshold, fold_change_threshold) {
  if (low < low_threshold && other > (low * fold_change_threshold)) {
    return("Specific")
  } else {
    return("Not Specific")
  }
}

# Identify specific drugs for each cell line
specific_drugs_mosue <- data_mouse %>%
  mutate(
    SpecificToMOSE = mapply(is_specific, MOSE, ID8, MoreArgs = list(low_threshold = low_threshold, fold_change_threshold = fold_change_threshold)),
    SpecificToID8 = mapply(is_specific, ID8, MOSE, MoreArgs = list(low_threshold = low_threshold, fold_change_threshold = fold_change_threshold))
  )

mose.specific.drug = specific_drugs_mosue %>% filter(SpecificToMOSE == 'Specific') %>% pull(Drug)


arranged_drugs <- specific_drugs_mosue %>%
  mutate(SpecificOrder = case_when(
    SpecificToMOSE == "Specific" ~ 1,
    SpecificToID8 == "Specific" ~ 2,
    TRUE ~ 3
  )) %>%
  arrange(SpecificOrder, Drug) %>%
  select(-SpecificOrder)  # Remove the helper column used for sorting

# Print the arranged dataframe
print(arranged_drugs)

write.csv(arranged_drugs,'3rd/viability_3rd.csv')





```````````````Analysis: human cell line`````````````````````````````````````````````
#if the drug is also useful in human cell line
data_human = data_3rd %>% select(-ID8,-MOSE)

human.dt = data_human %>% filter(Drug %in% mose.specific.drug)



#Complexheatmap
{
  heatmap.df = merge(arranged_drugs,data_human,by = 'Drug',all.x = T,sort = F)
  
  heatmap.df = heatmap.df %>% filter(SpecificToMOSE == 'Specific' | SpecificToID8== 'Specific') %>% select(Drug,MOSE,ID8,column_names)
  rownames(heatmap.df) = NULL
  heatmap.df = heatmap.df %>% distinct(Drug,.keep_all = T)
  heatmap.df = column_to_rownames(heatmap.df,'Drug')
  
  
 
  
  mat = as.matrix(heatmap.df)
  mat_scaled = t(apply(mat, 1, scale))
  colnames(mat_scaled) <- colnames(heatmap.df)
  
  
  
  # Create a vector for the annotation
  arranged_drugs =  arranged_drugs %>% distinct(Drug,.keep_all = T)
  
  heatmap.df1 = arranged_drugs %>% filter(SpecificToMOSE == 'Specific' | SpecificToID8== 'Specific')
  annotation_vector <- ifelse(heatmap.df1$SpecificToMOSE == "Specific", "MOSE",
                              ifelse(heatmap.df1$SpecificToID8 == "Specific", "ID8",'Not Specific'))
  
  # Create the annotation
  left_annotation <- rowAnnotation(
    Specificity = annotation_vector,
    col = list(Specificity = c("MOSE" = "#FF9874", "ID8" = "#7C93C3"))
  )
  
  
  cell.color = c('#C0C78C','#FF8A8A','#AD49E1','#A1D6B2','#D1E9F6','#D0B8A8',
                 '#606676','#55AD9B')
  names(cell.color) = colnames(heatmap.df)
  
  Species.color = c('#9BB0C1','#7469B6')
  names(Species.color) = c('Mouse','Human')
  
  
  top_annotation = HeatmapAnnotation(
    Cell = factor(colnames(heatmap.df),levels = colnames(heatmap.df)),
    Species = c(rep('Mouse',2),rep('Human',6)),
    col = list(Cell = cell.color,Species = Species.color)
  )
  
  # Create the heatmap with the left annotation
  cell_size <- unit(5, "mm")  # Adjust the size as needed
  
  library(circlize)
  heatmap <- Heatmap(mat_scaled,
                     name = "Scaled Viability",
                     col = colorRamp2(c(-2, 0, 2), c("#1679AB", "white", "#C80036")), # Adjust the color scale as needed
                     row_names_side = "left",
                     column_names_side = "top",
                     show_row_dend = F,
                     show_column_dend = F,
                     left_annotation = left_annotation,
                     top_annotation = top_annotation,
                     row_gap = unit(1, "mm"),
                     column_gap = unit(1, "mm"),
                     cluster_columns = F,cluster_rows = F,
                     width = cell_size * ncol(mat_scaled),
                     height = cell_size * nrow(mat_scaled))
  draw(heatmap)
  pdf('heatmap_3nd.pdf',width = 11.008,height = 9)
  draw(heatmap)
  dev.off()
  
}







