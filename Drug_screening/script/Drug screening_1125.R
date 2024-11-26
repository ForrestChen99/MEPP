getwd()
setwd("/Users/forrest/Desktop/Rtest/Drug_screening")
####241125
#nature compounds 里面删掉kras重新进行分析
#实际上thioguanine和fk228都是在epi里面得到的，所以对nature compounds里面的没有影响。重新做nature compounds的图和表就好了



library(tidyverse)

select = dplyr::select
filter = dplyr::filter

library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
##MOSE drug screening

df = read.csv('MOSE_RAW.csv')
vector <- as.vector(t(df))
print(vector)


map = read.csv('MAP_19.csv')
map_ID = map[seq(1,nrow(map),2),]
ID = as.vector(t(map_ID))
print(ID)

map_drug = map[seq(2,nrow(map),2),]
drug = as.vector(t(map_drug))
print(drug)


MOSE = data.frame(ID = ID,
                  Drug = drug,
                  MOSE = vector)

value_hpj = read.csv('HPJ_378.csv')
hpj = as.vector(t(value_hpj))
MOSE$M378 = hpj

MOSE <- MOSE %>% mutate_all(~replace(., . == 999, NA))

##delete row of NA value
MOSE = na.omit(MOSE)
MOSE = MOSE %>% filter_all(all_vars(. != 'Empty'))


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
df = MOSE
specific_drugs <- df %>%
  mutate(
    SpecificToMOSE = mapply(is_specific, MOSE, M378, MoreArgs = list(low_threshold = low_threshold, fold_change_threshold = fold_change_threshold)),
    SpecificToM378 = mapply(is_specific, M378, MOSE, MoreArgs = list(low_threshold = low_threshold, fold_change_threshold = fold_change_threshold))
  )

write.csv(specific_drugs,'241125/nature_specific_0.75_1.5_epi.csv')

arranged_drugs_nat <- specific_drugs %>%
  mutate(SpecificOrder = case_when(
    SpecificToMOSE == "Specific" ~ 1,
    SpecificToM378 == "Specific" ~ 2,
    TRUE ~ 3
  )) %>%
  arrange(SpecificOrder, Drug) %>%
  select(-SpecificOrder)  # Remove the helper column used for sorting



{
  ##Complexhetamap
  heatmap.df = arranged_drugs_nat %>% filter(SpecificToMOSE == 'Specific' | SpecificToM378 == 'Specific') %>% select(Drug,MOSE,M378)
  rownames(heatmap.df) = NULL
  heatmap.df = column_to_rownames(heatmap.df,'Drug')
  
  mat = as.matrix(heatmap.df)
  mat_scaled = apply(mat, 2, scale)# 2 scale by column
  colnames(mat_scaled) <- colnames(heatmap.df)
  rownames(mat_scaled) <- rownames(heatmap.df)
  
  
  # Create a vector for the annotation
  heatmap.df1 = arranged_drugs_nat %>% filter(SpecificToMOSE == 'Specific'|SpecificToM378 == 'Specific')
  annotation_vector <- ifelse(heatmap.df1$SpecificToMOSE == "Specific", "MOSE",
                              ifelse(heatmap.df1$SpecificToM378 == "Specific", "M378", "Not Specific"))
  
  # Create the annotation
  left_annotation <- rowAnnotation(
    Specificity = annotation_vector,
    col = list(Specificity = c("MOSE" = "#BE2490", "M378" = "#FF9999", "Not Specific" = "grey"))
  )
  
  
  # Create the heatmap with the left annotation
  cell_size <- unit(5, "mm")  # Adjust the size as needed
  
  heatmap <- Heatmap(mat_scaled,
                     name = "Viability",
                     col = colorRamp2(c(0, 0.5, 1), c("#1679AB", "white", "#C80036")), # Adjust the color scale as needed
                     row_names_side = "left",
                     column_names_side = "top",
                     show_row_dend = F,
                     show_column_dend = F,
                     left_annotation = left_annotation,
                     row_gap = unit(1, "mm"),
                     column_gap = unit(1, "mm"),
                     cluster_columns = F,cluster_rows = F,
                     width = cell_size * ncol(mat_scaled),
                     height = cell_size * nrow(mat_scaled))
  
  # Draw the heatmap
  pdf('heatmap_nat.pdf',width = 7,height = 13)
  draw(heatmap)
  dev.off()
  
  
}







````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````








