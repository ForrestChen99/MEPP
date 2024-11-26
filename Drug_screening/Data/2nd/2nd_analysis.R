#DRUG SCREENING 2ND

setwd("/Users/forrest/Desktop/Rtest/Drug_screening")

library(tidyverse)

select = dplyr::select
filter = dplyr::filter

library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ComplexHeatmap)



###MOSE / S378 / F378
mose_2nd = read.csv('2nd/mose_2nd.csv')
vector <- as.vector(t(mose_2nd))
print(vector)

s378_data = read.csv('2nd/378s_2nd.csv')
s378_via = as.vector(t(s378_data))
print(s378_via)

f378_data = read.csv('2nd/378f_2nd.csv')
f378_via = as.vector(t(f378_data))
print(f378_via)




###MAP DATA
map_2nd = read.csv('2nd/Drugmap_2nd.csv')

map_ID = map_2nd[seq(1,nrow(map_2nd),2),]
ID = as.vector(t(map_ID))
print(ID)

map_drug = map_2nd[seq(2,nrow(map_2nd),2),]
drug = as.vector(t(map_drug))
print(drug)


data_2nd = data.frame(ID = ID,
                      Drug = drug,
                      mose = vector,
                      s378 = s378_via,
                      f378 = f378_via)

data_2nd = data_2nd %>% filter_all(all_vars(. != 'EMPTY'))



```````````````Analysis`````````````````````````````````````````````
# Define thresholds
low_threshold <- 0.75
fold_change_threshold <- 1.5
# Function to determine specificity
is_specific <- function(low, other1, other2, low_threshold, fold_change_threshold) {
  if (low < low_threshold && other1 > (low * fold_change_threshold) && other2 > (low * fold_change_threshold)) {
    return("Specific")
  } else {
    return("Not Specific")
  }
}

# Identify specific drugs for each cell line
df = data_2nd
specific_drugs <- df %>%
  mutate(
    SpecificToMOSE = mapply(is_specific, mose, s378, f378, MoreArgs = list(low_threshold = low_threshold, fold_change_threshold = fold_change_threshold)),
    SpecificTos378 = mapply(is_specific, s378, mose, f378, MoreArgs = list(low_threshold = low_threshold, fold_change_threshold = fold_change_threshold)),
    SpecificTof378 = mapply(is_specific, f378, s378, mose, MoreArgs = list(low_threshold = low_threshold, fold_change_threshold = fold_change_threshold))
  )

arranged_drugs <- specific_drugs %>%
  mutate(SpecificOrder = case_when(
    SpecificToMOSE == "Specific" ~ 1,
    SpecificTos378 == "Specific" ~ 2,
    SpecificTof378 == "Specific" ~ 3,
    TRUE ~ 4
  )) %>%
  arrange(SpecificOrder, Drug) %>%
  select(-SpecificOrder)  # Remove the helper column used for sorting

# Print the arranged dataframe
print(arranged_drugs)

write.csv(arranged_drugs,'viability_2nd.csv')

````````````````````````````````````````````````````````


##Complexhetamap
heatmap.df = arranged_drugs %>% filter(SpecificToMOSE == 'Specific' | SpecificTos378== 'Specific' |SpecificTof378 == 'Specific') %>% select(Drug,mose,s378,f378)
rownames(heatmap.df) = NULL
heatmap.df = column_to_rownames(heatmap.df,'Drug')

mat = as.matrix(heatmap.df)
mat_scaled = t(apply(mat, 1, scale))
colnames(mat_scaled) <- colnames(heatmap.df)



# Create a vector for the annotation
heatmap.df1 = arranged_drugs %>% filter(SpecificToMOSE == 'Specific' | SpecificTos378== 'Specific' |SpecificTof378 == 'Specific')
annotation_vector <- ifelse(heatmap.df1$SpecificToMOSE == "Specific", "mose",
                            ifelse(heatmap.df1$SpecificTos378 == "Specific", "s378",
                                   ifelse(heatmap.df1$SpecificTof378 == "Specific", "f378", "Not Specific")))

# Create the annotation
left_annotation <- rowAnnotation(
  Specificity = annotation_vector,
  col = list(Specificity = c("mose" = "#BE2490", "s378" = "#007AB5", "f378" = "#FF9999", "Not Specific" = "grey"))
)


# Create the heatmap with the left annotation
cell_size <- unit(5, "mm")  # Adjust the size as needed

library(circlize)
heatmap <- Heatmap(mat_scaled,
                   name = "Scaled Viability",
                   col = colorRamp2(c(-1.5, 0, 1.5), c("#1679AB", "white", "#C80036")), # Adjust the color scale as needed
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
draw(heatmap)
pdf('heatmap_2nd.pdf',width = 7,height = 13)
draw(heatmap)
dev.off()










































