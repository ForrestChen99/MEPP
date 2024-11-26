getwd()
setwd("/Users/forrest/Desktop/Rtest/Drug_screening")

library(tidyverse)

select = dplyr::select
filter = dplyr::filter

library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ComplexHeatmap)

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
                Viability = vector)

value_stomach = read.csv('KRAS_STOMACH.csv')
value = as.vector(t(value_stomach))
print(value)

value_hpj = read.csv('HPJ_378.csv')
hpj = as.vector(t(value_hpj))
MOSE$KRAS = value
MOSE$HPJ = hpj

MOSE <- MOSE %>% mutate_all(~replace(., . == 999, NA))

##delete row of NA value
MOSE = na.omit(MOSE)
MOSE = MOSE %>% filter_all(all_vars(. != 'Empty'))



MOSE$FC = MOSE$KRAS/MOSE$Viability
MOSE$FC1 = MOSE$HPJ/MOSE$Viability
MOSE$FC2 = MOSE$Viability/MOSE$KRAS
MOSE$FC3 = MOSE$HPJ/MOSE$KRAS
MOSE$FC4 = MOSE$Viability/MOSE$HPJ
MOSE$FC5 = MOSE$KRAS/MOSE$HPJ

##analysis

threshold = 0.75
summary(MOSE)

mose1 = MOSE %>% filter(Viability < threshold)


filter1 = MOSE %>% filter(Viability < threshold & FC >2) %>% select(Drug)
filter2 =  MOSE %>% filter(Viability < threshold & FC >1.5) %>% select(Drug)
filter3 = MOSE %>% filter(Viability < threshold & FC >1.2) %>% select(Drug)
filter4 = MOSE %>% filter(Viability < threshold & FC >1) %>% select(Drug)


dt_plot = MOSE %>% arrange(Viability)
dt_plot = dt_plot %>% rownames_to_column('Rank')


dt_plot$Rank = as.numeric(dt_plot$Rank)
dt_plot$Label = ifelse(dt_plot$FC>2,'FC>2','FC≤2')

#write.csv(dt_plot,'MOSE_DATA.csv')


ggplot(dt_plot, aes(x = order(Rank), y = Viability, color = FC)) +
  geom_point(aes(color= Label)) + 
  scale_color_manual(values = c("FC≤2" = "grey", "FC>2" = "red")) +
  xlab("Rank") + 
  ylab("Viability") + 
  ggtitle("Drug Viability: MOSE vs KRAS") +
  theme_minimal()


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
df = MOSE
specific_drugs <- df %>%
  mutate(
    SpecificToMOSE = mapply(is_specific, MOSE, KRAS, M378, MoreArgs = list(low_threshold = low_threshold, fold_change_threshold = fold_change_threshold)),
    SpecificToKRAS = mapply(is_specific, KRAS, MOSE, M378, MoreArgs = list(low_threshold = low_threshold, fold_change_threshold = fold_change_threshold)),
    SpecificToM378 = mapply(is_specific, M378, MOSE, KRAS, MoreArgs = list(low_threshold = low_threshold, fold_change_threshold = fold_change_threshold))
  )

# Print the dataframe with specificity information
print(specific_drugs)

specific_drugs = specific_drugs %>% select(-c('FC','FC1','FC2','FC3','FC4','FC5'))
write.csv(specific_drugs,'specific_0.75_1.5.csv')

arranged_drugs <- specific_drugs %>%
  mutate(SpecificOrder = case_when(
    SpecificToMOSE == "Specific" ~ 1,
    SpecificToKRAS == "Specific" ~ 2,
    SpecificToM378 == "Specific" ~ 3,
    TRUE ~ 4
  )) %>%
  arrange(SpecificOrder, Drug) %>%
  select(-SpecificOrder)  # Remove the helper column used for sorting

# Print the arranged dataframe
print(arranged_drugs)




##Complexhetamap
heatmap.df = arranged_drugs %>% filter(SpecificToMOSE == 'Specific' | SpecificToKRAS== 'Specific' |SpecificToM378 == 'Specific') %>% select(Drug,MOSE,KRAS,M378)
rownames(heatmap.df) = NULL
heatmap.df = column_to_rownames(heatmap.df,'Drug')

mat = as.matrix(heatmap.df)
mat_scaled = t(apply(mat, 1, scale))
colnames(mat_scaled) <- colnames(heatmap.df)

Heatmap(mat_scaled,cluster_rows = F,cluster_columns = F,
        left_annotation = )


# Create a vector for the annotation
heatmap.df1 = arranged_drugs %>% filter(SpecificToMOSE == 'Specific' | SpecificToKRAS== 'Specific' |SpecificToM378 == 'Specific')
annotation_vector <- ifelse(heatmap.df1$SpecificToMOSE == "Specific", "MOSE",
                            ifelse(heatmap.df1$SpecificToKRAS == "Specific", "KRAS",
                                   ifelse(heatmap.df1$SpecificToM378 == "Specific", "M378", "Not Specific")))

# Create the annotation
left_annotation <- rowAnnotation(
  Specificity = annotation_vector,
  col = list(Specificity = c("MOSE" = "#BE2490", "KRAS" = "#007AB5", "M378" = "#FF9999", "Not Specific" = "grey"))
)


# Create the heatmap with the left annotation
cell_size <- unit(5, "mm")  # Adjust the size as needed

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

# Draw the heatmap
pdf('heatmap1.pdf',width = 7,height = 13)
draw(heatmap)
dev.off()


````````````````````````````Epigenetics``````````````````````````````````````````````````````````````
epi_map = read.csv('epi_layout.csv')
mose_epi= read.csv('mose_epi.csv')
m378_epi = read.csv('M378_epi.csv')

mose.epi <- as.vector(t(mose_epi))
m378.epi = as.vector(t(m378_epi))



epi_map1 = epi_map[seq(2,nrow(epi_map),2),]
drug = as.vector(t(epi_map1))
print(drug)

epi_ID = epi_map[seq(1,nrow(epi_map),2),]
ID = as.vector(t(epi_ID))
print(ID)


EPI = data.frame(ID = ID,
                  Drug = drug,
                  MOSE = mose.epi,
                 M378 = m378.epi)

EPI = EPI %>% filter_all(all_vars(. != 'Empty'))


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
df = EPI
specific_drugs_epi <- df %>%
  mutate(
    SpecificToMOSE = mapply(is_specific, MOSE, M378, MoreArgs = list(low_threshold = low_threshold, fold_change_threshold = fold_change_threshold)),
    SpecificToM378 = mapply(is_specific, M378, MOSE, MoreArgs = list(low_threshold = low_threshold, fold_change_threshold = fold_change_threshold))
  )

write.csv(specific_drugs,'specific_0.75_1.5_epi.csv')

arranged_drugs_epi <- specific_drugs_epi %>%
  mutate(SpecificOrder = case_when(
    SpecificToMOSE == "Specific" ~ 1,
    SpecificToM378 == "Specific" ~ 2,
    TRUE ~ 3
  )) %>%
  arrange(SpecificOrder, Drug) %>%
  select(-SpecificOrder)  # Remove the helper column used for sorting

write.csv(arranged_drugs_epi,'specific_0.75_1.5_epi.csv')


{
  ##Complexhetamap
  heatmap.df = arranged_drugs_epi %>% filter(SpecificToMOSE == 'Specific' | SpecificToM378 == 'Specific') %>% select(Drug,MOSE,M378)
  rownames(heatmap.df) = NULL
  heatmap.df = column_to_rownames(heatmap.df,'Drug')
  
  mat = as.matrix(heatmap.df)
  #mat_scaled = apply(mat, 2, scale)# 2 scale by column
  #colnames(mat_scaled) <- colnames(heatmap.df)
  

  # Create a vector for the annotation
  heatmap.df1 = arranged_drugs_epi %>% filter(SpecificToMOSE == 'Specific'|SpecificToM378 == 'Specific')
  annotation_vector <- ifelse(heatmap.df1$SpecificToMOSE == "Specific", "MOSE",
                              ifelse(heatmap.df1$SpecificToM378 == "Specific", "M378", "Not Specific"))
  
  # Create the annotation
  left_annotation <- rowAnnotation(
    Specificity = annotation_vector,
    col = list(Specificity = c("MOSE" = "#BE2490", "M378" = "#FF9999", "Not Specific" = "grey"))
  )
  
  
  # Create the heatmap with the left annotation
  cell_size <- unit(5, "mm")  # Adjust the size as needed
  
  heatmap <- Heatmap(mat,
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
  pdf('heatmap_epi.pdf',width = 7,height = 13)
  draw(heatmap)
  dev.off()
  
  
}










