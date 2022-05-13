#Author: Sriramteja Veerisetti 
#BF 528 
#Biologist Role in Project 3
#Individual Project 

library(tidyverse)
library(dplyr)

#'@details The purpose here is to filter the deseq .csv file so garner only 
#'significant genes 
deseq_Ahr <- read.csv("/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_3_Biologist_Role/DESeqAhRresults.csv")
colnames(deseq_Ahr)[1] <- 'gene'

deseq_CAR_PXR <- read.csv("/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_3_Biologist_Role/DESeqCAR-PXRresults.csv")
colnames(deseq_CAR_PXR )[1] <- 'gene'

deseq_PPARA <- read.csv("/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_3_Biologist_Role/DESeqPPARAresults.csv")
colnames(deseq_PPARA)[1] <- 'gene'


#'@details The Wang et. al paper filtered the DE genes based on p value that was less than 0.05.
#'as well as a fold change that is greater than 1.5. They considered these 
#'genes to be differentially expressed. Figure 4 of the paper shows this.  
deseq_Ahr_filtered <- deseq_Ahr %>%
  filter(pvalue < 0.05) %>%
  filter(log2FoldChange > 1.5)

deseq_CAR_PXR_filtered <- deseq_CAR_PXR %>%
  filter(pvalue < 0.05) %>%
  filter(log2FoldChange > 1.5)

deseq_PPARA_filtered <- deseq_PPARA %>%
  filter(pvalue < 0.05) %>%
  filter(log2FoldChange > 1.5)

#'Save the results to a .csv file 
deseq_Ahr_filtered_new <- write.csv (deseq_Ahr_filtered, "/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_3_Biologist_Role/deseq_Ahr_filtered.csv")
deseq_Ahr_filtered_z <- read_csv("/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_3_Biologist_Role/deseq_Ahr_filtered.csv")
deseq_Ahr_filtered_z <- deseq_Ahr_filtered_z[, -c(1)]
deseq_Ahr_filtered_z <- write_csv(deseq_Ahr_filtered_z, "/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_3_Biologist_Role/deseq_Ahr_filtered_z.csv")

deseq_CAR_PXR_filtered_new <- write.csv (deseq_CAR_PXR_filtered, "/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_3_Biologist_Role/deseq_CAR_PXR_filtered.csv")
deseq_CAR_PXR_filtered_z  <- read_csv("/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_3_Biologist_Role/deseq_CAR_PXR_filtered.csv")
deseq_CAR_PXR_filtered_z <- deseq_CAR_PXR_filtered_z[, -c(1)]
deseq_CAR_PXR_filtered_z <- write_csv(deseq_CAR_PXR_filtered_z, "/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_3_Biologist_Role/deseq_CAR_PXR_zfiltered.csv")
#deseq_CAR_PXR_filtered_z <- deseq_CAR_PXR_filtered_z[, -c(1)]
#deseq_CAR_PXR_filtered_z <- write_csv(deseq_CAR_PXR_filtered_z, "/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_3_Biologist_Role/deseq_CAR_PXR_zfiltered.csv")

deseq_PPARA_filtered_new <- write.csv (deseq_PPARA_filtered, "/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_3_Biologist_Role/deseq_PPARA_filtered.csv")
deseq_PPARA_filtered_z  <- read_csv("/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_3_Biologist_Role/deseq_PPARA_filtered.csv")
deseq_PPARA_filtered_z <- deseq_PPARA_filtered_z[, -c(1)]
deseq_PPARA_filtered_z <- write_csv(deseq_PPARA_filtered_z, "/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_3_Biologist_Role/deseq_PPARA_filtered_z_filtered.csv")


deseq_Ahr_filtered <- deseq_Ahr_filtered[, 2]
#deseq_CAR_PXR_filtered <- as.data.frame(deseq_CAR_PXR_filtered)
deseq_CAR_PXR_filtered <- deseq_CAR_PXR_filtered[, 2]
deseq_PPARA_filtered <- deseq_PPARA_filtered[, 2]
deseq_combo_list <- list(deseq_Ahr_filtered, deseq_CAR_PXR_filtered, deseq_PPARA_filtered)
deseq_combo <- deseq_combo_list %>% reduce(full_join, by = 'gene')
deseq_combo <- as.matrix(deseq_combo)
deseq_combo <- as.data.frame(deseq_combo)

deseq_combo <- write_csv(deseq_combo, "/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_3_Biologist_Role/deseq_combo.csv" )



#'@details Here we want to filter the counts matrix. Our goal here is to 
#'add an average count column and coefficient of variation matrix   
AhR_deseq_ncounts <- read.csv("/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_3_Biologist_Role/AhR_deseq_norm_counts.csv")
colnames(AhR_deseq_ncounts)[1] <- 'gene'

CARPXR_deseq_ncounts <- read.csv("/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_3_Biologist_Role/CARPXR_deseq_norm_counts.csv")
colnames(CARPXR_deseq_ncounts)[1] <- 'gene'

PPARA_deseq_ncounts <- read.csv("/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_3_Biologist_Role/PPARA_deseq_norm_counts.csv")
colnames(PPARA_deseq_ncounts)[1] <- 'gene'

df_list <- list(AhR_deseq_ncounts, CARPXR_deseq_ncounts, PPARA_deseq_ncounts)

combination <- df_list %>% reduce(full_join, by = 'gene')

#'Append averages to a separate column
combination$gene <- NULL
combination$average_ncounts <- rowMeans(subset(combination, select = 2:6))
sdeviation <- apply(combination, 1, sd)
combination$sdeviation <- sdeviation
combination <- mutate(combination, coefficient_of_variation = sdeviation/average_ncounts)


#' @details Usually coefficient of variation (cv) values that are less than 1 are
#'significant because they are considered to be low-variance. Low variance 
#'meaning that they do not disperse from the variation much. 
combination <- combination %>%
  filter(coefficient_of_variation < 1)

combination <- combination %>%
  filter(average_ncounts > 100)

#' @details We need to convert the dataframe into a matrix in order to 
#' develop a heatmap. First we will have to get rid of the columns of the 
#' df that deal with the cv, sd, and average counts. 

combination <- combination[, -c(19:21)]
combination_matrix <- as.matrix(combination)
deseq_filtered_heatmap <- heatmap(combination_matrix, cexRow = 0.7, cexCol = 0.7, 
                                  col = terrain.colors(256))

deseq_filtered_heatmap
