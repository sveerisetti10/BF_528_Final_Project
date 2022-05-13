#Author: Sriramteja Veerisetti 
#BF 528 
#Analyst Role in Project 3
#Individual Project 

#'Note from Sri: Make sure to run install the limma package if you do not have 
#'on your computer! 
#BiocManager::install("limma")
library(tidyverse)
library(dplyr)
library(limma)
library(reshape2)


############## Microarray Differential Expression with Limma ##################

#'@details We can utilize the next few steps for each chemical we decide to 
#'utilize. For example, we will only be using the chemical: 3-METHYLCHOLANTHRENE
#'for the next few steps. 

data <- read_csv('/project/bf528/project_3/groups/group_6_mic_info.csv')

#Here we need to load in the full RNA normalized matrix from all of the experiments
rma <- read.table('/project/bf528/project_3/samples/liver-normalization-rma.txt',
                  sep='\t',
                  as.is=TRUE,
                  header=TRUE,
                  row.names=1,)

#'We will start by designing a model matrix so that we can utilize 
#'it in the lmfit() function 
design <- model.matrix(
  ~factor(
    data$chemical, 
    level = c('Control', '3-METHYLCHOLANTHRENE')
  )
)

colnames(design) <- c('chemical', '3-METHYLCHOLANTHRENE')

#'Here we will only ask for the data within the data matrix to pertain to 
#'control as well as the 3-METHYLCHOLANTHRENE chemical. 
rma.specificiations <- rma[paste0('X', data$array_id[data$chemical == 'Control' | data$chemical == '3-METHYLCHOLANTHRENE'])]

#Here we want to run the limma function. This is the basic outline on how to do this
fit_limma <- lmFit(rma.specificiations, design)
fit_limma <- eBayes(fit_limma)
t_limma <- topTable(fit_limma, coef=2, n=nrow(rma.specificiations), adjust='BH')
t_limma <- tibble::rownames_to_column(t_limma, "Gene")
colnames(t_limma)[6] <- "padj"
t_limma <- t_limma[with(t_limma, order(-padj)), ]


#Write out the results to a csv file 
t_limma_written <- write_csv(t_limma, "/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_1_Analyst_Role/limma_3_METHYLCHOLANTHRENE_.csv")


#'@details We will only be using the chemical: FLUCONAZOLE
#'for the next few steps. 

data <- read_csv('/project/bf528/project_3/groups/group_6_mic_info.csv')

#Here we need to load in the full RNA normalized matrix from all of the experiments
rma <- read.table('/project/bf528/project_3/samples/liver-normalization-rma.txt',
                  sep='\t',
                  as.is=TRUE,
                  header=TRUE,
                  row.names=1,)

design <- model.matrix(
  ~factor(
    data$chemical, 
    level = c('Control', 'FLUCONAZOLE')
  )
)

colnames(design) <- c('Intercept', 'FLUCONAZOLE')
rma.specificiations <- rma[paste0('X', data$array_id[data$chemical == 'Control' | data$chemical == 'FLUCONAZOLE'])]

fit_limma <- lmFit(rma.specificiations, design)
fit_limma <- eBayes(fit_limma)
t_limma <- topTable(fit_limma, coef=2, n=nrow(rma.specificiations), adjust='BH')
t_limma <- tibble::rownames_to_column(t_limma, "Gene")
colnames(t_limma)[6] <- "padj"
t_limma <- t_limma[with(t_limma, order(-padj)), ]

#Write out the results to a csv file 
t_limma_written <- write_csv(t_limma, "/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_1_Analyst_Role/limma_FLUCONAZOLE_.csv")



#'@details We will only be using the chemical: PIRINIXIC_ACID
#'for the next few steps. 

data <- read_csv('/project/bf528/project_3/groups/group_6_mic_info.csv')

#Here we need to load in the full RNA normalized matrix from all of the experiments
rma <- read.table('/project/bf528/project_3/samples/liver-normalization-rma.txt',
                  sep='\t',
                  as.is=TRUE,
                  header=TRUE,
                  row.names=1,)

design <- model.matrix(
  ~factor(
    data$chemical, 
    level = c('Control', 'PIRINIXIC_ACID')
  )
)

colnames(design) <- c('Intercept', 'PIRINIXIC_ACID')
rma.specificiations <- rma[paste0('X', data$array_id[data$chemical == 'Control' | data$chemical == 'PIRINIXIC_ACID'])]

fit_limma <- lmFit(rma.specificiations, design)
fit_limma <- eBayes(fit_limma)
t_limma <- topTable(fit_limma, coef=2, n=nrow(rma.specificiations), adjust='BH')
t_limma <- tibble::rownames_to_column(t_limma, "Gene")


#'@details next we need to find the top ten DE genes from each of the analyses 
#'that were conducted. We will order the padj column so that we get the lowest 
#'padj values at the top because the smaller the padj value, the more significant 
#'the results are. 

colnames(t_limma)[6] <- "padj"
t_limma <- t_limma[with(t_limma, order(-padj)), ]


#Write out the results to a csv file 
t_limma_written <- write_csv(t_limma, "/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_1_Analyst_Role/limma_PIRINIXIC_ACID_.csv")


#'@details Now that we have the differential expression results from the 
#'Limma process, we can now filter the files so that significant genes with a 
#'p-adj value < 0.5 are included. 

t_methyl <- read_csv("/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_1_Analyst_Role/limma_3_METHYLCHOLANTHRENE_.csv")
#We change the name to padj for convenience
filtered_t_methyl <- t_methyl %>%
  filter(padj < 0.05)

t_limma_written <- write_csv(filtered_t_methyl, "/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_1_Analyst_Role/limma_3_METHYLCHOLANTHRENE_filtered.csv")


t_flucon <- read_csv("/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_1_Analyst_Role/limma_FLUCONAZOLE_.csv")
filtered_t_flucon <- t_flucon %>%
  filter(padj < 0.05)

t_limma_written <- write_csv(filtered_t_flucon, "/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_1_Analyst_Role/limma_flucon_filtered.csv")


t_pirinixic <- read_csv("/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_1_Analyst_Role/limma_PIRINIXIC_ACID_.csv")
filtered_t_pirinixic <- t_pirinixic %>%
  filter(padj < 0.05)

t_limma_written <- write_csv(filtered_t_pirinixic, "/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_1_Analyst_Role/limma_pirinixic_filtered.csv")


#'@details Now that we have the differential expression results from the 
#'Limma process, we can now take the DE results and produce histograms
#'and also volcano plots
#'
#'This process is done three times, once for every single filtered csv file 
#'that deals with its respected chemical 

methyl_info <- read_csv("/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_1_Analyst_Role/limma_3_METHYLCHOLANTHRENE_filtered.csv")

methyl_histogram <- 
  ggplot(data = methyl_info, aes(x = logFC)) +
  geom_histogram(color = '#8EF02B', fill = '#9B82FF') +
  ggtitle('logFC of Significant Genes vs. count') + 
  theme_bw() +
  theme(legend.position = "bottom") 

methyl_histogram

methyl_volcano_plot <- 
  ggplot(data = methyl_info, aes(x = logFC, y = -log(P.Value))) +
  geom_point(color = '#8EF02B', fill = '#9B82FF') +
  ggtitle('logFC vs. -log(P.Value)') + 
  theme_bw() +
  theme(legend.position = "bottom") 

methyl_volcano_plot


flucanazole_info <- read_csv("/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_1_Analyst_Role/limma_flucon_filtered.csv")

flucanazole_histogram <- 
  ggplot(data = flucanazole_info, aes(x = logFC)) +
  geom_histogram(color = '#8EF02B', fill = '#9B82FF') +
  ggtitle('logFC vs. count') + 
  theme_bw() +
  theme(legend.position = "bottom") 

flucanazole_histogram


flucanazole_volcano_plot <- 
  ggplot(data = flucanazole_info, aes(x = logFC, y = -log(P.Value))) +
  geom_point(color = '#8EF02B', fill = '#9B82FF') +
  ggtitle('logFC vs. -log(P.Value)') + 
  theme_bw() +
  theme(legend.position = "bottom")

flucanazole_volcano_plot 


pirinixic_acid_info <- read_csv("/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_1_Analyst_Role/limma_pirinixic_filtered.csv")

pirinixic_acid_histogram <- 
  ggplot(data = pirinixic_acid_info, aes(x = logFC)) +
  geom_histogram(color = '#8EF02B', fill = '#9B82FF') +
  ggtitle('logFC vs. count') + 
  theme_bw() +
  theme(legend.position = "bottom")

pirinixic_acid_histogram


pirinixic_acid_volcano_plot <- 
  ggplot(data = pirinixic_acid_info, aes(x = logFC, y = -log(P.Value))) +
  geom_point(color = '#8EF02B', fill = '#9B82FF') +
  ggtitle('logFC vs. -log(P.Value)') + 
  theme_bw() +
  theme(legend.position = "bottom")

pirinixic_acid_volcano_plot


############ Concordance between microarray and RNA-Seq DE Genes ############ 

map_probe_id_file <- read_csv("/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_1_Analyst_Role/Part_6/refseq_affy_map.csv")

#'@details Next, we want to map the Affymetrix Probe IDs that are from the 
#'microarray analysis section to the refSeq identifies that are used by the 
#'RNA-Seq analysis. It is important that "common" variables are formed between
#'the dataframes so that they can be merged together. 

#Read in all of the DESeq Files that were previously produced
deseq_AhR <- read_csv('/projectnb2/bf528/users/lava-lamp/BF528_Project_3/analysis/DESeqAhRresults.csv')
deseq_car_pxr <- read_csv('/projectnb2/bf528/users/lava-lamp/BF528_Project_3/analysis/DESeqCAR-PXRresults.csv')
deseq_PPARA <- read_csv('/projectnb2/bf528/users/lava-lamp/BF528_Project_3/analysis/DESeqPPARAresults.csv')

#'Filter the DESeq files so only significant genes remain.This means that we 
#'need to take find pvalue < 0.05

deseq_AhR <- deseq_AhR %>%
  filter(pvalue < 0.05)

deseq_car_pxr <- deseq_car_pxr %>%
  filter(pvalue < 0.05)

deseq_PPARA <- deseq_PPARA %>%
  filter(pvalue < 0.05)


#'Next, we need to change the name 'probeid' from the map_probe_id_file so that
#'it becomes 'Gene' because we can combine the results from 
#'limma with the map via "Gene". We can then take that merged dataframe and then 
#'merge with the DESeq results via the "REFSEQ ID"

colnames(map_probe_id_file)[2] <- "Gene"
colnames(deseq_AhR)[1] <- 'REFSEQ'
colnames(deseq_car_pxr)[1] <- 'REFSEQ'
colnames(deseq_PPARA)[1] <- 'REFSEQ'

combination_1 <- merge(methyl_info, map_probe_id_file, by = 'Gene')
combination_1 <- merge(combination_1, deseq_AhR, by = 'REFSEQ')

combination_2 <- merge(flucanazole_info, map_probe_id_file, by = 'Gene')
combination_2 <- merge(combination_2, deseq_car_pxr, by = 'REFSEQ')

combination_3 <- merge(pirinixic_acid_info, map_probe_id_file, by = 'Gene')
combination_3 <- merge(combination_3, deseq_PPARA, by = 'REFSEQ')


#'@details Next, we have to actually perform the calculations in order to find 
#'the concordance values. The wang paper states that the best way to 
#'calculate this is via a formula: 
#'
#'@param n1 This is the first gene set from the microarray data
#'@param n2 This is the second gene set DESeq data
#'@param nx This is the number of 'true' overlapping genes 
#'@param N This is the total number of genes in the genome. In this case 
#'it would be 25,225
#'
#'@formula The formula to actually find the concorodance is: 
#'
#'                 nx = ((N*n0) - (n1*n2)) / (n0+N-n1-n2)
#'
#'@param n0 There is a formula to calculate this
#'
#'
#'                n0 = nx + (((n1-nx)(n2-nx)) / n0 + N - n1 - n2)    

n0_first_set <- nrow(combination_1)
n1_first_set <- nrow(methyl_info)
n2_first_set <- nrow(deseq_AhR)

n0_second_set <- nrow(combination_2)
n1_second_set <- nrow(flucanazole_info)
n2_second_set <- nrow(deseq_car_pxr)

n0_third_set <- nrow(combination_3)
n1_third_set <- nrow(pirinixic_acid_info)
n2_third_set <- nrow(deseq_PPARA)

N = 25225

#'@details Now that we have the variables needed in order to find the 
#'"true" overlapping genes we can go ahead and calculate the nx values for 
#'each set. The nx value represents the number of concordant pairs

nx_first_set <- ((N*n0_first_set) - (n1_first_set*n2_first_set)) / (n0_first_set+N-n1_first_set-n2_first_set)

nx_second_set <- ((N*n0_second_set) - (n1_second_set*n2_second_set)) / (n0_second_set+N-n1_second_set-n2_second_set)

nx_third_set <- ((N*n0_third_set) - (n1_third_set*n2_third_set)) / (n0_third_set+N-n1_third_set-n2_third_set)


#'@details Finally, we can calculate the concordance. According to the Wang paper
#'the formula for the concordance value is 
#'
#'                   concordance = 2(nx) / (n1+n2)  
    

concordance_1 <- (2*nx_first_set) / (n1_first_set - n2_first_set)
concordance_2 <- (2*nx_second_set) / (n1_second_set - n2_second_set)
concordance_3 <- (2*nx_third_set) / (n1_third_set - n2_third_set)
Concordance_Values <- c(concordance_1, concordance_2, concordance_3)

#There are the DEG values for the microarray analysis
Number_of_DEGs <- c(n1_first_set, n1_second_set, n1_third_set)

#These are the DEG values for the RNA-Seq analysis
DEGs_sequence <- c(n2_first_set, n2_second_set, n2_second_set)

concordance_DEGs_df <- data.frame(Concordance_Values, Number_of_DEGs)

concordance_DEGs_sequence_df <- data.frame(Concordance_Values, DEGs_sequence)

#'@details We can now go ahead and plot the concordance values against the 
#'number of DEGs from both the microarray analysis as well as the RNA-Seq 
#'analysis  

concordance_DEGs_plot <- 
  ggplot(concordance_DEGs_df, aes(x = Number_of_DEGs, y = Concordance_Values))+
  geom_point(color = '#8EF02B', fill = '#9B82FF')+
  theme_bw()+
  geom_smooth(method=lm, se=FALSE)+
  labs(xlab = 'Quantity of DEGs from Microarray Analysis')+
  theme(legend.position = "bottom")


concordance_DEGs_plot

concordance_DEGs_sequence_plot <- 
  ggplot(concordance_DEGs_sequence_df, aes(x = DEGs_sequence, y = Concordance_Values))+
  geom_point(color = '#8EF02B', fill = '#9B82FF')+
  theme_bw()+
  geom_smooth(method=lm, se=FALSE)+
  theme(legend.position = "bottom")

concordance_DEGs_sequence_plot

#'@details Next we need to split the differential expressed genes from both 
#'methodologies, DESeq and Limma, into two different groups: 'above-median' 
#'and also 'below-median'.  
#'
#'We need to find the median of both the AveExpr row in the limma results 
#'and the baseMean limma results 

#Here we can calculate the median via the median() function
median_deseq_AhR <- median(deseq_AhR$baseMean)

#Here we want all the baseMean values that are greater than the median value 
#of the baseMean column 
above_median_deseq_AhR <- deseq_AhR %>%
  filter(baseMean > median_deseq_AhR )

#Here we want all the baseMean values that are less than the median value 
#of the baseMean column 
below_median_deseq_AhR <- deseq_AhR %>%
  filter(baseMean < median_deseq_AhR )

#Median value, above, and below filtering for 3_METHYLCHOLANTHRENE
median_limma_methyl <- median(methyl_info$AveExpr)

above_median_limma_methyl <- methyl_info %>%
  filter(AveExpr > median_limma_methyl)

below_median_limma_methyl <- methyl_info %>%
  filter(AveExpr < median_limma_methyl)

#Median value, above, and below filtering for DESeq genes in deseq_car_pxr
median_deseq_car_pxr <- median(deseq_car_pxr$baseMean)

above_median_deseq_car_pxr <- deseq_car_pxr %>%
  filter(baseMean > median_deseq_car_pxr)

below_deseq_car_pxr <- deseq_car_pxr %>%
  filter(baseMean < median_deseq_car_pxr)

#Median value, above, and below filtering for limma genes in flucanazole
median_limma_flucanazole_info <- median(flucanazole_info$AveExpr)

above_median_limma_flucanazole_info <- flucanazole_info %>%
  filter(AveExpr > median_limma_flucanazole_info)

below_median_limma_flucanazole_info <- flucanazole_info %>%
  filter(AveExpr < median_limma_flucanazole_info)


#Median value, above, and below filtering for DESeq genes in deseq_PPARA
median_deseq_PPARA <- median(deseq_PPARA$baseMean)

above_median_deseq_PPARA <- deseq_PPARA %>%
  filter(baseMean > median_deseq_PPARA)

below_median_deseq_PPARA <- deseq_PPARA %>%
  filter(baseMean < median_deseq_PPARA)


#Median value, above, and below filtering for DESeq genes in deseq_PPARA
median_limma_pirinixic_acid_info <- median(pirinixic_acid_info$AveExpr)

above_median_limma_pirinixic_acid_info <- pirinixic_acid_info %>%
  filter(AveExpr > median_limma_pirinixic_acid_info) 

below_median_limma_pirinixic_acid_info <- pirinixic_acid_info %>%
  filter(AveExpr < median_limma_pirinixic_acid_info)

#'@details Now that we have split all of the DE genes into above-median 
#'and below-median groups, we now have to compute the concordance of each 
#'group of genes. All of the calculations below will be for the 
#'above-genes 
#'

combo_above_1 <- merge(above_median_limma_methyl, map_probe_id_file, by = "Gene")
combo_above_1 <- merge(combo_above_1, above_median_deseq_AhR, by = 'REFSEQ')

combo_above_2 <- merge(above_median_limma_flucanazole_info, map_probe_id_file, by = 'Gene')
combo_above_2 <- merge(combo_above_2, above_median_deseq_car_pxr, by = 'REFSEQ')

combo_above_3 <- merge(above_median_limma_pirinixic_acid_info, map_probe_id_file, by = "Gene")
combo_above_3 <- merge(combo_above_3, above_median_deseq_PPARA, by = 'REFSEQ')

n0_first_set_above <- nrow(combo_above_1)
n1_first_set_above <- nrow(above_median_limma_methyl)
n2_first_set_above <- nrow(above_median_deseq_AhR)

n0_second_set_above <- nrow(combo_above_2)
n1_second_set_above <- nrow(above_median_limma_flucanazole_info)
n2_second_set_above <- nrow(above_median_deseq_car_pxr)

n0_third_set_above <- nrow(combo_above_3)
n1_third_set_above <- nrow(above_median_limma_pirinixic_acid_info)
n2_third_set_above <- nrow(above_median_deseq_PPARA)

nx_first_set_above <- ((N*n0_first_set_above) - (n1_first_set_above*n2_first_set_above)) / (n0_first_set_above+N-n1_first_set_above-n2_first_set_above)
nx_second_set_above <- ((N*n0_second_set_above) - (n1_second_set_above*n2_second_set_above)) / (n0_second_set_above+N-n1_second_set_above-n2_second_set_above)
nx_third_set_above <- ((N*n0_third_set_above) - (n1_third_set_above*n2_third_set_above)) / (n0_third_set_above+N-n1_third_set_above-n2_third_set_above)

concordance_1_above <- (2*nx_first_set_above) / (n1_first_set_above - n2_first_set_above)
concordance_2_above <- (2*nx_second_set_above) / (n1_second_set_above - n2_second_set_above)
concordance_3_above <- (2*nx_third_set_above) / (n1_third_set_above - n2_third_set_above)


#'@details Now that we have calculate the concordance values for the above 
#'genes we now have to compute the concordance for the below genes. The 
#'calculations below will deal with the below median genes. 

combo_below_1 <- merge(below_median_limma_methyl, map_probe_id_file, by = "Gene")
combo_below_1 <- merge(combo_below_1, below_median_deseq_AhR, by = 'REFSEQ')

combo_below_2 <- merge(below_median_limma_flucanazole_info, map_probe_id_file, by = 'Gene')
combo_below_2 <- merge(combo_below_2, below_deseq_car_pxr, by = 'REFSEQ')

combo_below_3 <- merge(below_median_limma_pirinixic_acid_info, map_probe_id_file, by = "Gene")
combo_below_3 <- merge(combo_below_3, below_median_deseq_PPARA, by = 'REFSEQ')

n0_first_set_below <- nrow(combo_below_1)
n1_first_set_below <- nrow(below_median_limma_methyl)
n2_first_set_below <- nrow(below_median_deseq_AhR)

n0_second_set_below <- nrow(combo_below_2)
n1_second_set_below <- nrow(below_median_limma_flucanazole_info)
n2_second_set_below <- nrow(below_deseq_car_pxr)

n0_third_set_below <- nrow(combo_below_3)
n1_third_set_below <- nrow(below_median_limma_pirinixic_acid_info)
n2_third_set_below <- nrow(below_median_deseq_PPARA)


nx_first_set_below <- ((N*n0_first_set_below) - (n1_first_set_below*n2_first_set_below)) / (n0_first_set_below+N-n1_first_set_below-n2_first_set_below)
nx_second_set_below <- ((N*n0_second_set_below) - (n1_second_set_below*n2_second_set_below)) / (n0_second_set_below+N-n1_second_set_below-n2_second_set_below)
nx_third_set_below <- ((N*n0_third_set_below) - (n1_third_set_below*n2_third_set_below)) / (n0_third_set_below+N-n1_third_set_below-n2_third_set_below)


concordance_1_below <- (2*nx_first_set_below) / (n1_first_set_below - n2_first_set_below)
concordance_2_below <- (2*nx_second_set_below) / (n1_second_set_below - n2_second_set_below)
concordance_3_below <- (2*nx_third_set_below) / (n1_third_set_below - n2_third_set_below)


#'@details Now that we have all of the concordances for total the DE genes, 
#'those that fall below the median value, and those that are above the median 
#'value, we can produce a combined bar plot with all this information! First, 
#'we need to combine the data that we have produced 

combined_3_methyl <- c(concordance_1, concordance_1_above, concordance_1_below)
combined_fluconazole <- c(concordance_2, concordance_2_above, concordance_2_below)
combined_pirinixic_acid <- c(concordance_3, concordance_3_above, concordance_3_below)

combined_bar_plot <- as.data.frame(rbind(combined_3_methyl, combined_fluconazole, combined_pirinixic_acid))
MOA <- c('AhR/3-methylchol', 'car_pxr/fluconazole', 'PPARA/pirinixic_acid')
row.names(combined_bar_plot) <- MOA
concordance_sections <- c('All DE Genes', "Above Median", "Below Median")
colnames(combined_bar_plot) <- c('All DE Genes', "Above Median", "Below Median")

combined_bar_plot <- tibble::rownames_to_column(combined_bar_plot, 'row_names')

combined_bar_plot_melt <- melt(combined_bar_plot, 'row_names')

concordance_bar_plot <- 
  ggplot(combined_bar_plot_melt, aes(x = row_names, y = value)) +
  geom_bar(aes(fill = variable), stat = 'identity') +
  theme(legend.position = "bottom") +
  theme_bw()

concordance_bar_plot


