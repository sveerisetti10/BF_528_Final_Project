#Author: Sriramteja Veerisetti 
#BF 528 
#Analyst Role in Project 3
#Microarray Differential Expression with Limma 
#Concordance between microarray and RNA-Seq DE Genes 
#Individual Project 

#'Note from Sri: Make sure to run install the limma package if you do not have 
#'on your computer! 
#BiocManager::install("limma")
library(tidyverse)
library(dplyr)
library(limma)


#'@details We can utilize the next few steps for each chemical we decide to 
#'utilize. For example, we will only be using the chemical: 3-METHYLCHOLANTHRENE
#'for the next few steps. 
#'
#'

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
#'
#'
#'We will start by designing a model matrix so that we can utilize 
#'it in the lmfit() function 

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
#'
#'
#'We will start by designing a model matrix so that we can utilize 
#'it in the lmfit() function 

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
#'
#'This process is done three times, once for every single filtered csv file 
#'that deals with its respected chemical 

methyl_info <- read_csv("/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_1_Analyst_Role/limma_3_METHYLCHOLANTHRENE_filtered.csv")
methyl_histogram <- 
  ggplot(data = methyl_info, aes(x = logFC)) +
  geom_histogram(color = '#8EF02B', fill = '#9B82FF') +
  ggtitle('logFC vs. count') + 
  theme_bw() +
  theme(legend.position = "bottom") 
  
  
flucanazole_info <- read_csv("/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_1_Analyst_Role/limma_flucon_filtered.csv")
flucanazole_histogram <- 
  ggplot(data = flucanazole_info, aes(x = logFC)) +
  geom_histogram(color = 'red', fill = 'blue') +
  ggtitle('logFC vs. count') + 
  theme_bw() +
  theme(legend.position = "bottom") 


pirinixic_acid_info <- read_csv("/usr4/bf528/sv/Documents/BF_528_Individual_Project/Project_1_Analyst_Role/limma_pirinixic_filtered.csv")
pirinixic_acid_histogram <- 
  ggplot(data = pirinixic_acid_info, aes(x = logFC)) +
  geom_histogram(color = 'blue', fill = 'red') +
  ggtitle('logFC vs. count') + 
  theme_bw() +
  theme(legend.position = "bottom") 

pirinixic_acid_histogram
