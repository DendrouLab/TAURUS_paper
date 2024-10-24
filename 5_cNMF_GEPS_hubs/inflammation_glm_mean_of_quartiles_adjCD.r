library(dplyr)
library(tidyverse)
library(ggplot2)
library(viridis)
library(stringr)
library(optparse)
library(rstatix)
library(reshape2)
library(cowplot)
library(tidystats)
library(purrr)
library(data.table)
library(HardyWeinberg)
library(igraph)
library(pheatmap)
library(ComplexHeatmap)
library(stringi)
library(ggpubr)
library(matrixStats)
library(lme4)
library(tibble)

program_hub <- read.csv('q0.75_partitions_labelled.txt', sep = '\t')

###load in programs - with expression 
filenames <- list.files(recursive = FALSE, pattern = "_program_hub_discovery.csv")
df_list <- lapply(filenames, read.csv, header = TRUE)
final <- bind_rows(df_list)
final <- final %>% dplyr::select(sample_id, program, q0.99, q0.95, q0.75, q0.5, q0.25)
final <- final %>% rowwise() %>% dplyr::mutate(Mean_score = mean(c(q0.99,q0.95,q0.75,q0.5,q0.25)))
final$program <- stri_replace_last_fixed(final$program, '_', '_cNMF_Usage_')

###read in metadata
metadata <- read.csv('metadata.csv')

###merge into the metadata - for CD
metadata <- merge(x=metadata,y=final,by="sample_id",all.x=TRUE)
metadata <- metadata %>% filter(Disease == "CD")
metadata$Inflammation_revised <- "Inflamed"
metadata$Inflammation_revised[metadata$MM_scaled < 6.5] <- "Non_inflamed"
metadata <- merge(x=metadata,y=program_hub,by="program",all.x=TRUE)

disease_duration <- read.csv('../Data/Disease_duration.csv')
metadata <- merge(metadata, disease_duration, by = "Patient", all.x = TRUE)

##some soupy programs will have NA in bucket column so lets remove this
metadata <- metadata[!is.na(metadata$bucket),]

###loop through the bucket, get dotplot, get table
df_fdr_list = list()
df_nofdr_list = list()

for(i in unique(metadata$bucket)){
    y = toString(i)
    bucket_isolated <- metadata %>% filter(bucket == i)
    ###print it out individually
    pdf(paste0('Sensitivity_analysis_adj_CD/mean_of_quartiles/',y,'_boxplots_inflammation.pdf'), width=15, height=6)
    print(ggboxplot(bucket_isolated, "program", "Mean_score", color = "Inflammation_revised", palette = c("#00AFBB", "#E7B800")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
    dev.off()

    test_list = list()
    for (j in unique(bucket_isolated$program)){
        z = toString(j)
        GEP_isolated <- bucket_isolated %>% filter(program == j)
        ###do the actual GLM model here
        result <- lmerTest::lmer(Mean_score ~ Inflammation_revised + Treatment + Age + Gender + Site + Disease_duration + (1 | Patient), data = GEP_isolated, control = lmerControl(optimizer = "bobyqa"))
        result <- as.data.frame(summary(result)$coefficients)
        result <- tibble::rownames_to_column(result, "rownames")
        result$program <- z
        result$bucket <- y
        test_list[[j]] <- result
    }

    result_bucket <- bind_rows(test_list, .id = "column_label")
    df_nofdr_list[[i]] <- result_bucket

    ###collate the inflammation specific pvalues and apply fdr
    result_bucket <- result_bucket %>% filter(rownames == "Inflammation_revisedNon_inflamed")
    result_bucket$FDR <- p.adjust(result_bucket$`Pr(>|t|)`, method = "fdr")
    df_fdr_list[[i]] <- result_bucket
}

result_across_buckets <- bind_rows(df_fdr_list, .id = "column_label")
write.csv(result_across_buckets, 'Sensitivity_analysis_adj_CD/mean_of_quartiles/inflammation_tx_adjusted_CD.csv', row.names=FALSE)
