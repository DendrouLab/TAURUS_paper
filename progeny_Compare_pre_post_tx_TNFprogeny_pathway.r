library(ggplot2)
library(tidyverse)
library(lme4)

data <- read.csv('progeny_processed.csv')

#subset_data <- read.csv('/Users/tomthomas30/Desktop/Rebuttal/Abundance/Treatment_Remission/Longitudinal/R_CD_Remission_metadata.csv')
#subset_data <- read.csv('/Users/tomthomas30/Desktop/Rebuttal/Abundance/Treatment_Remission/Longitudinal/NR_CD_Remission_metadata.csv')
subset_data <- read.csv('/Users/tomthomas30/Desktop/Rebuttal/Abundance/Treatment_Remission/Longitudinal/R_UC_Remission_metadata.csv')
#subset_data <- read.csv('/Users/tomthomas30/Desktop/Rebuttal/Abundance/Treatment_Remission/Longitudinal/NR_UC_Remission_metadata.csv')


###get 75th percentile by group - minor
data <- data %>% filter(sample_id %in% subset_data$SampleID3) #subset samples
data_summary <- data %>% group_by(minor, sample_id) %>% summarize(TNFa_score = quantile(TNFa, c(0.75)))
###make metadata
###make metadata - FOR UC PAIRED, SUBSTITUTE THE X FOR BARCODE
data <- subset(data, select = -c(barcode, final_analysis, minor, major, sub_bucket, bucket,Androgen,EGFR,Estrogen,Hypoxia,JAK.STAT,MAPK, NFkB, PI3K, TGFb, TNFa, Trail, VEGF, WNT, p53))
data <- data[!duplicated(data[ , c("sample_id")]), ]
#merge response status into score
response <- read.csv('/Users/tomthomas30/Desktop/Rebuttal/Longitudinal_GEX/Additional_response.csv')
data <- merge(data, response, by = "Patient", all.x = TRUE)
###merge metadata into score
data_summary  <- merge(data_summary, data, by = "sample_id", all.x = TRUE)
df = list()
df_result = list()

to_loop <- unique(data_summary$minor)

for (i in to_loop){
    data_mock <- data_summary %>% filter(minor == i)
    ###plot the figure
    p <- ggplot(data=data_mock, aes_string(x = 'Treatment', y = 'TNFa_score', fill='Treatment')) + geom_dotplot(binaxis='y', stackdir='center')
    p <- p + stat_summary(fun.y=median, geom="point", shape=18, size=3, color="red")+ theme_classic()
    p <- p + ggtitle(i) + theme(plot.title = element_text(hjust = 0.5))
    df[[i]] <- p
    ###save the data
    result <- lmerTest::lmer(TNFa_score ~ Treatment + (1 | Patient), data = data_mock, control = lmerControl(optimizer = "bobyqa"))
    result <- as.data.frame(summary(result)$coefficients)
    result <- result[-1,]
    rownames(result) = paste0(i,"_",rownames(result))
    df_result[[i]] <- result
}

pdf("processed/pre_post_tx/UCresp_minor.pdf")
df
dev.off()
df_total <- bind_rows(df_result, .id = "column_label")
df_total$FDR <- p.adjust(df_total$`Pr(>|t|)`, method = "fdr", n = length(df_total$`Pr(>|t|)`))
write.csv(df_total, 'processed/pre_post_tx/UCresp_minor.csv')
