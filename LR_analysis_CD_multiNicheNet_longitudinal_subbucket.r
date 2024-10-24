###we are going to run MultiNicheNet analysis
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(multinichenetr)
library(Seurat)
library(tidyverse)

setwd("/gpfs3/well/dendrou/shared/loom/TAURUS/TAURUS_analyses/multinichenetr/")
print("reading in Seurat object")
#read in Seurat object
data <- readRDS('/gpfs3/well/dendrou/shared/loom/TAURUS/TAURUS_analyses/UCell/adata_raw_annotated_Seurat.rds')

##add in disease duration, additional response data, inflammation score factors etc.
diseaseduration <- read_csv('/gpfs3/well/dendrou/shared/loom/TAURUS/TAURUS_analyses/multinichenetr/Disease_duration.csv')
response <- read_csv('/gpfs3/well/dendrou/shared/loom/TAURUS/TAURUS_analyses/multinichenetr/Additional_response.csv')
data@meta.data$Rownames <- rownames(data@meta.data)
data@meta.data <- dplyr::left_join(data@meta.data, diseaseduration, by = "Patient")
data@meta.data <- dplyr::left_join(data@meta.data, response, by = "Patient")
data@meta.data <- data@meta.data %>% tibble::column_to_rownames(var = "Rownames")

data@meta.data  <- data@meta.data  %>% mutate(Group =case_when(MM_scaled <= 6.5 ~ "Non_inflamed", MM_scaled > 6.5 ~ "Inflamed"))
data@meta.data$Group <- factor(data@meta.data$Group, levels = c("Non_inflamed", "Inflamed"))

##filter object so its comparison appropriate
keep_samples <- read_csv('/gpfs3/well/dendrou/shared/loom/TAURUS/TAURUS_analyses/multinichenetr/Longitudinal_CD/Longitudinal_CD_metadata.csv')
data <- subset(data, subset = SampleID2 %in% keep_samples$SampleID2)

###lets also ensure SampleID is patient_site_treatment
data$Patient_Site_Treatment <- paste0(data$Patient,"_",data$Site,"_",data$Treatment)
data$Remission_Treatment <- paste0(data$Remission,"_",data$Treatment)

#we need to homogenise the minor populations now
data@meta.data$minor <- as.character(data@meta.data$minor)
data@meta.data$minor[data@meta.data$minor == "Ileal_BEST4_OTOP2"] <- "BEST4_OTOP2"
data@meta.data$minor[data@meta.data$minor == "Non_ileal_BEST4_OTOP2"] <- "BEST4_OTOP2"

data@meta.data$minor[data@meta.data$minor == "Ileal_EEC"] <- "EEC"
data@meta.data$minor[data@meta.data$minor == "Non_ileal_EEC"] <- "EEC"

data@meta.data$minor[data@meta.data$minor == "Ileal_enterocyte"] <- "Enterocyte"
data@meta.data$minor[data@meta.data$minor == "Non_ileal_enterocyte"] <- "Enterocyte"

data@meta.data$minor[data@meta.data$minor == "Ileal_goblet"] <- "Goblet"
data@meta.data$minor[data@meta.data$minor == "Non_ileal_goblet"] <- "Goblet"

data@meta.data$minor[data@meta.data$minor == "Ileal_paneth"] <- "Paneth"
data@meta.data$minor[data@meta.data$minor == "Non_ileal_paneth"] <- "Paneth"

data@meta.data$minor[data@meta.data$minor == "Ileal_stem_LGR5pos"] <- "Stem_LGR5pos"
data@meta.data$minor[data@meta.data$minor == "Non_ileal_stem_LGR5pos"] <- "Stem_LGR5pos"

data@meta.data$minor[data@meta.data$minor == "Ileal_tuft"] <- "Tuft"
data@meta.data$minor[data@meta.data$minor == "Non_ileal_tuft"] <- "Tuft"

data@meta.data$minor[data@meta.data$minor == "Non_ileal_M_like"] <- "M_like"
data@meta.data$minor <- as.factor(data@meta.data$minor)

#Alternatively, we can just re-label epithelium if moving at sub_bucket level
data@meta.data$sub_bucket <- as.character(data@meta.data$sub_bucket)
data@meta.data$sub_bucket[data@meta.data$sub_bucket == "Ileal_epithelium"] <- "Epithelium"
data@meta.data$sub_bucket[data@meta.data$sub_bucket == "Non_ileal_epithelium"] <- "Epithelium"
data@meta.data$sub_bucket <- as.factor(data@meta.data$sub_bucket)

#convert to SCE object
data = Seurat::as.SingleCellExperiment(data, assay = "RNA")

###QC checking
cells_by_sample <- as.data.frame(table(SummarizedExperiment::colData(data)$sub_bucket, SummarizedExperiment::colData(data)$sub_bucket))
cells_by_condition <- as.data.frame(table(SummarizedExperiment::colData(data)$sub_bucket, SummarizedExperiment::colData(data)$Remission_Treatment))

##set up MultiNicheNet  
lr_network = readRDS("lr_network.RDS")
lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor) %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor))
ligand_target_matrix = readRDS("ligand_target_matrix.RDS")
colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()
data = alias_to_symbol_SCE(data, "human") %>% makenames_SCE()

#Step 0.3: Prepare settings of the MultiNicheNet cell-cell communication analysis
SummarizedExperiment::colData(data)$Patient_Site_Treatment = SummarizedExperiment::colData(data)$Patient_Site_Treatment %>% make.names()
SummarizedExperiment::colData(data)$sample_id = SummarizedExperiment::colData(data)$sample_id %>% make.names()
SummarizedExperiment::colData(data)$Disease = SummarizedExperiment::colData(data)$Disease %>% make.names()
SummarizedExperiment::colData(data)$Remission_Treatment = SummarizedExperiment::colData(data)$Remission_Treatment %>% make.names()

#Step 0.3: Prepare settings of the MultiNicheNet cell-cell communication analysis
sample_id = "Patient_Site_Treatment"
group_id = "Remission_Treatment"
celltype_id = "sub_bucket"
covariates = NA 
batches = NA

#All vs all analysis
senders_oi = SummarizedExperiment::colData(data)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(data)[,celltype_id] %>% unique()

print("Running abundance")
#Step 1: Extract cell type abundance and expression information from receiver and sender cell types, and link this expression information for ligands of the sender cell types to the corresponding receptors of the receiver cell types
min_cells = 10

###get abundance info
abundance_expression_info = get_abundance_expression_info(sce = data, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, min_cells = min_cells, senders_oi = senders_oi, receivers_oi = receivers_oi, lr_network = lr_network, batches = batches)

###Step 2: Perform genome-wide differential expression analysis of receiver and sender cell types to define DE genes between the conditions of interest. Based on this analysis, we can define the logFC/p-value of ligands in senders and receptors in receivers, and define the set of affected target genes in the receiver.
#Define the contrasts and covariates of interest for the DE analysis.
contrasts_oi = c("'(Responder_Post-Responder_Pre)-(Non_responder_Post-Non_responder_Pre)','(Non_responder_Post-Non_responder_Pre)-(Responder_Post-Responder_Pre)'")
contrast_tbl = tibble(contrast =
                        c("(Responder_Post-Responder_Pre)-(Non_responder_Post-Non_responder_Pre)", "(Non_responder_Post-Non_responder_Pre)-(Responder_Post-Responder_Pre)"),
                      group = c("Responder_Post","Non_responder_Post")) 

 


print("Running DEG")
#Perform the DE analysis for each cell type.
DE_info = get_DE_info(sce = data, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, batches = batches, covariates = covariates, contrasts_oi = contrasts_oi, min_cells = min_cells)

empirical_pval = FALSE
if(empirical_pval == TRUE){
  DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
} 

empirical_pval = FALSE
if(empirical_pval == FALSE){
  celltype_de = DE_info$celltype_de$de_output_tidy
} else {
  celltype_de = DE_info_emp$de_output_tidy_emp %>% dplyr::select(-p_val, -p_adj) %>% dplyr::rename(p_val = p_emp, p_adj = p_adj_emp)
}

print("Combining DEG info")
#Combine DE information for ligand-senders and receptors-receivers (similar to step1 - abundance_expression_info$sender_receiver_info)
sender_receiver_de = combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)

#Step 3: Predict NicheNet ligand activities and NicheNet ligand-target links based on these differential expression results
#Define the parameters for the NicheNet ligand activity analysis
logFC_threshold = 0.50
p_val_threshold = 0.05
fraction_cutoff = 0.05

#p_val_adj = TRUE - if false it counts the threshold towards unadjusted p values
p_val_adj = TRUE

#select which top n of the predicted target genes will be considered (here: top 250 targets per ligand).
top_n_target = 250

verbose = TRUE
cores_system = 1
n.cores = min(cores_system, sender_receiver_de$receiver %>% unique() %>% length()) # use one core per receiver cell type

print("Running NicheNet")
#Run the NicheNet ligand activity analysis
ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(get_ligand_activities_targets_DEgenes(
  receiver_de = celltype_de,
  receivers_oi = receivers_oi,
  ligand_target_matrix = ligand_target_matrix,
  logFC_threshold = logFC_threshold,
  p_val_threshold = p_val_threshold,
  p_val_adj = p_val_adj,
  top_n_target = top_n_target,
  verbose = verbose, 
  n.cores = n.cores
)))

#Step 4: Use the information collected above to prioritize all sender-ligand—receiver-receptor pairs.
prioritizing_weights_DE = c("de_ligand" = 1,
                         "de_receptor" = 1)
prioritizing_weights_activity = c("activity_scaled" = 2)

prioritizing_weights_expression_specificity = c("exprs_ligand" = 2,
                         "exprs_receptor" = 2)

prioritizing_weights_expression_sufficiency = c("frac_exprs_ligand_receptor" = 1)

prioritizing_weights_relative_abundance = c( "abund_sender" = 0,
                         "abund_receiver" = 0)

prioritizing_weights = c(prioritizing_weights_DE, 
                         prioritizing_weights_activity, 
                         prioritizing_weights_expression_specificity,
                         prioritizing_weights_expression_sufficiency, 
                         prioritizing_weights_relative_abundance)

sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(data) %>% tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group")
}

path = "./Longitudinal_CD/sub_bucket/"

list(abundance_expression_info = abundance_expression_info, grouping_tbl = grouping_tbl, sender_receiver_de = sender_receiver_de, ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes, contrast_tbl = contrast_tbl, sender_receiver_tbl = sender_receiver_tbl, celltype_de = celltype_de) %>% saveRDS(paste0(path,"intermediary_output_pre_prioritisation.rds"))


print("Running prioritisation")
##Run the prioritization
prioritization_tables = suppressMessages(generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  prioritizing_weights = prioritizing_weights,
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender
))

print("Run correlation of LR and target expression")
#Step 5: Add information on prior knowledge and expression correlation between LR and target expression.
lr_target_prior_cor = lr_target_prior_cor_inference(prioritization_tables$group_prioritization_tbl$receiver %>% unique(), abundance_expression_info, celltype_de, grouping_tbl, prioritization_tables, ligand_target_matrix, logFC_threshold = logFC_threshold, p_val_threshold = p_val_threshold, p_val_adj = p_val_adj)

print("Saving multiNicheNet results")

multinichenet_output = list(
    celltype_info = abundance_expression_info$celltype_info,
    celltype_de = celltype_de,
    sender_receiver_info = abundance_expression_info$sender_receiver_info,
    sender_receiver_de =  sender_receiver_de,
    ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
    prioritization_tables = prioritization_tables,
    grouping_tbl = grouping_tbl,
    lr_target_prior_cor = lr_target_prior_cor
  ) 
multinichenet_output = make_lite_output(multinichenet_output)

save = TRUE
if(save == TRUE){
  saveRDS(multinichenet_output, paste0(path, "multinichenet_output_longitudinal_CD.rds"))

}

print("Completed")
