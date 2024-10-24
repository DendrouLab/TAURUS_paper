# Libraries
library(MAST)
library(SingleCellExperiment)
library(ggpubr)
library(qs)
library(optparse)
library(data.table)
library(dittoSeq)
library(EnhancedVolcano)


# Options
option_list = list(
  make_option(c("--filename"), type="character", default="/well/cartography/projects/analysis/302_external_oxford_datasets/TAURUS/data/subsetted_colon_final_analysis_celltype/colon_ABCA8pos_WNT2Bpos_FOShi_fibroblast_sce.qs"),
  make_option(c("--celltype"), type="character", default="ABCA8pos_WNT2Bpos_FOShi_fibroblast"),
  make_option(c("--out_dir"), type="character", default="/well/cartography/projects/analysis/302_external_oxford_datasets/TAURUS/revisions/results/MAST/CD_colon_samples/final_analysis")
  )
#
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


message(paste0("Outdir: ", opt$out_dir))
message("Runnin Treatment analysis")

# Parameters
freq_expressed <- 0.1
FCTHRESHOLD <- log2(1.2)
fcthreshold <- 0.25
message("LOG FC Threshold: 0.25")

# Import data
sce <- qread(opt$filename)
sce 
# Sanity check
message("Making Sure we are running CD only and which Remisison status ")
print(table(sce$Disease))
print(table(sce$Remission))

# # MAST
celltype_ofinterest <- as.character(opt$celltype)
print(celltype_ofinterest)

message("Starting MAST")

sca <- SceToSingleCellAssay(sce)
rm(sce)
gc()

## Prelim CDR
colData(sca)$nGeneOn <- colSums(assay(sca)>0)

## Remove invariant genes
sca <- sca[freq(sca)>0,]

## Recalculate CDR
cdr2 <-colSums(assay(sca)>0)
colData(sca)$cngeneson <- scale(cdr2)
message("Done re-calculating CDR")


expressed_genes <- freq(sca) > freq_expressed
message("# of Expressed genes: ", length(which(expressed_genes==TRUE)))

sca <- sca[expressed_genes,]

# Save Results
savefile<-file.path(opt$out_dir, paste0("MAST", celltype_ofinterest, ".Rdata"))

Treatment <-factor(colData(sca)$Treatment)

Treatment <-relevel(Treatment ,"Pre")

colData(sca)$Treatment <- Treatment

Remission <-factor(colData(sca)$Remission)

Remission <-relevel(Remission ,"Responder")

colData(sca)$Remission <- Remission

colData(sca)$Site <-factor(colData(sca)$Site)

message("Starting zlm")
# A note on the nAGQ=0 argument for fitting the discrete component of MAST's hurdle model:
# 
# Fitting the model involves optimizing an objective function, the Laplace approximation to the deviance, with respect to the
# parameters. The Laplace approximation to the deviance requires determining the conditional modes of the random effects.
# These are the values that maximize the conditional density of the random effects, given the model parameters and the data.
# This is done using Penalized Iteratively Reweighted Least Squares (PIRLS). In most cases PIRLS is fast and stable. It is simply
# a penalized version of the IRLS algorithm used in fitting GLMs.
#
# The distinction between the "fast" (nAGQ=0) and "slow" (nAGQ=1) algorithms in lme4 (to which glmer belongs)
# is whether the fixed-effects parameters are optimized in PIRLS or a nonlinear optimizer.

#zlmCond <- zlm(~Disease + cngeneson + (1|Patient), sca,method = 'glmer', ebayes = FALSE)
# dint include sample_id as semed like over kill

message("Formula:~Remission+ Treatment + Age + Gender + Site+ Disease_duration + cngeneson + Remission:Treatment + (1|Patient/sample_id)")

zlmCond <- zlm(~Remission+ Treatment + Age + Gender + Site + Disease_duration + cngeneson + Remission:Treatment + (1|Patient/sample_id), sca, method = 'glmer', ebayes = FALSE, fitArgsD = list(nAGQ = 0))

message("Done zlm")

#Run likelihood ratio test for the condition coefficient.
summary_1 <- summary(zlmCond, doLRT = "RemissionNon_responder:TreatmentPost")


message("Saving zlm")

save(sca,zlmCond,summary_1,file=savefile)

# Plot inflammation_final
message("Plotting Treatment")

summaryDt <- summary_1$datatable
fcHurdle <- merge(summaryDt[contrast=="RemissionNon_responder:TreatmentPost" & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryDt[contrast=="RemissionNon_responder:TreatmentPost" & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleSig<- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig, fdr)
setorder(fcHurdle, fdr)

message("Saving fcHurdleSig and fchurdle for treatment")
write.csv(fcHurdle, file.path(opt$out_dir, paste0('MAST_fcHurdle_summary_remissionNR_treatmentpost_with_site', celltype_ofinterest,'.csv')),row.names=FALSE)
write.csv(fcHurdleSig, file.path(opt$out_dir, paste0('MAST_fcHurdleSig_summary_remissionNR_treatmentpost_with_site', celltype_ofinterest,'.csv')),row.names=FALSE)

## Volcano
rownames(fcHurdle) <- fcHurdle$primerid
plot_volcano <- EnhancedVolcano(fcHurdle,
                lab = rownames(fcHurdle),
                x = 'coef',
                y = 'fdr', 
                FCcutoff = fcthreshold, 
                pCutoff = 5e-02, 
                title = paste0("fcHurdle Remission*Treatment: ", celltype_ofinterest), 
                subtitle = "Non_responder(Post-Pre) (right), Responder(Post-Pre) (left)")
ggsave(plot = plot_volcano, filename = file.path(opt$out_dir,"figures", paste0("EnhancedVolcano_remission_NR_treatment_post_wih_site", celltype_ofinterest,"_FC", fcthreshold, ".png")))

