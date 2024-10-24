library(tidyverse)
library(MASC)
library(lme4)
library(ggrepel)

## Adapted from the MASC R package (Fonseka, Rao, et al. Sci Trans Med 2017)
MASC.me <- function(dataset, cluster, contrast, random_effects = NULL, fixed_effects = NULL,
                 verbose = FALSE, save_models = FALSE, save_model_dir = NULL, weights = NULL) {


  # Generate design matrix from cluster assignments
  cluster <- as.character(cluster)
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, dataset)

  # Convert cluster assignments to string
  cluster <- as.character(cluster)
  # Prepend design matrix generated from cluster assignments
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, dataset)
  # Create output list to hold results
  res <- vector(mode = "list", length = length(unique(cluster)))
  names(res) <- attributes(designmat)$dimnames[[2]]

  # Create model formulas
  if (!is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0(c(paste0(fixed_effects, collapse = " + "),
                          paste0("(1|", random_effects, ")", collapse = " + ")),
                        collapse = " + ")
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
    }
  } else if (!is.null(fixed_effects) && is.null(random_effects)) {
    model_rhs <- paste0(fixed_effects, collapse = " + ")
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
      # For now, do not allow models without mixed effects terms
      stop("No random effects specified")
    }
  } else if (is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0("(1|", random_effects, ")", collapse = " + ")
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
    }
  } else {
    model_rhs <- "1" # only includes intercept
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
      stop("No random or fixed effects specified")
    }
  }

  # Initialize list to store model objects for each cluster
  cluster_models <- vector(mode = "list",
                           length = length(attributes(designmat)$dimnames[[2]]))
  names(cluster_models) <- attributes(designmat)$dimnames[[2]]

  # Run nested mixed-effects models for each cluster
  for (i in seq_along(attributes(designmat)$dimnames[[2]])) {
    test_cluster <- attributes(designmat)$dimnames[[2]][i]
    if (verbose == TRUE) {
      message(paste("Creating logistic mixed models for", test_cluster))
    }
    null_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ 1 + "),
                                   model_rhs), collapse = ""))
    full_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ ", contrast, " + "),
                                   model_rhs), collapse = ""))
    # Run null and full mixed-effects models
    null_model <- lme4::glmer(formula = null_fm, data = dataset,
                              family = binomial, nAGQ = 1, verbose = 0,
                              control = glmerControl(optimizer = "bobyqa"), weights = weights)
 #   print(summary(null_model))
    full_model <- lme4::glmer(formula = full_fm, data = dataset,
                              family = binomial, nAGQ = 1, verbose = 0,
                              control = glmerControl(optimizer = "bobyqa"), weights = weights)
 #   print(summary(full_model))
    flush.console()
    model_lrt <- anova(null_model, full_model)
    # calculate confidence intervals for contrast term beta
    contrast_lvl2 <- paste0(contrast, levels(dataset[[contrast]])[2])
    contrast_ci <- confint.merMod(full_model, method = "Wald",
                                  parm = contrast_lvl2)
    # Save model objects to list
    cluster_models[[i]]$null_model <- null_model
    cluster_models[[i]]$full_model <- full_model
    cluster_models[[i]]$model_lrt <- model_lrt
    cluster_models[[i]]$confint <- contrast_ci
  }

  # Organize results into output dataframe
  output <- data.frame(cluster = attributes(designmat)$dimnames[[2]],
                       size = colSums(designmat))
  output$model.pvalue <- sapply(cluster_models, function(x) x$model_lrt[["Pr(>Chisq)"]][2])
  output[[paste(contrast_lvl2, "OR", sep = ".")]] <- sapply(cluster_models, function(x) exp(fixef(x$full)[[contrast_lvl2]]))
  output[[paste(contrast_lvl2, "OR", "95pct.ci.lower", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "2.5 %"]))
  output[[paste(contrast_lvl2, "OR", "95pct.ci.upper", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "97.5 %"]))

  # Return MASC results and save models if specified
  if (save_models == TRUE) {
    saveModelObj(cluster_models, save_dir = save_model_dir)
    return(output)
  } else {
    return(output)
  }
}


meta_data <- read.table("cell_frequency.tsv", sep = "\t", header = T)
meta_data <- meta_data %>% filter(Disease == "CD")

meta_data <- meta_data %>% mutate(Group =case_when(MM_scaled <= 6.5 ~ "NI", MM_scaled > 6.5 ~ "I"))
meta_data <- meta_data %>% drop_na(Group)
meta_data$Group <- factor(meta_data$Group, levels = c("NI", "I"))
meta_data <- meta_data %>% filter(sub_bucket == "Ileal_epithelium")

var_name = "Group"

var_list  = unique(meta_data$sub_bucket)
plot_list = list()
temp <- data.frame()

###add disease duration and site
disease_duration <- read_csv('Disease_duration.csv')
disease_duration <- as.data.frame(disease_duration)
meta_data <- left_join(meta_data, disease_duration, by = "Patient")

meta_data$nUMI <- scale(log2(meta_data$total_counts))

for (ll in var_list){
    a = meta_data %>% filter(sub_bucket == ll)
    masc_df_status <- data.frame(a, clus = a$final_analysis)
    masc_res <- MASC.me(masc_df_status, factor(masc_df_status$clus),
                contrast = var_name,
                random_effects = c("Patient/sample_id"), #add Patient
                fixed_effects = c("Age","Gender","Treatment","Disease_duration","pct_counts_mt"), ##n_genes_by_counts
                verbose = TRUE,
                save_models = F) %>%
        dplyr::mutate(BH = p.adjust(model.pvalue, method = "BH")) %>%
        dplyr::arrange(model.pvalue)
    
    masc_res$sub_bucket <- ll
    
    temp = rbind(temp,masc_res)

    b <- length(unique(a$final_analysis))

    g1 <- ggplot(data = masc_res, aes(x = GroupI.OR, y = -log10(BH))) +
            theme_classic() + ylab("-log(p-value)") + xlab("Odds ratio non-inflamed vs. inflamed") +
            geom_vline(xintercept = 1, linetype = "dashed", color = "darkblue") +
            geom_errorbarh(aes(xmin=GroupI.OR.95pct.ci.lower, xmax=GroupI.OR.95pct.ci.upper), col = "darkgrey") +
            geom_point() + 
            geom_hline(yintercept = -log10(.05), linetype = "dashed", color = "darkred") + 
            scale_x_log10(breaks = c(0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8)) + geom_label_repel(size = 3, label=rownames(masc_res))
            
    plot_list[[ll]] = g1
}

pdf("MASC_inflammation_iCD.pdf", onefile = T, width = 30, height = 20)
for (ll in var_list) {
    print(plot_list[[ll]])
}
dev.off()

write.csv(temp, 'MASC_inflammation_iCD.csv')