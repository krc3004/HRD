#######################################################################################################################################
# Fig1_FthroughI_Fig3_AthroughD.R
# Authors: Chirag Krishna, Robert Samstein, Nadeem Riaz

# plots key survival curves for immune checkpoint blockade (ICB)-treated patients with mutations in HR genes vs those without 
# Two HR gene sets are presented. For Fig1 F/G, HR_Gene_Set is BRCA1, BRCA2, PALB2, CHEK2, ATM
# For Fig 1 H/I, HR_Alternatibe_Gene_Set is BRCA1, BRCA2, PALB2, RAD50, RAD51, RAD51B, RAD51C, RAD51D, MRE11A, NBN

# Fig 1 FGHI, Fig 3 ABCD 
#######################################################################################################################################

## load required packages
library(data.table)
library(survival)
library(survminer)
library(survMisc)
library(ggplot2)

## load clinical data
impact_icb = as.data.frame(fread("~/Documents/updated_HR_figs/code_repo/IMPACT_ICB_treated.txt"))
impact_icb$Sample2 = paste("Pat", seq(nrow(impact_icb)), sep = "")
# impact_icb = impact_icb[,which(colnames(impact_icb) != "Sample")]
# colnames(impact_icb)[ncol(impact_icb)] = "Sample"
# impact_icb = impact_icb[,c("Sample", colnames(impact_icb)[2:(ncol(impact_icb)-1)])]

## df is the clinical data
## cancers can be any, but for paper should be "pan" or "HR-associated"
## mutation_type should be one of "BRCA1_Mutation", "BRCA2_Mutation", "HR_Gene_Set", "HR_Alternative_Gene_Set"
## plot_title can be anything
HR_survival = function(df, cancers, mutation_type, plot_title){

  ## first subset the clinical data on cancer types and mutation type
  if(cancers == "pan"){surv_df = df[,c("Sample", "Cancer_Type", "Time_to_death", "Alive_Dead", mutation_type)]}
  if(cancers == "HR-associated"){surv_df = df[which(df$Cancer_Type %in% c("Breast Cancer", "Ovarian Cancer", "Prostate Cancer", "Pancreatic Cancer")),
                                              c("Sample", "Cancer_Type", "Time_to_death", "Alive_Dead", mutation_type)]}
  if(! cancers %in% c("pan", "HR-associated")){surv_df = df[which(df$Cancer_Type %in% cancers),c("Sample", "Cancer_Type", "Time_to_death", "Alive_Dead", mutation_type)]}
  
  ## plot survival curves and log-rank p value
  HR_fit = do.call(survfit, list(formula = Surv(surv_df$Time_to_death, surv_df$Alive_Dead)~surv_df[,mutation_type]))
  print(ggsurvplot(HR_fit, pval = TRUE, font.x = c(14, "bold", "black"),
                 font.y = c(14, "bold", "black"), font.tickslab = c(12, "bold", "black"), size = 1.5, palette = c("blue", "red"), data = surv_df, title = plot_title))
  
  ## print number of patients with and without mutation
  print(table(surv_df[,mutation_type]))
}

## Fig 1F: HR Gene Set, HR-associated cancers
HR_survival(impact_icb, "HR-associated", "HR_Gene_Set", "HR Gene Set HR-associated cancers: Fig 1F")

## Fig 1G: HR Gene Set, pan cancer
HR_survival(impact_icb, "pan", "HR_Gene_Set", "HR Gene Set pan cancer: Fig 1G")

## Fig 1H: HR Alternative Gene Set, HR-associated cancers
HR_survival(impact_icb, "HR-associated", "HR_Alternative_Gene_Set", "HR Alternative Gene Set HR-associated cancers: Fig 1H")

## Fig 1I: HR Alternative Gene Set, HR-associated cancers
HR_survival(impact_icb, "pan", "HR_Alternative_Gene_Set", "HR Alternative Gene Set pan cancer: Fig 1I")

## Fig 3A: BRCA1 vs WT, HR-associated cancers
HR_survival(impact_icb, "HR-associated", "BRCA1_Mutation", "BRCA1 HR-associated cancers: Fig 3A")

## Fig 3B: BRCA2 vs WT, HR-associated cancers
HR_survival(impact_icb, "HR-associated", "BRCA2_Mutation", "BRCA2 pan cancer: Fig 3B")

## Fig 3C: BRCA1 vs WT, pan cancer
HR_survival(impact_icb, "pan", "BRCA1_Mutation", "BRCA1 pan cancer: Fig 3C")

## Fig 3D: BRCA2 vs WT, pan cancer
HR_survival(impact_icb, "pan", "BRCA2_Mutation", "BRCA2 pan cancer: Fig 3D")
