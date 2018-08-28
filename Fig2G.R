#######################################################################################################
# Fig2G.R
# Authors: Chirag Krishna, Robert Samstein, Nadeem Riaz

# compares indels, SNVs, and neopeptides for two sets of comparisons:
# BRCA1 triple negative vs all triple negative; BRCA2 non-triple negative vs all non-triple negative
#######################################################################################################

## load required packages
library(data.table)
library(ggplot2)

## load TCGA breast cancer data
bdat = as.data.frame(fread("~/Documents/updated_HR_figs/code_repo/TCGA_BreastCancer_data.txt"))

## set up data frames for grouped boxplots
## BRCA1 indels
b1TN_indel = as.data.frame(cbind(rep("Indels"), c(log10(bdat[which(bdat$BRCA1_triple_negative == 1),]$Indels+1), log10(bdat[which(bdat$all_triple_negative == 1),]$Indels+1))))
b1TN_indel$comp = c(rep("BRCA1 TN", length(bdat[which(bdat$BRCA1_triple_negative == 1),]$Id)), rep("All TN", length(bdat[which(bdat$all_triple_negative == 1),]$Id)))
## BRCA1 SNVs
b1TN_snv = as.data.frame(cbind(rep("SNVs"), c(log10(bdat[which(bdat$BRCA1_triple_negative == 1),]$SNVs+1), log10(bdat[which(bdat$all_triple_negative == 1),]$SNVs+1))))
b1TN_snv$comp = c(rep("BRCA1 TN", length(bdat[which(bdat$BRCA1_triple_negative == 1),]$Id)), rep("All TN", length(bdat[which(bdat$all_triple_negative == 1),]$Id)))
## BRCA1 Neopeptides
b1TN_neopep = as.data.frame(cbind(rep("Neopeptides"), c(log10(bdat[which(bdat$BRCA1_triple_negative == 1),]$Neopeptides+1), log10(bdat[which(bdat$all_triple_negative == 1),]$Neopeptides+1))))
b1TN_neopep$comp = c(rep("BRCA1 TN", length(bdat[which(bdat$BRCA1_triple_negative == 1),]$Id)), rep("All TN", length(bdat[which(bdat$all_triple_negative == 1),]$Id)))

## combine all data frames
b1TN_comps = do.call(rbind, list(b1TN_indel, b1TN_snv, b1TN_neopep))
b1TN_comps$V2 = as.numeric(as.character(b1TN_comps$V2))
b1TN_comps$comp = factor(b1TN_comps$comp, levels = c("BRCA1 TN", "All TN"))

## needed for boxplots
quant = function(x){
  r = quantile(x, probs =c(0.05,0.25, 0.5, 0.75,0.95))
  names(r) = c("ymin", "lower", "middle", "upper", "ymax")
  return(r)
}

## plot
ggplot(b1TN_comps, aes(V1, V2)) +
  stat_summary(fun.data = quant, aes(fill = comp), geom = "boxplot", position = "dodge", color = "black", lwd = 2, outlier.shape=NA) +
  scale_x_discrete() +
  ylim(0,4.2) +
  scale_fill_manual(values = c("red", "blue")) +
  theme(axis.ticks = element_line(color = "black", size = 2)) +
  theme(axis.ticks.length = unit(.08, "in")) +
  theme(axis.text = element_text(color = "black", face="bold", size=16)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(legend.title=element_blank()) + 
  theme(legend.text=element_text(size=12, face = "bold")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position="none")

## print p values for all comparisons
## BRCA1 indels
wilcox.test(b1TN_comps[which(b1TN_comps$comp == "BRCA1 TN" & b1TN_comps$V1 == "Indels"),]$V2, 
            b1TN_comps[which(b1TN_comps$comp == "All TN" & b1TN_comps$V1 == "Indels"),]$V2)
## BRCA1 SNVs
wilcox.test(b1TN_comps[which(b1TN_comps$comp == "BRCA1 TN" & b1TN_comps$V1 == "SNVs"),]$V2, 
            b1TN_comps[which(b1TN_comps$comp == "All TN" & b1TN_comps$V1 == "SNVs"),]$V2)
## BRCA1 neopeptides
wilcox.test(b1TN_comps[which(b1TN_comps$comp == "BRCA1 TN" & b1TN_comps$V1 == "Neopeptides"),]$V2, 
            b1TN_comps[which(b1TN_comps$comp == "All TN" & b1TN_comps$V1 == "Neopeptides"),]$V2)


## do the same for BRCA2 NTN vs all NTN comparisons

## set up data frames for grouped boxplots
## BRCA2 indels
b2NTN_indel = as.data.frame(cbind(rep("Indels"), c(log10(bdat[which(bdat$BRCA2_non_triple_negative == 1),]$Indels+1), log10(bdat[which(bdat$all_non_triple_negative == 1),]$Indels))))
b2NTN_indel$comp = c(rep("BRCA2 NTN", length(bdat[which(bdat$BRCA2_non_triple_negative == 1),]$Id)), rep("All NonTN", length(bdat[which(bdat$all_non_triple_negative == 1),]$Id)))
## BRCA2 SNVs
b2NTN_snv = as.data.frame(cbind(rep("SNVs"), c(log10(bdat[which(bdat$BRCA2_non_triple_negative == 1),]$SNVs+1), log10(bdat[which(bdat$all_non_triple_negative == 1),]$SNVs))))
b2NTN_snv$comp = c(rep("BRCA2 NTN", length(bdat[which(bdat$BRCA2_non_triple_negative == 1),]$Id)), rep("All NonTN", length(bdat[which(bdat$all_non_triple_negative == 1),]$Id)))
## BRCA2 Neopeptides
b2NTN_neopep = as.data.frame(cbind(rep("Neopeptides"), c(log10(bdat[which(bdat$BRCA2_non_triple_negative == 1),]$Neopeptides+1), log10(bdat[which(bdat$all_non_triple_negative == 1),]$Neopeptides))))
b2NTN_neopep$comp = c(rep("BRCA2 NTN", length(bdat[which(bdat$BRCA2_non_triple_negative == 1),]$Id)), rep("All NonTN", length(bdat[which(bdat$all_non_triple_negative == 1),]$Id)))

## combine all data frames
b2NTN_comps = do.call(rbind, list(b2NTN_indel, b2NTN_snv, b2NTN_neopep))
b2NTN_comps$V2 = as.numeric(as.character(b2NTN_comps$V2))
b2NTN_comps$comp = factor(b2NTN_comps$comp, levels = c("BRCA2 NTN", "All NonTN"))

## plot
ggplot(b2NTN_comps, aes(V1, V2)) +
  stat_summary(fun.data = quant, aes(fill = comp), geom = "boxplot", position = "dodge", color = "black", lwd = 2, outlier.shape=NA) +
  scale_x_discrete() +
  ylim(0,4.2) +
  scale_fill_manual(values = c("orange", "darkgreen")) +
  theme(axis.ticks = element_line(color = "black", size = 2)) +
  theme(axis.ticks.length = unit(.08, "in")) +
  theme(axis.text = element_text(color = "black", face="bold", size=16)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(legend.title=element_blank()) + 
  theme(legend.text=element_text(size=12, face = "bold")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.line = element_line(colour = 'black', size = 2)) +
  theme(legend.position="none")

## print p values for all comparisons
## BRCA2 indels
wilcox.test(b2NTN_comps[which(b2NTN_comps$comp == "BRCA2 NTN" & b2NTN_comps$V1 == "Indels"),]$V2, 
            b2NTN_comps[which(b2NTN_comps$comp == "All NonTN" & b2NTN_comps$V1 == "Indels"),]$V2)
## BRCA2 SNVs
wilcox.test(b2NTN_comps[which(b2NTN_comps$comp == "BRCA2 NTN" & b2NTN_comps$V1 == "SNVs"),]$V2, 
            b2NTN_comps[which(b2NTN_comps$comp == "All NonTN" & b2NTN_comps$V1 == "SNVs"),]$V2)
## BRCA2 Neopeptides
wilcox.test(b2NTN_comps[which(b2NTN_comps$comp == "BRCA2 NTN" & b2NTN_comps$V1 == "Neopeptides"),]$V2, 
            b2NTN_comps[which(b2NTN_comps$comp == "All NonTN" & b2NTN_comps$V1 == "Neopeptides"),]$V2)


wilcox.test(bdat[which(bdat$LST >= 15),]$CYT, bdat[which(bdat$LST < 15),]$CYT)
qs3 = quantile(bdat$Abs_Exposure_Sig.3, na.rm = TRUE)
bdat$s3h = ifelse(bdat$Abs_Exposure_Sig.3 >= qs3[4], 1, 0)
wilcox.test(bdat[which(bdat$s3h == 1),]$CYT, bdat[which(bdat$s3h == 0),]$CYT)


fcna_b12 = as.data.frame(c(as.numeric(res[which(res$Sample %in% bdat[which(bdat$BRCA2_pathogenic == 1),]$Id),]$facets_CNA_Genome),
                           as.numeric(res[which(res$Sample %in% bdat[which(bdat$BRCA1_pathogenic == 1),]$Id),]$facets_CNA_Genome)))
colnames(fcna_b12) = c("val")
fcna_b12 = as.data.frame(cbind(fcna_b12, c(rep("BRCA2", nrow(res[which(res$Sample %in% bdat[which(bdat$BRCA2_pathogenic == 1),]$Id),])), rep("BRCA1", nrow(res[which(res$Sample %in% bdat[which(bdat$BRCA1_pathogenic == 1),]$Id),])))))
colnames(fcna_b12) = c("val", "group")
fcna_b12$dum = rep("comp1", nrow(fcna_b12))

wilcox.test(fcna_b12[which(fcna_b12$group == "BRCA1"),]$val, fcna_b12[which(fcna_b12$group == "BRCA2"),]$val)

fcna_b12_TN = as.data.frame(c(as.numeric(res[which(res$Sample %in% bdat[which(bdat$BRCA1_triple_negative == 1),]$Id),]$facets_CNA_Genome),
                              as.numeric(res[which(res$Sample %in% bdat[which(bdat$all_triple_negative == 1),]$Id),]$facets_CNA_Genome)))
colnames(fcna_b12_TN) = c("val")
fcna_b12_TN = as.data.frame(cbind(fcna_b12_TN, c(rep("BRCA1 TN", nrow(res[which(res$Sample %in% bdat[which(bdat$BRCA1_triple_negative == 1),]$Id),])), rep("All TN", nrow(res[which(res$Sample %in% bdat[which(bdat$all_triple_negative == 1),]$Id),])))))
colnames(fcna_b12_TN) = c("val", "group")
fcna_b12_TN$dum = rep("comp2", nrow(fcna_b12_TN))

wilcox.test(fcna_b12_TN[which(fcna_b12_TN$group == "BRCA1 TN"),]$val, fcna_b12_TN[which(fcna_b12_TN$group == "All TN"),]$val)

fcna_b12_NTN = as.data.frame(c(as.numeric(res[which(res$Sample %in% bdat[which(bdat$BRCA2_non_triple_negative == 1),]$Id),]$facets_CNA_Genome),
                               as.numeric(res[which(res$Sample %in% bdat[which(bdat$all_non_triple_negative == 1),]$Id),]$facets_CNA_Genome)))
colnames(fcna_b12_NTN) = c("val")
fcna_b12_NTN = as.data.frame(cbind(fcna_b12_NTN, c(rep("BRCA2 NTN", nrow(res[which(res$Sample %in% bdat[which(bdat$BRCA2_non_triple_negative == 1),]$Id),])), rep("All NTN", nrow(res[which(res$Sample %in% bdat[which(bdat$all_non_triple_negative == 1),]$Id),])))))
colnames(fcna_b12_NTN) = c("val", "group")
fcna_b12_NTN$dum = rep("comp3", nrow(fcna_b12_NTN))

wilcox.test(fcna_b12_NTN[which(fcna_b12_NTN$group == "BRCA2 NTN"),]$val, fcna_b12_NTN[which(fcna_b12_NTN$group == "All NTN"),]$val)

fcna_combs = do.call(rbind, list(fcna_b12_TN, fcna_b12_NTN, fcna_b12))
fcna_combs$group = factor(fcna_combs$group, levels = c("BRCA1 TN", "All TN", "BRCA2 NTN", "All NTN", "BRCA1", "BRCA2"))
fcna_combs$dum = factor(fcna_combs$dum, levels = c("comp1", "comp2", "comp3"))

ggplot(fcna_combs, aes(x=dum, y=val, fill=group)) + 
  geom_boxplot() + 
  scale_x_discrete() +
  ylim(0,1) +
  scale_fill_manual(values = c("red", "blue", "orange", "darkgreen", "pink", "purple")) +
  #annotate("text", x=1, y=1000, label= "p = 0.003", fontface = "bold", size = 5) + 
  #annotate("text", x=2, y=1000, label= "p = 0.17", fontface = "bold", size = 5) + 
  #annotate("text", x=3, y=1000, label= "p = 0.06", fontface = "bold", size = 5) + 
  theme(axis.text = element_text(color = "black", face="bold", size=16)) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank()) +
  theme(legend.title=element_blank()) + 
  theme(legend.text=element_text(size=12, face = "bold")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


