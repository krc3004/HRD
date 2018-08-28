##########################################################################################
# Fig2E.R
# Authors: Chirag Krishna, Robert Samstein, Nadeem Riaz

# plot heatmaps in figure 2E, comparing immune microenvironment composition of BRCA1-mutant 
# and BRCA2-mutant samples as measured by ssGSEA scores
##########################################################################################

## load required packages
library(data.table)
library(pheatmap)
library(RColorBrewer)

## load TCGA breast cancer data
bdat = as.data.frame(fread("~/Documents/updated_HR_figs/code_repo/TCGA_BreastCancer_data.txt"))

## subset on BRCA1 and BRCA2-mutant samples
hmap = bdat[which(bdat$BRCA1_pathogenic == 1 | bdat$BRCA2_pathogenic == 1),]

## set rownames of matrix (need for plotting heatmap)
rownames(hmap) = hmap$Id

## keep columns with cell type
hmap = hmap[,10:29]

## scale across both BRCA1 and BRCA2-mutant samples
## plotting separately, but normalizing across both types of samples
for(i in seq(ncol(hmap))){
  hmap[,i] = (hmap[,i] - mean(hmap[,i]))/sd(hmap[,i])
}

## transpose for plotting
hmap = as.data.frame(t(hmap))

## split BRCA2s and BRCA1s for separate heatmaps
hmap_brca1 = hmap[,which(colnames(hmap) %in% bdat[which(bdat$BRCA1_pathogenic == 1),]$Id)]
hmap_brca2 = hmap[,which(colnames(hmap) %in% bdat[which(bdat$BRCA2_pathogenic == 1),]$Id)]

## reorder rownames so they appear in same exact order as in paper
correct_row_order = c("CD8 T cells", "Tcm cells", "Tem cells", "Th1 cells", "Th17 cells", "Th2 cells", "Treg cells", 
                  "aDC", "DC", "Eosinophils", "iDC", "Macrophages", "Mast cells", "Neutrophils", 
                  "NK CD56bright cells", "NK CD56dim cells", "NK cells", "pDC", "Tfh cells", "Tgd cells")
hmap_brca1 = hmap_brca1[correct_row_order,]
hmap_brca2 = hmap_brca2[correct_row_order,]

## set color range
hmap_color_range = seq(-6, 4, by = 1)

## plot heatmaps
## BRCA1
pheatmap(hmap_brca1, scale = "none", cluster_rows = FALSE, cluster_cols = TRUE, treeheight_col = 0, treeheight_row = 0, 
         color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(hmap_color_range)), 
         breaks = hmap_color_range, show_colnames = TRUE, main = "BRCA1")

## BRCA2
pheatmap(hmap_brca2, scale = "none", cluster_rows = FALSE, cluster_cols = TRUE, treeheight_col = 0, treeheight_row = 0, 
         color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(hmap_color_range)), 
         breaks = hmap_color_range, show_colnames = TRUE, main = "BRCA2")