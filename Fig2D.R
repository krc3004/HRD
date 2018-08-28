##########################################################################################
# Fig2D.R
# Authors: Chirag Krishna, Robert Samstein, Nadeem Riaz

# plot PCAs in figure 2D, comparing immune microenvironment composition of BRCA1-mutant 
# and BRCA2-mutant samples as measured by ssGSEA scores
##########################################################################################

## load required packages
library(data.table)

## load TCGA breast cancer data
bdat = as.data.frame(fread("~/Documents/updated_HR_figs/code_repo/TCGA_BreastCancer_data.txt"))

## df is the TCGA table, which contains the sample IDs, BRCA1/2 mutation status, and ssGSEA scores (among other columns)
## cell_type can be one of "innate", "tcell", "both" (combines innate and t cell)
## plot_title can be anything
HR_PCA = function(df, cell_type, plot_title){
  
  ## first subset the table on all of the BRCA1 and BRCA2-mutant cases
  mutants = df[which(df$BRCA1_pathogenic == 1 | df$BRCA2_pathogenic == 1),]
  
  ## set rownames, useful for plotting later
  rownames(mutants) = mutants$Id
  
  ## pick the cell type for PCA
  if(cell_type == "innate"){col_selection = colnames(df)[8:21]} ## gets the innate cell types
  if(cell_type == "tcell"){col_selection = colnames(df)[c(19,22:29)]} ## gets the tcell types. Note gamma delta is also innate!
  if(cell_type == "both"){col_selection = colnames(df)[8:29]} ## all cells in the table (combines innate and tcell)
  
  ## now compute PCs
  cell_pca = prcomp(mutants[,which(colnames(mutants) %in% col_selection)], scale. = TRUE)
  
  ## pull out the PCs
  cell_pca_PCs = as.data.frame(cell_pca$x)
  
  ## assign IDs (helps with plotting)
  cell_pca_PCs$id = rownames(cell_pca_PCs)
  
  ## assign BRCA1 or BRCA2 mutation status based on ID (helps with plotting)
  cell_pca_PCs$status = apply(cell_pca_PCs, 1, function(x) ifelse(x["id"] %in% df[which(df$BRCA1_pathogenic == 1),]$Id, "BRCA1", "BRCA2"))
  
  ## assign colors based on BRCA1 or BRCA2 mutation status
  ## BRCA2 is blue, BRCA1 is red
  cell_pca_PCs$color = apply(cell_pca_PCs, 1, function(x) ifelse(x["status"] == "BRCA2", "blue", "red"))
  
  ## plot PC1 vs PC2
  ## can also plot any other pair of PCs or set of PCs e.g. cell_pca_PCS[,1:4] for top 4 PCs
  plot(cell_pca_PCs[,1:2], pch = 21, cex = 2.2, bg = cell_pca_PCs$color, col = "black", upper.panel = NULL, xaxt = "n", yaxt = "n", xlab = "PC1", ylab = "PC2", main = plot_title)
}

## Innate
HR_PCA(bdat, "innate", "Innate")

## Tcell
HR_PCA(bdat, "tcell", "TCell")

## Both
HR_PCA(bdat, "both", "Both")
