if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install core processing packages
BiocManager::install("oligo") 
BiocManager::install("affy")
BiocManager::install(c("AnnotationDbi", "hgu133plus2.db"))

BiocManager::install("affycoretools")
library(affycoretools)

# Install visualization packages
BiocManager::install("ComplexHeatmap")
# or install.packages("pheatmap") if you prefer the pheatmap package

library(oligo)
library(Biobase) # required for the exprs() function
library(hgu133plus2.db)
library(AnnotationDbi)


# Set your working directory to the folder containing your .cel files
setwd(r"(C:\Users\art27\OneDrive\Documents\Bioinformatics Code\Data\Stroke)") 

# Get a list of all .cel files in the directory
celFiles <- list.celfiles(full.names = TRUE, listGzipped = TRUE)


# Read the .cel files into an AffyBatch object
rawData <- read.celfiles(celFiles)


# Perform RMA normalization to obtain a gene expression matrix
expressionData <- rma(rawData)

expressionData <- annotateEset(expressionData, hgu133plus2.db)


# Extract the expression matrix (log2 expression values)
expressionMatrix <- exprs(expressionData)


ids <- mapIds(hgu133plus2.db, keys = rownames(expressionMatrix), 
              column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
keep <- !is.na(ids)
mat <- expressionMatrix[keep, ]
ids <- ids[keep]
mat_gene <- aggregate(mat, by = list(ids), FUN = mean)
rownames(mat_gene) <- mat_gene$Group.1
mat_gene <- as.matrix(mat_gene[, -1]) # Remove the Group column


#Rename columns
colnames(mat_gene) <- substr(colnames(mat_gene), 1, nchar(colnames(mat_gene)) - 7)


# This highlights relative changes rather than absolute expression values
scaledMatrix <- t(scale(t(mat_gene))) 


# Calculate variance for each gene
geneVars <- apply(scaledMatrix, 1, var)
# Select top 50 most variable genes
topGenes <- names(sort(geneVars, decreasing = TRUE))[1:50]
heatmapData <- scaledMatrix[topGenes, ]

library(pheatmap)

# Create the heatmap
pheatmap(heatmapData, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         show_rownames = TRUE, # Hide individual gene names if too many
         main = "Heatmap of Top 50 DEGS in Stroke")

