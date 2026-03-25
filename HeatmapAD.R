library(pheatmap)
library(DESeq2) # Highly recommended for normalization
library(org.Hs.eg.db)
library(dplyr)
library(stringr)

counts_data <- read.table(r"(C:\Users\art27\OneDrive\Documents\Bioinformatics Code\Data\AD\GSE104704_RNA-Seq_Table.txt\GSE104704_raw_counts_GRCh38.p13_NCBI (2).tsv)", header=TRUE, sep="\t", row.names=1)

# Auto replace Id numbers with gene symbols
current_ids <- rownames(counts_data)

gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys = current_ids, 
                       column = "SYMBOL", 
                       keytype = "ENTREZID", 
                       multiVals = "first")


#Handle duplicates
counts_df <- as.data.frame(counts_data)
counts_df$Symbol <- gene_symbols


counts_final <- counts_df %>%
  group_by(Symbol) %>%
  summarise(across(everything(), sum)) %>%
  filter(!is.na(Symbol)) # Remove rows where mapping failed


# Set unique symbols as row names
counts_final <- as.data.frame(counts_final)
rownames(counts_final) <- counts_final$Symbol
counts_final$Symbol <- NULL # Remove the extra column


# Removed outlier group
counts_final = subset(counts_final, select = -c(GSM2806239, GSM2806221, GSM2806231) )

# Rename young, avoiding regex like the plague, yes I know this is inefficient
counts_final <- counts_final %>% rename(Young_1 = GSM2806214)
counts_final <- counts_final %>% rename(Young_2 = GSM2806215)
counts_final <- counts_final %>% rename(Young_3 = GSM2806216)
counts_final <- counts_final %>% rename(Young_4 = GSM2806217)
counts_final <- counts_final %>% rename(Young_5 = GSM2806218)
counts_final <- counts_final %>% rename(Young_6 = GSM2806219)
counts_final <- counts_final %>% rename(Young_7 = GSM2806220)
# Rename old, avoiding regex like the plague, yes I know this is inefficient
counts_final <- counts_final %>% rename(Old_1 = GSM2806222)
counts_final <- counts_final %>% rename(Old_2 = GSM2806223)
counts_final <- counts_final %>% rename(Old_3 = GSM2806224)
counts_final <- counts_final %>% rename(Old_4 = GSM2806225)
counts_final <- counts_final %>% rename(Old_5 = GSM2806226)
counts_final <- counts_final %>% rename(Old_6 = GSM2806227)
counts_final <- counts_final %>% rename(Old_7 = GSM2806228)
counts_final <- counts_final %>% rename(Old_8 = GSM2806229)
counts_final <- counts_final %>% rename(Old_9 = GSM2806230)
# Rename old_diseased, avoiding regex like the plague, yes I know this is inefficient
counts_final <- counts_final %>% rename(AD_1 = GSM2806232)
counts_final <- counts_final %>% rename(AD_2 = GSM2806233)
counts_final <- counts_final %>% rename(AD_3 = GSM2806234)
counts_final <- counts_final %>% rename(AD_4 = GSM2806235)
counts_final <- counts_final %>% rename(AD_5 = GSM2806236)
counts_final <- counts_final %>% rename(AD_6 = GSM2806237)
counts_final <- counts_final %>% rename(AD_7 = GSM2806238)
counts_final <- counts_final %>% rename(AD_8 = GSM2806240)
counts_final <- counts_final %>% rename(AD_9 = GSM2806241)
counts_final <- counts_final %>% rename(AD_10 = GSM2806242)
counts_final <- counts_final %>% rename(AD_11 = GSM2806243)


# Log 2 normalization
# Add 1 avoids log(0) errors
log_counts <- log2(counts_final + 1)


top_genes <- head(order(rowVars(as.matrix(log_counts)), decreasing=TRUE), 50)
plot_matrix <- log_counts[top_genes, ]

# Z-score scale
pheatmap(plot_matrix, 
         scale = "row",          # Normalizes each gene's expression
         clustering_distance_rows = "correlation", 
         main = "Heatmap of Top 50 DEGS in Alzheimer's")

