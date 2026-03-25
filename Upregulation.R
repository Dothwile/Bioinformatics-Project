DE_results <- read.table(file="GSE191230.top.table_primary_vs_metastasis.tsv", sep="\t", header=T) # read in the table
sig_down_genes <- DE_results %>% filter(pvalue < 0.05 & log2FoldChange > 1) %>% pull(Symbol) # filter for significant genes
write(sig_down_genes, file = "top_down_genes.tsv", sep="\n") # output a list of significant genes