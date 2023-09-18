#Get human proteins interacting with viral proteins
#Analysis rerun to be done with a taxa based analysis 
#Additionally, will be making changes to do the SASA based analysis 
#######################################################################################################
library(tidyverse)
library(data.table)
prediction_file <- fread('/ix/djishnu/Priyamvada/viral-human-ppi/data_from_Le/v2h_binary_list_20230908_human_viral_gene_names.txt')
#subset to genes and virus taxon 
genes_interface <- prediction_file[, c("Human_gene_name", "Viral_taxonomy_id")]
genes_interface <- genes_interface %>% group_by(Viral_taxonomy_id) %>% unique()
virus_taxas <- unique(genes_interface$Viral_taxonomy_id) #389 total unique taxa
#get only viral_taxa with more than 10 human gene interactors 
n_genes <- 10 
unique_genes_df <- data.frame(table(genes_interface$Viral_taxonomy_id))
colnames(unique_genes_df) <- c('viral_scientific_name', 'no_unique-genes')
##keep only virus with n_unique_genes >= 10
virus_to_keep <- as.character(unique_genes_df[which(unique_genes_df$`no_unique-genes` >= 10), 1])
genes_interface <- genes_interface[which(genes_interface$Viral_taxonomy_id %in% virus_to_keep),]
genes_interface_list <- split(genes_interface, f = genes_interface$Viral_taxonomy_id)
for (virus in 1:length(virus_to_keep)) {
  df = genes_interface_list[[virus_to_keep[virus]]]
  #file_name = str_replace_all(virus_to_keep[virus], '/', '_._')
  write.table(df, file = paste0("/ix/djishnu/Priyamvada/viral-human-ppi/taxa_based_analysis/Genes_for_virus", '/', virus_to_keep[virus], ".txt"), sep = '\t', col.names = F, row.names = F, quote = F)
}
