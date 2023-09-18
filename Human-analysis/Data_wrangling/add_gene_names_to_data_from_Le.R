#Add gene names to the predicted dataset 
library(data.table)
library(biomaRt)
library(tidyverse)
pred_data <- fread("/ix/djishnu/Priyamvada/viral-human-ppi/v2h_binary_list_with_ires_iseq_category_seq_20230908.txt")
human_proteins <- unique(pred_data$Human_protein_name)
viral_proteins <- unique(pred_data$Viral_protein_name)
write.table(human_proteins, file = '/ix/djishnu/Priyamvada/viral-human-ppi/data_from_Le/v2h_binary_list_20230908_human_proteins.txt', col.names = F, 
            row.names = F, quote = F, sep = '\t')
write.table(viral_proteins, file = '/ix/djishnu/Priyamvada/viral-human-ppi/data_from_Le/v2h_binary_list_20230908_viral_proteins.txt', col.names = F, 
            row.names = F, quote = F, sep = '\t')
human_protein_map <- fread('/ix/djishnu/Priyamvada/viral-human-ppi/data_from_Le/uniprot_human_protein_name_to_gene_map.tsv.gz')
colnames(human_protein_map) <- c('Human_protein_name', 'Human_gene_name')
viral_protein_map <- fread('/ix/djishnu/Priyamvada/viral-human-ppi/data_from_Le/uniprot_viral_protein_name_to_gene_map.tsv.gz')
colnames(viral_protein_map) <- c('Viral_protein_name', 'Viral_gene_name')
#add viral and human gene names 
pred_data <- left_join(pred_data, human_protein_map, by = 'Human_protein_name')
pred_data <- left_join(pred_data,viral_protein_map, by = 'Viral_protein_name')
write.table(pred_data, file = '/ix/djishnu/Priyamvada/viral-human-ppi/data_from_Le/v2h_binary_list_20230908_human_viral_gene_names.txt', col.names = T, 
            row.names = F, quote = F, sep = '\t')