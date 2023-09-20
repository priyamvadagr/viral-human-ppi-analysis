#map snps to genes#
library(tidyverse)
library(data.table)
library(stringr)
options(warn = 1)
for (i in 1:22) {
chr <- i
bisque_output_dir <- '/ix/djishnu/Priyamvada/viral-human-ppi/1000_genome/bisque/output_files/hg38'
bisque_output_files  <- list.files(bisque_output_dir, pattern = paste0("1000G.EUR.hg38", ".", chr, ".bim_[0-9]", ".txt.bisque.out.txt.bisque.all.annot.out")) #all annotations added to bisque file 
uniprot_file_dir <- "/ix/djishnu/Priyamvada/viral-human-ppi/1000_genome/bisque/output_files/hg38/uniprot_gene_mapping_output/uniprot_id_mapping_out"
uniprot_file_name <- paste0(uniprot_file_dir, '/', 'chr.', chr, '.gene.names.tsv.gz')
uniprot_gene_names <- fread(uniprot_file_name)
colnames(uniprot_gene_names) <- c('HUMAN_PROTEIN', 'GENE')
output_dir <- "/ix/djishnu/Priyamvada/viral-human-ppi/1000_genome/bisque/output_files/hg38/processed_with_gene/updated_analysis_sasa"
output_file_name <- paste0('Chr.', chr,  '.coding.variants.txt')
chr_bisque_df <- data.frame('HUMAN_PROTEIN' = character(), 'VIRAL_PROTEIN' = character(),
                            'AA_POS' = numeric(), 'SNP_POS' = character(), 'EFFECT' = character(),
                            'RSID' = character(), 'INTERFACE_ANNOT' = character(), 'SASA_ANNOT' = character(),
                            'GRANTHAM_ANNOT' = character(), 'DOOLITTLE_ANNOT' = character())
for (file in 1:length(bisque_output_files)) {
  print(bisque_output_files[file])
  bisque_output_df <- fread(paste0(bisque_output_dir, '/', bisque_output_files[file]), header = T)
  #bisque_output_df[bisque_output_df == " "] <- NA
  bisque_output_df <- drop_na(bisque_output_df)
  #colnames(bisque_output_df) <- c('SNP_ID', 'From', 'PROTEIN_POS', 'SUBSTITUTION', 'rsid')
  bisque_output_df$HUMAN_PROTEIN <- str_remove(bisque_output_df$HUMAN_PROTEIN, '-[0-9]') #had removed when mapping ids in uniprot 
  bisque_output_df <- left_join(bisque_output_df, uniprot_gene_names, by = 'HUMAN_PROTEIN', relationship = 'many-to-many')
  #colnames(bisque_output_df)[6] <- 'GENE'
  chr_bisque_df <- rbind(chr_bisque_df, bisque_output_df)
}
write.table(chr_bisque_df, paste0(output_dir, '/', output_file_name), col.names = T, quote = F, row.names = F, sep = '\t')
}

