#Get coding snps for all genes interacting with a protein for one viral taxa#
library(tidyverse)
library(data.table)
library(stringr)
virus_gene_dir <- '/ix/djishnu/Priyamvada/viral-human-ppi/taxa_based_analysis/High_confidence_interactions_only/Genes_for_virus' #directory with the gene list for all virus 
virus_gene_files <- list.files(virus_gene_dir, pattern = '.txt')
snp_annot_dir <- '/ix/djishnu/Priyamvada/viral-human-ppi/1000_genome/bisque/output_files/hg38/processed_with_gene/updated_analysis_sasa' #directory with annotated snps
snp_annot_files <- list.files(snp_annot_dir, pattern = 'coding.variants.txt')
output_dir <- '/ix/djishnu/Priyamvada/viral-human-ppi/taxa_based_analysis/High_confidence_interactions_only/Genes_for_virus/snp_to_gene_annotated'
output_dir_gene_no_snps <- '/ix/djishnu/Priyamvada/viral-human-ppi/taxa_based_analysis/Genes_for_virus/snp_to_gene_annotated/no_gene_mapped'
output_file_names <- paste0(str_remove(virus_gene_files, '.txt'), '_snps_to_gene_mapped.txt')
output_file_names_not_mapped <- paste0(str_remove(virus_gene_files, '.txt'), 'genes_no_snps.txt')
for (vir_fam in 1:length(virus_gene_files)) {
  vir_fam_genes_df <- fread(paste0(virus_gene_dir, '/', virus_gene_files[vir_fam]), header = F, col.names = c('GENE', 'TAXA_ID'), sep = '\t')
  vir_fam_gene_snp_df <- data.frame('HUMAN_PROTEIN' = character(), 'VIRAL_PROTEIN' = character(),
                                    'AA_POS' = numeric(), 'SNP_POS' = character(), 'EFFECT' = character(),
                                    'RSID' = character(), 'INTERFACE_ANNOT' = character(), 'SASA_ANNOT' = character(),
                                    'GRANTHAM_ANNOT' = character(), 'DOOLITTLE_ANNOT' = character(), 'GENE' = character(), 'TAXA_ID' = character())
  for (chr in 1:length(snp_annot_files)) {
    annot_df <- fread(paste0(snp_annot_dir, '/', snp_annot_files[chr]))
    gene_snp_chr_df <- left_join(vir_fam_genes_df, annot_df, by = 'GENE')
    gene_snp_chr_df <- gene_snp_chr_df %>% drop_na()
    vir_fam_gene_snp_df <- rbind(vir_fam_gene_snp_df, gene_snp_chr_df)
  }
  genes_not_mapped <- setdiff(vir_fam_genes_df$GENE, vir_fam_gene_snp_df$GENE)
  #vir_fam_gene_snp_df <- vir_fam_gene_snp_df[, c('GENE', 'VIRUS', "SNP_ID", 'RSID', 'ANNOTATION', 'UNIPROT_ID')]
  write.table(vir_fam_gene_snp_df, file = paste0(output_dir, '/', output_file_names[vir_fam]), col.names = T, row.names = F, sep = '\t', quote = F)
  write.table(genes_not_mapped, file = paste0(output_dir_gene_no_snps, '/', output_file_names_not_mapped[vir_fam]), row.names = T, quote = F, sep = '\t')
}
