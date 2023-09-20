#get Fst values for SNPs - SNPs in the form of chr:pos:alt:ref#
library(stringr)
library(Biostrings)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(data.table)
args <- commandArgs(trailingOnly = T)
print(args)
pops <- args[1]
fst_results_dir <- args[2]
fst_results_files <- list.files(paste0(fst_results_dir, '/', pops)) 
viral_snp_genes_dir <- args[3]
viral_snp_genes_files <- paste0(args[4], '_snps_to_gene_mapped.txt')
output_dir <- args[5] 

#Get data frame with snps mapped to gene
  viral_fam_name <- str_remove(viral_snp_genes_files, '_snps_to_gene_mapped.txt')
  viral_fam_snp_genes_df <- fread(paste0(viral_snp_genes_dir, '/', viral_snp_genes_files), header = T)
  colnames(viral_fam_snp_genes_df)[which(colnames(viral_fam_snp_genes_df) == 'SNP_POS')] <- 'SNP'
  viral_fam_snp_genes_df <- unique(viral_fam_snp_genes_df)
  viral_fam_snp_genes_df_snp_id_reversed <- viral_fam_snp_genes_df
  viral_fam_snp_genes_df_snp_id_reversed$SNP <- paste0('chr', viral_fam_snp_genes_df_snp_id_reversed$SNP)
  str_sub(viral_fam_snp_genes_df_snp_id_reversed$SNP, -3) <- Biostrings::reverse(str_sub(viral_fam_snp_genes_df_snp_id_reversed$SNP, -3)) #flipping because alt and ref do not match for certain rsids 
  viral_fam_snp_genes_df$SNP <- paste0('chr', viral_fam_snp_genes_df$SNP)
  viral_fam_snp_genes_df <- rbind(viral_fam_snp_genes_df, viral_fam_snp_genes_df_snp_id_reversed)
  snp_fst_df <- data.frame('HUMAN_PROTEIN' = character(), 'VIRAL_PROTEIN' = character(),
                           'AA_POS' = numeric(), 'SNP_POS' = character(), 'EFFECT' = character(),
                           'RSID' = character(), 'INTERFACE_ANNOT' = character(), 'SASA_ANNOT' = character(),
                           'GRANTHAM_ANNOT' = character(), 'DOOLITTLE_ANNOT' = character(), 'GENE' = character(),
                           'TAXA_ID' = character(), 'FST_VAL' = numeric())
  for (file in 1:length(fst_results_files)) {
    chr <- fread(paste0(fst_results_dir, '/', pops, '/', fst_results_files[file]))
    chr_snp_fst_df <- left_join(viral_fam_snp_genes_df, chr, by = 'SNP')
    chr_snp_fst_df <- drop_na(chr_snp_fst_df)
    snp_fst_df <- rbind(snp_fst_df, chr_snp_fst_df)
  }
  if(dir.exists(paste0(output_dir, '/', pops)) == FALSE) {
    dir.create(paste0(output_dir, '/', pops))
    print("Population directory does not exist")} else {
      print("Population directory exists")
    } 
  write.table(snp_fst_df, file = paste0(output_dir, '/', pops, '/', viral_fam_name , '_', pops, '_', '_fst_vals.txt'), quote = F, sep = '\t', col.names = T, row.names = F)


