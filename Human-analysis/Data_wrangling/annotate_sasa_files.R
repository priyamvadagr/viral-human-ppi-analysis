library(tidyverse)
library(stringr)
library(data.table)
library(xlsx)
library(Biostrings)
all_sasa_files <- data.frame(file_names = list.files('/ix/djishnu/Priyamvada/viral-human-ppi/all_sasa'))
all_sasa_files$PPI_pair <- str_remove(all_sasa_files$file_names, '_sasa.txt')
all_sasa_files <- separate(all_sasa_files, PPI_pair, into = c("Prot_1", "Prot_2"), sep = "_")
all_sasa_files$PPI_pair <- str_remove(all_sasa_files$file_names, '_sasa.txt')
all_pred_files <- fread('/ix/djishnu/Priyamvada/viral-human-ppi/data_from_Le/v2h_binary_list_20230908_human_viral_gene_names.txt')
#all file names do not follow a pattern - changing to make them Human_protein_name_Viral_protein_name 
all_sasa_files <- all_sasa_files %>% mutate('Prot_1_type' = case_when(Prot_1 %in% all_pred_files$Viral_protein_uid ~ 'Viral_Protein',
                                            Prot_1 %in% all_pred_files$Human_protein_uid ~ 'Human_Protein'))
all_sasa_files <- all_sasa_files %>% mutate('Prot_2_type' = case_when(Prot_2 %in% all_pred_files$Viral_protein_uid ~ 'Viral_Protein',
                                                                      Prot_2 %in% all_pred_files$Human_protein_uid ~ 'Human_Protein'))
all_sasa_files <- all_sasa_files %>% 
  mutate('New_PPI' = case_when(Prot_1_type == 'Human_Protein' & Prot_2_type == 'Viral_Protein' ~ as.character(PPI_pair),
                               Prot_1_type == 'Viral_Protein' & Prot_2_type == 'Human_Protein' ~ paste0(Prot_2, '_', Prot_1)))
in_dir <- '/ix/djishnu/Priyamvada/viral-human-ppi/all_sasa'
out_dir <- '/ix/djishnu/Priyamvada/viral-human-ppi/all_sasa_annotated'
aa_max_asa <- read.xlsx('/ix/djishnu/Priyamvada/viral-human-ppi/AA_Max_ASA.xlsx', sheetIndex = 1, header = F)
aa_max_asa <- aa_max_asa[, c(1, 2)]
colnames(aa_max_asa) <- c('AMINO_ACID', 'MAX_ASA')
aa_code <- read.xlsx('/ix/djishnu/Priyamvada/viral-human-ppi/AA_Max_ASA.xlsx', sheetIndex = 2, header = T)
aa_code$Amino.acid <- str_to_title(aa_code$Amino.acid)
aa_code$Three.letter.code <- str_to_upper(aa_code$Three.letter.code)
colnames(aa_code) <- c('AMINO_ACID', 'AA_RESIDUE_3_CODE', 'AA_CODE')
#convert amino acid names to three letter codes 
aa_max_asa <- left_join(aa_max_asa, aa_code, by = 'AMINO_ACID')
aa_max_asa <- aa_max_asa[, c('AA_RESIDUE_3_CODE', 'AA_CODE', 'MAX_ASA')]
#using threshold of 30% to define accessible residues 
for (file in 1:length(all_sasa_files$file_names)) {
  sasa_file <- fread(paste0(in_dir, '/', all_sasa_files$file_names[file]))
  colnames(sasa_file) <- c('Prot', 'AA_RESIDUE_ID', 'AA_RESIDUE_3_CODE', 'DELTA_SASA', 'SASA')
  sasa_file <- left_join(sasa_file, aa_max_asa, by = 'AA_RESIDUE_3_CODE')
  sasa_file$REL_ASA <- sasa_file$SASA/sasa_file$MAX_ASA
  sasa_file <- sasa_file %>% mutate('INTERFACE_ANNOT' = case_when(DELTA_SASA >= 1 ~ 'I',
                                                                  DELTA_SASA < 1 ~ 'NI'))
  sasa_file <- sasa_file %>% mutate('SURFACE_ANNOT' = case_when(REL_ASA >= 0.3 ~ 'S',
                                                                REL_ASA < 0.3 ~ 'NS'))
  write.table(sasa_file, file = paste0(out_dir, '/', all_sasa_files$New_PPI[file]), 
              col.names = T, row.names = F, quote = F, sep = '\t')
}
