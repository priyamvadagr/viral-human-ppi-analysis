#process all_preds_file to add Surface variant informations 
library(tidyverse)
library(data.table)
library(stringr)
library(glue)
preds_file <- fread('/ix/djishnu/Priyamvada/viral-human-ppi/data_from_Le/v2h_binary_list_20230908_human_viral_gene_names.txt')
preds_file$PPI_pair <- paste0(preds_file$Human_protein_uid, '_', preds_file$Viral_protein_uid) #SASA files saved in this order
Viral_Interface_Surgace_annotation <- c() #vectors to be added to the final pred file 
Human_Interface_Surface_annotation <- c()
for (PPI in 1:length(preds_file$PPI_pair)) {
  PPI_name <- preds_file$PPI_pair[PPI]
  if (file.exists(paste0('/ix/djishnu/Priyamvada/viral-human-ppi/all_sasa_annotated', '/',  PPI_name))) {
  SASA_file <- fread(paste0('/ix/djishnu/Priyamvada/viral-human-ppi/all_sasa_annotated', '/',  PPI_name))
  Viral_prot <- preds_file$Viral_protein_uid[PPI]
  Viral_interface_res <- data.frame(str_split(preds_file$Viral_ires[PPI], ','))
  colnames(Viral_interface_res) <- "AA_RESIDUE_ID"
  SASA_file$AA_RESIDUE_ID <- as.character(SASA_file$AA_RESIDUE_ID)
  SASA_virus <- SASA_file[Prot %in% Viral_prot, ]
  Viral_interface_res <- left_join(Viral_interface_res, SASA_virus, by = 'AA_RESIDUE_ID')
  Viral_interface_surface_annot <- glue_collapse(Viral_interface_res$SURFACE_ANNOT, ',')
  Viral_Interface_Surgace_annotation <- c(Viral_Interface_Surgace_annotation, Viral_interface_surface_annot)
  Human_prot <- preds_file$Human_protein_uid[PPI]
  Human_interface_res <- data.frame(str_split(preds_file$Human_ires[PPI], ','))
  colnames(Human_interface_res) <- 'AA_RESIDUE_ID'
  SASA_human <- SASA_file[Prot %in% Human_prot,]
  Human_interface_res <- left_join(Human_interface_res, SASA_human, by = 'AA_RESIDUE_ID')
  Human_interface_surface_annot <- glue_collapse(Human_interface_res$SURFACE_ANNOT, ',')
  Human_Interface_Surface_annotation <- c(Human_Interface_Surface_annotation, Human_interface_surface_annot)
  }
  else {
    Viral_Interface_Surgace_annotation <- c(Viral_Interface_Surgace_annotation, '-')
    Human_Interface_Surface_annotation <- c(Human_Interface_Surface_annotation, '-')
  }
}
preds_file$Viral_interface_res_surface_annot <- Viral_Interface_Surgace_annotation
preds_file$Human_interface_res_surface_annot <- Human_Interface_Surface_annotation
write.table(preds_file, file = '/ix/djishnu/Priyamvada/viral-human-ppi/data_from_Le/v2h_binary_list_20230908_surface_annot.txt', 
            col.names = T, row.names = F, sep = '\t', quote = F)