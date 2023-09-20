#Compare FSt for top 10 or 5 percent
library(stringr)
library(Biostrings)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggprism)
library(data.table)
source('/ix/djishnu/Priyamvada/viral-human-ppi/viral-human-ppi-R-scripts/ppi-virus-human-interface/functions_analysis_Fst.R')
args <- commandArgs(trailingOnly = T)
print(args) 
fst_analysis_dir <- args[1]
pops <- args[2]
percent <- as.numeric(args[3])*100
output_dir <- args[4]
fst_cut_off_df <- read.table(paste0(args[5], '/', pops, '/', pops, '_coding.variants..fst.cut.off.txt'), header = T)
sig_cut_off <- as.numeric(args[6])
viral_fams <- list.files(paste0(fst_analysis_dir, '/', pops), pattern = '.txt') 
pred_file <- fread('/ix/djishnu/Priyamvada/viral-human-ppi/data_from_Le/v2h_binary_list_20230908_human_viral_gene_names.txt')
#removing interactions with no prediction interfaces
pred_file <- pred_file[Human_ires != '-',]
pred_file <- pred_file[Viral_ires != '-',]
ppi_annot_dir <- '/ix/djishnu/Priyamvada/viral-human-ppi/all_sasa_annotated'
interface_df <- data.frame('Top_I' = numeric(), 'Top_NI' = numeric(), 'I_p_val' = numeric(), 'Virus' = numeric())
grantham_df <- data.frame('Top_High_Grantham' = numeric(), 'Top_Low_Grantham' = numeric(), 'Grantham_p_val' = numeric(), 'Virus' = numeric())
doolittle_df <- data.frame('Top_High_Doolittle' = numeric(), 'Top_Low_Doolittle' = numeric(), 'Doolittle_p_val' = numeric(), 'Virus' = numeric())
surface_df <- data.frame('Top_S' = numeric(), 'Top_NS' = numeric(), 'Surface_p_val' = numeric(), 'Virus' = numeric())
for (vir in 1:length(viral_fams)) {
  #get the protein interaction information for viral taxa 
  vir_name <- as.numeric(str_remove(viral_fams[vir], paste0('_', pops, '_', '_fst_vals.txt')))
  pred_vir <- pred_file[Viral_taxonomy_id == vir_name,]
  if (nrow(pred_vir) > 1) {
  ppi_pairs <- paste0(pred_vir$Human_protein_uid, '_'      , pred_vir$Viral_protein_uid)
  ppi_anot_df <- create_human_prot_annot(ppi_lists = ppi_pairs, ppi_dir = ppi_annot_dir)
  fst_res_file <- list.files(viral_fams[vir], pattern = '.txt')
  fst_res_df <- fread(paste0(fst_analysis_dir, '/', pops, '/', viral_fams[vir]))
  fst_res_df <- fst_res_df[, c('HUMAN_PROTEIN', 'RSID', 'AA_POS', 'INTERFACE_ANNOT', 'SASA_ANNOT', 'GRANTHAM_ANNOT', 'DOOLITTLE_ANNOT', 'WEIR_AND_COCKERHAM_FST')]
  fst_res_df <- process_fst_df(fst_df = fst_res_df)
  fst_cut_off <- fst_cut_off_df[which(fst_cut_off_df$Cut_off_threshold == percent),]$Fst_val
  top_fst_df <- fst_res_df[which(fst_res_df$WEIR_AND_COCKERHAM_FST > fst_cut_off),]
  interface_df <- rbind(interface_df, prop_test_annot(annot_df = ppi_anot_df, annot_col_name = 'interface', top_fst_df = top_fst_df, vir = vir_name))
  grantham_df <- rbind(grantham_df, prop_test_annot(annot_df = ppi_anot_df, annot_col_name = 'grantham', top_fst_df = top_fst_df, vir = vir_name))
  doolittle_df <- rbind(doolittle_df, prop_test_annot(annot_df = ppi_anot_df, annot_col_name = 'doolittle', top_fst_df = top_fst_df, vir = vir_name))
  surface_df <- rbind(surface_df, prop_test_annot(annot_df = ppi_anot_df, annot_col_name = 'surface', top_fst_df = top_fst_df, vir = vir_name))
  } else {
    print(paste0(vir_name, 'has no predictions'))
  }
} 
if(dir.exists(paste0(output_dir, '/', pops)) == FALSE) {
  dir.create(paste0(output_dir, '/', pops))
  print("Pop directory does not exist")} else {
    print("Pop fam directory exists")
  }  
if(dir.exists(paste0(output_dir, '/', pops, '/', percent)) == FALSE) {
    dir.create(paste0(output_dir, '/', pops, '/', percent))
    print("Percentage directory does not exist")} else {
      print("Percentage fam directory exists")
    }
interface_df <- interface_df %>% mutate('Significance' = case_when(I_p_val < 0.05 & Top_I >= sig_cut_off ~ "***",
                                                                   TRUE ~ ""))
grantham_df <- grantham_df %>% mutate('Significance' = case_when(Grantham_p_val < 0.05 & Top_High_Grantham >= sig_cut_off ~ "***",
                                                                   TRUE ~ ""))
doolittle_df <- doolittle_df %>% mutate('Significance' = case_when(Doolittle_p_val < 0.05 & Top_High_Doolittle >= sig_cut_off ~ "***",
                                                                   TRUE ~ ""))
surface_df <- surface_df %>% mutate('Significance' = case_when(Surface_p_val < 0.05 & Top_S >= sig_cut_off ~ "***",
                                                                 TRUE ~ ""))
interface_p_vals <- interface_df
interface_p_vals$log_p_val <- -log10(interface_p_vals$I_p_val) 
order_of_labels <- as.character(interface_p_vals[order(interface_p_vals$log_p_val, decreasing = T),]$Virus)
virus_df <-  interface_df$Virus
pred_file_categories <- pred_file[Viral_taxonomy_id %in% virus_df, c('Viral_taxonomy_id', 'organism', 'Category')]
pred_file_categories <- unique(pred_file_categories)
colnames(pred_file_categories) <- c('Virus', 'organism', 'Category')
colors_categories <- c("Human_viruses" = '#1d2f6f',"Non-human_viruses_that_do_not_infect_humans" = '#1f7a8c', "Non-human_viruses_that_infect_humans" = '#bfdbf7')
analysis_df_list <- list(interface_df, grantham_df, doolittle_df, surface_df)
analysis_names <- c('INTERFACE', 'GRANTHAM', 'DOOLITTLE', 'SURFACE')
names(analysis_df_list) <- analysis_names
for (analysis in 1:length(analysis_names)) {
  df <- analysis_df_list[[analysis_names[analysis]]]
  write.table(df, file = paste0(output_dir, '/', pops, '/', percent, analysis_names[analysis], '_p_values.txt'), col.names = T, row.names = F, quote = F, sep = '\t')
  colnames(df)[3] <- 'p_val'
  df$log_p_val <- -log10(df$p_val)
  df <- left_join(df, pred_file_categories, by = 'Virus')
  df$Virus <- as.character(df$Virus)
  plot_stars <- data.frame('group1' = df$Virus, 'group2' = df$Virus, 'p.adj' = df$Significance, y.position = -log10(df$p_val))
  virus_names <- df$organism[match(order_of_labels, df$Virus)]
  png(paste0(output_dir, '/', pops, '/', percent, '/', analysis_names[analysis], '_p_values_hist.png'), height = 1000, width = 1044)
  print(ggplot(df) +
          geom_bar(aes(x=Virus, y=log_p_val, fill=Category), stat = 'identity', color = 'black', size = 0.5) +
          scale_fill_manual(values = colors_categories) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.75), 
                plot.title = element_text(face = 'bold'), strip.text = element_text(size = 14, colour = "black"),
                axis.title = element_text(face = 'bold', size = 14, color = 'black'), strip.background = element_blank(),
                strip.text.x = element_text(size = 14, color = 'black', face = 'bold'),
                axis.text.x = element_text(size = 14, color = 'black', face = 'bold'),
                axis.text.y = element_text(size = 14, color = 'black', face = 'bold')) +
          #geom_hline(yintercept=1.3, linetype="dashed", 
          #color = "red", size=0.25) + 
          add_pvalue(plot_stars, bracket.size = 0, fontface = 'bold', label.size = 6, tip.length = 0) +
          scale_x_discrete(guide = guide_axis(angle = 90),  limits = as.character(order_of_labels), labels = virus_names) +
          labs(title = paste0("Fst P-values for", " ", pops, "FOR", " ", analysis_names[analysis]), x = "VIRUS", y = "NEGATIVE LOG P-VAL"))
  dev.off()
}


