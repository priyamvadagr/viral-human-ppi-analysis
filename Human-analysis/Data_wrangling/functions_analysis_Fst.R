#function to get all ppi pair for a virus and create one annote file : if residue is considered I for one but NI for another we will be keeping it as I 
create_human_prot_annot <- function(ppi_lists, ppi_dir) {
  annot_df <- data.frame()
  total_length <- c()
  for (ppi in 1:length(ppi_lists)) {
    human_prot <- str_split(ppi_lists[ppi], '_')[[1]][1]
    annot_file <- fread(paste0(ppi_dir, '/', ppi_lists[ppi]))
    annot_file <- annot_file[Prot == human_prot,]
    annot_file <- annot_file[, c('Prot','AA_RESIDUE_ID', 'INTERFACE_ANNOT', 'SURFACE_ANNOT', 'GRANTHAM_CATEGORY', 'DOOLITTLE_CATEGORY')]
    prot_length <- nrow(annot_file)
    names(prot_length) <- human_prot
    annot_df <- rbind(annot_df, annot_file)
    total_length <- c(total_length, prot_length)
  }
  annot_df <- annot_df[!duplicated(annot_df[,c('Prot','AA_RESIDUE_ID', 'INTERFACE_ANNOT')]),]
  #check if any residue and protein is repeated - if yes check for Interface annot 
  rows_to_remove <- annot_df %>% 
    group_by(Prot, AA_RESIDUE_ID) %>%
    filter(n() > 1)
   x <- annot_df %>% 
    group_by(Prot, AA_RESIDUE_ID) %>%
    filter(n() > 1) %>% filter(INTERFACE_ANNOT == 'I') %>% ungroup()
  annot_df <- anti_join(annot_df, rows_to_remove)
  annot_df <- rbind(annot_df, x)
  colnames(annot_df) <- c('Prot','AA_RESIDUE_ID', 'INTERFACE_ANNOT', 'SURFACE_ANNOT', 'GRANTHAM_ANNOT', 'DOOLITTLE_ANNOT')
  return(annot_df)
}

#process fst df 
#keeping unique protein and rsid pairs 
process_fst_df <- function(fst_df) {
    fst_res_df <- fst_res_df[!duplicated(fst_res_df[, c('HUMAN_PROTEIN', 'RSID', 'AA_POS', 'INTERFACE_ANNOT')]),]
    rows_to_remove <- fst_res_df %>% 
    group_by(HUMAN_PROTEIN, RSID) %>%
    filter(n() > 1)
    x <- fst_res_df %>% 
            group_by(HUMAN_PROTEIN, RSID) %>%
              filter(n() > 1) %>% filter(INTERFACE_ANNOT == 'I') %>% ungroup()
    fst_res_df <- anti_join(fst_res_df, rows_to_remove)
    fst_res_df <- rbind(fst_res_df, x)
    return(fst_res_df)
}

#function to get prop.test 
prop_test_annot <- function(annot_df, annot_col_name, top_fst_df, vir) {
  if (annot_col_name == 'interface') {
    annot_df <- annot_df %>% mutate('INTERFACE_ANNOT_MODIFIED' = case_when((INTERFACE_ANNOT == 'I' & SURFACE_ANNOT == 'S') ~ 'I',
                                                                                (INTERFACE_ANNOT == 'I' & SURFACE_ANNOT == 'NS') ~ 'NI',
                                                                                (INTERFACE_ANNOT == 'NI' & SURFACE_ANNOT == 'NS') ~ 'NI',
                                                                                (INTERFACE_ANNOT == 'NI' & SURFACE_ANNOT == 'S') ~ 'NI'))
    total_NI <- length(which(annot_df$INTERFACE_ANNOT_MODIFIED == 'NI')) #for interface analysis I == I & S and NI == all other
    total_I <- length(which(annot_df$INTERFACE_ANNOT_MODIFIED == 'I'))
    top_fst_df <- top_fst_df %>% mutate('INTERFACE_ANNOT_MODIFIED' = case_when((INTERFACE_ANNOT == 'I' & SASA_ANNOT == 'S') ~ 'I',
                                                                             (INTERFACE_ANNOT == 'I' & SASA_ANNOT == 'NS') ~ 'NI',
                                                                             (INTERFACE_ANNOT == 'NI' & SASA_ANNOT == 'NS') ~ 'NI',
                                                                             (INTERFACE_ANNOT == 'NI' & SASA_ANNOT == 'S') ~ 'NI'))
    top_I <- length(which(top_fst_df$INTERFACE_ANNOT_MODIFIED == 'I'))
    top_NI <- length(which(top_fst_df$INTERFACE_ANNOT_MODIFIED == 'NI'))
    successes <- c(top_I, top_NI)
    trials <- c(total_I, total_NI)
    prop_test_p_val <- prop.test(successes, trials, alternative = "greater")$p.value
    result_df <- data.frame('Top_I' = top_I, 'Top_NI' = top_NI, 'I_p_val' = prop_test_p_val, 'Virus' = vir)
    return(result_df)
  }
  else if (annot_col_name == 'grantham') {
    total_NI <- length(which(annot_df$GRANTHAM_ANNOT == 'LOW')) #for interface analysis I == I & S and NI == all other
    total_I <- length(which(annot_df$GRANTHAM_ANNOT == 'HIGH'))
    top_I <- length(which(top_fst_df$GRANTHAM_ANNOT == 'HIGH'))
    top_NI <- length(which(top_fst_df$GRANTHAM_ANNOT == 'LOW'))
    successes <- c(top_I, top_NI)
    trials <- c(total_I, total_NI)
    prop_test_p_val <- prop.test(successes, trials, alternative = "greater")$p.value
    result_df <- data.frame('Top_High_Grantham' = top_I, 'Top_Low_Grantham' = top_NI, 'Grantham_p_val' = prop_test_p_val, 'Virus' = vir)
    return(result_df)

  } else if (annot_col_name == 'doolittle') {
    total_NI <- length(which(annot_df$DOOLITTLE_ANNOT == 'LOW')) #for interface analysis I == I & S and NI == all other
    total_I <- length(which(annot_df$DOOLITTLE_ANNOT == 'HIGH'))
    top_I <- length(which(top_fst_df$DOOLITTLE_ANNOT == 'HIGH'))
    top_NI <- length(which(top_fst_df$DOOLITTLE_ANNOT == 'LOW'))
    successes <- c(top_I, top_NI)
    trials <- c(total_I, total_NI)
    prop_test_p_val <- prop.test(successes, trials, alternative = "greater")$p.value
    result_df <- data.frame('Top_High_Doolittle' = top_I, 'Top_Low_Doolittle' = top_NI, 'Doolittle_p_val' = prop_test_p_val, 'Virus' = vir)
    return(result_df)
  } else if (annot_col_name == 'surface') {
    total_NI <- length(which(annot_df$SURFACE_ANNOT == 'NS')) #for interface analysis I == I & S and NI == all other
    total_I <- length(which(annot_df$SURFACE_ANNOT == 'S'))
    top_I <- length(which(top_fst_df$SASA_ANNOT == 'NS'))
    top_NI <- length(which(top_fst_df$SASA_ANNOT == 'NS'))
    successes <- c(top_I, top_NI)
    trials <- c(total_I, total_NI)
    prop_test_p_val <- prop.test(successes, trials, alternative = "greater")$p.value
    result_df <- data.frame('Top_S' = top_I, 'Top_NS' = top_NI, 'Surface_p_val' = prop_test_p_val, 'Virus' = vir)
    return(result_df)
  } 
}
#function to get odds.ratio