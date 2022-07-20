### RBP-network analysis
library(tidyverse)
library(dplyr)
library(ggplot2)
library(readxl)

#Wenhao's Output from 062221
master_table<- read.csv("/Users/katierothamel/Dropbox (Personal)/Rothamel,Katherine/Projects/RBP_Networks/20210223_RBP_Mix1_87_goodReps_borderline_scatter_volcano_overlap_edgelist_PPIannotation_full_LifeCycle_addECLIPinfo_sharedCLIPTargets_HumanProteinAtlasLocalization_addHydRaInfo_PSPscore_Disgenet_CellMapLocalizationNMF_SAFE_multiLog2FC_20210.csv")

#Removing self-loops 
master_table<- master_table %>% 
  filter(Prey != Bait)

#PARCLIP data exist from Mukherjee, 2018
PAR_RBPs<- read_excel("/Users/katierothamel/Dropbox (Personal)/Rothamel,Katherine/Citations/gky1185_supplemental_files/SupplementalTable1.xlsx")

#Add column that it has PAR-CLIP data
master_table<- master_table %>%
  mutate(PARCLIP = master_table$Prey %in% PAR_RBPs$RBP)

# C.GO et al, 2021 Biotinylation Interactions
biotin_local<- read_excel("/Users/katierothamel/Dropbox (Personal)/Rothamel,Katherine/Projects/RBP_Networks/Biotin_Localization.xlsx") %>% 
  rename(Prey= gene)

master_table_<- left_join(master_table, biotin_local, by = "Prey")

#David Ontology including BP,MF,CC,Pfam, KEGG_pathway 
anno_table<-read_excel("/Users/katierothamel/Dropbox (Personal)/Rothamel,Katherine/Projects/RBP_Networks/Old_Annotations/David_AnnotationofoPrey.xlsx") %>% 
  rename(Prey = ID)

combined_data<- left_join(master_table_, anno_table, by= "Prey")

prey_annotation <- combined_data %>% 
  select(Prey,prey_eclipped,Life_cycle_step, Prey_main_location, Prey_is_RBP, Prey_Disgenet_disease_.1, Prey_Disgenet_disease_.2, PARCLIP, GOTERM_BP_DIRECT, GOTERM_CC_DIRECT, GOTERM_MF_DIRECT, OMIM_DISEASE, PFAM, KEGG_PATHWAY, OMIM_DISEASE, Prey_SAFE_localization, Prey_MMF_localization)

write_csv(prey_annotation, "PreyAnnoTable_071321")

#Plotting interaction type for each unique bait, faceted on RNA-life cycle
combined_data$x <- factor(combined_data$x,                                    # Change ordering manually
                  levels = c("CPSF7","CSTF2","CSTF2T","PCF11","SYMPK","CASC3","CNOT2","DCP2","DGCR8","EXOSC10","EXOSC2","EXOSC9","LIN28B","PAN2","PARN","PUM2","TIA1","UPF2","XRN1","ALYREF","CHTOP","FYTTD1","NUP35","THOC1","XPO5","Bait","IGF2BP1","IGF2BP2","METTL3","BUD13","CACTIN","DHX38","GPKOW","HNRNPC","HNRNPL","LSM7","MFAP1","PRPF3","PRPF8","PTBP1","RBFOX2","RBM22","RBM8A","SART1","SF1","SF3A2","SF3B4","SNRNP200","SNRPA1","SNRPC","SNRPD3","SNU13","SRSF1","SRSF4","SRSF5","U2AF2","WBP11","LARP7","POLR2A","POLR2A_pSer2","POLR2A_pSer5","POLR2A_pSer7","TAF15","DDX3X","EIF3A","EIF4E","EIF4G2","MSI2","RPL18","RPL23A","RPLP0","RPS2","RPS24","RPS3","RPS5","UBAP2L","EWSR1"))

ggplot(combined_data, aes(x = Bait, fill = Interaction_type))+ 
  theme_minimal() +
  geom_bar(alpha=.5) +
  labs(y= "Number of Interactions", title = "Bait and Prey First Pass") +
  facet_wrap(~ Life_cycle_step, scales = "free_x") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text=element_text(size=rel(0.6))) +
  theme(legend.position="bottom")

 ggplot(combined_data, aes(x = Bait, fill = PPI_support))+ 
   theme_minimal() +
   geom_bar(alpha=.5) +
   labs(y= "Number of Interactions", title = "Bait and Prey First Pass") +
   facet_wrap(~ Life_cycle_step, scales = "free_x") +
   theme_minimal()+
   theme(axis.text.x = element_text(angle = 90)) +
   theme(axis.text=element_text(size=rel(0.6))) +
   theme(legend.position="bottom")
 
 
 #wrangling data into matrix (0,1) for interactions with Bait also seperating by interaction
 clean_long_master_table<- combined_data %>% 
   filter(Interaction_type == "Direct") %>%
   select(Prey, Bait)%>% 
   group_by(Bait) %>%
   mutate(Bait_ID = row_number()) %>%
   spread(Prey, Bait_ID)
 
 clean_long_master_table[is.na(clean_long_master_table)] = 0
 clean_long_master_table<-clean_long_master_table %>% 
   mutate_if(is.numeric, ~1 * (. > 0)) 
 
 #long data Prey-Prey Analysis log2FC
 direct_noR_clean_long_master_table<- combined_data %>%
   filter(Interaction_type != "RNA shielded") %>%
   select(Prey, Bait,log2FC.noR_for_direct.  ) %>% 
   group_by(Bait) %>%
   spread(Prey, log2FC.noR_for_direct.)
 
direct_noR_clean_long_master_table[is.na(direct_noR_clean_long_master_table)] = 0
 
direct_noR_log2corRBP<-cor(direct_noR_clean_long_master_table[,-1], method= "pearson")
 
write.csv(direct_noR_log2corRBP, "RNAmediated_direct_noR_log2corRBP_PreyPrey")
 

########################
direct_R_clean_long_master_table<- combined_data %>%
  filter(Interaction_type != "RNA shielded") %>%
   select(Prey, Bait,log2FC.R_for_direct.) %>% 
   group_by(Bait)%>%
   spread(Prey, log2FC.R_for_direct.)
 
direct_R_clean_long_master_table[is.na(direct_R_clean_long_master_table)] = 0
 
direct_R_log2corRBP<-cor(direct_R_clean_long_master_table[,-1], method= "pearson")
 
write.csv(direct_R_log2corRBP, "directRNAmediated_R_log2corRBP_PreyPrey")
 
 
#clean_long_master_table<-clean_long_master_table %>% 
  # mutate_if(is.numeric, ~1 * (. > 0)) 
 
log2corRBP<-cor(clean_long_master_table[,-1], method= "pearson")
 
write.csv(log2corRBP, "ALL_log2FC_preyprey_corRBP")

hist(corRBP)

#Bait-Bait Analysis 
clean_long_master_table_bait<- combined_data %>%
  select(Prey, Bait)%>% 
  group_by(Bait) %>%
  mutate(Bait_ID = row_number()) %>%
  spread(Bait, Bait_ID)

clean_long_master_table_bait[is.na(clean_long_master_table_bait)] = 0
clean_long_master_table_bait<-clean_long_master_table_bait %>% 
  mutate_if(is.numeric, ~1 * (. > 0)) 

corRBP_Bait<-cor(clean_long_master_table_bait[,-1])

write.csv(corRBP_Bait, "Bait_Bait_All_corRBP")


 



