### RBP-network analysis

#Manipulation of Bait-Target Master Table 

library(tidyverse)
library(dplyr)
library(ggplot2)
#install.packages("readxl")
library(readxl)

#read-in table 
master_table<- read_csv("/Users/katierothamel/Dropbox (Personal)/Rothamel,Katherine/Projects/RBP_Networks/RawData/20220417_RBP_fullData2_final_PPI_table_withIntensity_IP_DeepProteome_forStoichiometry.csv")

#No self-loppps
master_table<- master_table %>% filter(Prey != Bait) 

#ordering bait based of life-cyclestep
master_table <- master_table[order(master_table$Life_cycle_step),]

#write_csv(master_table, "20220613_noselfloop_master_table.csv")

#Plotting interaction type for each unique bait, faceted on RNA-life cycle

#Ordering Bargraph
master_table %>%  
ggplot(aes(x = fct_reorder(Bait, Life_cycle_step), fill = Interaction_type))+ 
  theme_minimal() +
  geom_bar(alpha=.5) +
  labs(y= "Number of Interactions", title = "Bait and Prey First Pass") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text=element_text(size=rel(0.6))) +
  theme(legend.position="bottom")

master_table %>%
  #filter(Interaction_type == "Direct") %>%
  ggplot( aes(x = Bait,fill =  PPI_support))+ 
  theme_minimal() +
  geom_bar(alpha=.5) +
  labs(y= "Number of Interactions", title = "Direct Bait and Prey Interactions") +
  facet_wrap(~ Life_cycle_step , scales = "free_x") +
  #facet_wrap(~PPI_support) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text=element_text(size=rel(0.6))) +
  theme(legend.position="bottom")

master_table %>%
  #filter(Interaction_type == "RNA mediated") %>%
  ggplot( aes(x = Bait, fill =  Interaction_type))+ 
  theme_minimal() +
  geom_bar(alpha=.5) +
  labs(y= "Number of Interactions", title = "RNA mediated Bait and Prey Interactions") +
  facet_wrap(~ Life_cycle_step , scales = "free_x") +
  #facet_wrap(~PPI_support) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text=element_text(size=rel(0.6))) +
  theme(legend.position="bottom")

#Is intensity score differ do to interaction type? We have +/- RNAse conditions. High Bait/Prey ratio means more bait compared to prey. 
library(ggrepel)

ggplot(master_table, aes(x= log(`bait_prey_intensity-ratio-IP_R`),y= log(`bait_prey_intensity-ratio-IP_noR`), color = Interaction_type, label= Prey))+ 
  geom_point() +
  theme_minimal() +
  geom_abline(intercept = 0, slope = 1) +
  theme(axis.text.x = element_text(angle = 90)) 

master_table %>% filter(Interaction_type == c("RNA shielded")) %>% 
ggplot(aes(x= log(`bait_prey_intensity-ratio-IP_R`),y= log(`bait_prey_intensity-ratio-IP_noR`), label= Prey))+ 
  theme_minimal() +
  geom_point(alpha=0.9,  color = "#619CFF") +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_text_repel(max.overlaps = 20)

###count the average number

meancounts<- master_table %>% 
  filter(Interaction_type =="RNA shielded") %>% 
  group_by(Bait) 

sum(meancounts$n)

length(which(meancounts$PPI_support == "False"))

  master_table %>% 
    gather(condition, intensity, bait_prey_intensity.ratio.IP_R:bait_prey_intensity.ratio.IP_noR) %>% 
  #  dplyr::filter(Prey %in% c("G3BP2")) %>% 
  ggplot( aes(x= Interaction_type, y= log(intensity), fill = Interaction_type)) +
    geom_boxplot(alpha=0.5) +
    theme_minimal()+
    #geom_hline(aes(yintercept = log(bait_intensity_noR)), linetype = "dashed", alpha= 0.5, color = "red") +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(legend.position="top") +
    facet_wrap(~condition)
  
  master_table %>% 
    gather(condition, intensity, bait_prey_intensity.ratio.IP_R:bait_prey_intensity.ratio.IP_noR) %>% 
    #  dplyr::filter(Prey %in% c("G3BP2")) %>% 
    ggplot( aes(x= Interaction_type, y= log(intensity), fill = Interaction_type)) +
    geom_boxplot(alpha=0.5) +
    theme_minimal()+
    #geom_hline(aes(yintercept = log(bait_intensity_noR)), linetype = "dashed", alpha= 0.5, color = "red") +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(legend.position="top") +
    facet_wrap(~condition) 
  
test<-   master_table %>%  
  filter(Bait %in% c("DLD", "DLST", "DLAT", "OGDH", "PDHA1", "PDHB", "PDHX")) 

# plotting no_R vs R
  master_table %>%  filter(Bait %in% c("DLD", "DLST", "DLAT", "OGDH", "PDHA1", "PDHB", "PDHX")) %>% 
    gather(condition, intensity, bait_prey_intensity.ratio.IP_R:bait_prey_intensity.ratio.IP_noR)%>% 
    ggplot(aes(x = Interaction_type, y= log(intensity), fill = Interaction_type)) +
    geom_boxplot(alpha=0.5) +
    theme_minimal()+
    facet_wrap(~condition) +
    #geom_hline(aes(yintercept = log(bait_intensity_noR)), linetype = "dashed", alpha= 0.5, color = "black") +
    #geom_hline(aes(yintercept = log(bait_intensity_R)), linetype = "dashed", alpha= 0.5, color = "gray") +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(legend.position="top") +
    labs(title= "PDH Interactors")
  
  master_table %>% 
    gather(condition, intensity, bait_prey_intensity.ratio.IP_R:bait_prey_intensity.ratio.IP_noR) %>% 
    filter(Bait %in% c("EXOSC10", "EXOSC2", "EXOSC9")) %>% 
    ggplot(aes(x = Interaction_type, y= log(intensity), fill = Interaction_type)) +
    geom_boxplot(alpha=0.5) +
    theme_minimal()+
    facet_wrap(~condition) +
    #geom_hline(aes(yintercept = log(bait_intensity_noR)), linetype = "dashed", alpha= 0.5, color = "black") +
    #geom_hline(aes(yintercept = log(bait_intensity_R)), linetype = "dashed", alpha= 0.5, color = "gray") +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(legend.position="top") +
    labs(title= "Exosome Interactors" )
  
  master_table %>%  
    gather(condition, intensity, bait_prey_intensity.ratio.IP_R:bait_prey_intensity.ratio.IP_noR) %>% 
   # filter(Bait %in% c("KARS1", "MARS1", "RARS1", "DARS")) %>% 
    ggplot(aes(x = Interaction_type, y= log(intensity), fill = Interaction_type)) +
    geom_violin(alpha=0.5) +
    geom_boxplot( width = .5, aplpha = 0.5) +
    theme_minimal()+
    facet_wrap(~condition) +
    #geom_hline(aes(yintercept = log(bait_intensity_noR)), linetype = "dashed", alpha= 0.5, color = "black") +
    #geom_hline(aes(yintercept = log(bait_intensity_R)), linetype = "dashed", alpha= 0.5, color = "gray") +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(legend.position="top") +
    labs(title= "Global Interactions")
  
  master_table %>% filter(Prey == "CPSF2") %>%
    ggplot( aes(x= log(bait_prey_ratio_R), y= log(bait_prey_ratio_noR), color = Interaction_type, label= Bait)) +
    geom_point(aes( size = 5)) +
    theme_minimal()+
    geom_text_repel()
    #geom_hline(yintercept = log(25000), linetype = "dashed", alpha= 0.5, color = "red") +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(legend.position="top") 
    
  
  library(ggrepel)
  
master_table %>% 
    ggplot(aes(x = log(bait_prey_intensity.ratio.IP_R), 
               y= log(bait_prey_intensity.ratio.IP_noR))) +
    geom_abline(slope = 1, intercept = 0) +
    geom_jitter(aes(color = Interaction_type)) +
    #geom_label_repel(aes(label = Prey)) +
    theme_minimal() 
    


master_table %>%
  gather(condition, intensity, bait_prey_intensity.ratio.IP_R:bait_prey_intensity.ratio.IP_noR) %>%
  ggplot(aes(x = Interaction_type, y= log(intensity), fill = Interaction_type))+
  #geom_hline(yintercept = 1) + geom_hline(yintercept = -1) +
  #geom_vline(xintercept = 1) + geom_vline(xintercept = -1) +
  #geom_abline(1, intercept = 0) +
  geom_boxplot(alpha= 0.5) +
  theme_minimal() +
  facet_wrap(~ Life_cycle_step + condition) 

master_table %>% 
  gather(condition, intensity, bait_prey_intensity.ratio.IP_R:bait_prey_intensity.ratio.IP_noR)  %>%
  ggplot(aes(x = Interaction_type, y= log(intensity), color = Interaction_type))+
  #geom_hline(yintercept = 1) + geom_hline(yintercept = -1) +
  #geom_vline(xintercept = 1) + geom_vline(xintercept = -1) +
  geom_abline(1, intercept = 0) +
  geom_point() +
  theme_minimal() 

  
#wrangling data into matrix (0,1) for interactions with Bait also separating by interaction
clean_long_master_table<- master_table %>% 
  filter(Interaction_type %in% c("RNA mediated", "Direct")) %>%
  select(Prey,Bait, bait_prey_intensity.ratio.IP_noR)%>% 
  unique() %>%
  group_by(Prey) %>% 
  tidyr::spread(Bait,bait_prey_intensity.ratio.IP_noR)

clean_long_master_table[is.na(clean_long_master_table)] = 0
clean_long_master_table<-clean_long_master_table %>% 
  mutate_if(is.numeric, ~1 * (. > 0)) 

corRBP<-cor(clean_long_master_table[,-1:-2])
corRBP[is.na(corRBP)] = 0
write.csv(corRBP, "corRBP_DirectRNAmediated_BaitBait")

length(unique(master_table$Bait))

RNAmed_long_master_table<- master_table %>% 
  filter(Interaction_type == "RNA mediated") %>%
  select(Prey, Bait, bait_prey_ratio_noR)%>% 
  group_by(Bait) %>% unique() %>%
  spread(Prey, bait_prey_ratio_noR)

RNAmed_long_master_table[is.na(RNAmed_long_master_table)] = 0
RNAmed_long_master_table<-RNAmed_long_master_table %>% 
  mutate_if(is.numeric, ~1 * (. > 0)) 

#annot_clean_long<- left_join(clean_long_master_table, master_table, by="Bait")


# heatmap and hierarchical clustering   
library(pheatmap)
library(RColorBrewer)

direct<-cor(direct_long_master_table[,-1])
RNAmed<-cor(RNAmed_long_master_table[,-1])
RNAmed[is.na(RNAmed)] = 0

direct_cor_matrix<-pheatmap(RNAmed, 
                          color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(1000), 
                          cellwidth = 1, 
                          cellheight = 1, 
                          cluster_cols=T, 
                          cluster_rows = T,  
                          show_rownames = F,
                          show_colnames = F)

RNAmed_cor_matrix<-pheatmap(RNAmed, 
                          color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(1000), 
                          cellwidth = 1, 
                          cellheight = 1, 
                          cluster_cols=T, 
                          cluster_rows = T,  
                          show_rownames = F,
                          show_colnames = F, cutree_rows = 12
                            )




(direct_cor_matrix$gtable$childrenOrder)



direct_hclust_gene <- hclust(dist(RNAmed), method = "complete")

direct_hclust_gene_ordered= data.frame(direct_hclust_gene$order)
direct_hclust_gene_ordered_bind = cbind(RNAmed, direct_hclust_gene_ordered[,1])

write.csv(direct_hclust_gene_ordered_bind, "test")


##################3
go_ontology<- read_excel("/Users/katierothamel/Desktop/Go_Ontology_Preys.xlsx")
go_ontology_long<- go_ontology %>% 
  gather(Column, Prey, -'High level GO category', -N) %>% na.omit()

write.csv(go_ontology_long, "go_ontology_long")

####
master_GO_table<-  left_join(master_table, go_ontology_long, by= "Prey")

ThreePrime_processing<-master_GO_table %>% filter(Life_cycle_step == "3' end processing")
write_csv(ThreePrime_processing, "ThreePrime_Processing")


ggplot(ThreePrime_processing, aes(x = Bait, fill = `High level GO category`))+ 
  theme_minimal() +
  geom_bar(alpha=.5) +
  labs(y= "Number of Interactions", title = "Bait and Prey First Pass") +
  #facet_wrap(~ `High level GO category` ) +
  theme(axis.text.x = element_text(angle = 90))


###############

#filtering master table by k means
kmean_queries<- read_csv("/Users/katierothamel/Dropbox (Personal)/Rothamel,Katherine/Projects/RBP_Networks/12kmeans_clusterz_preyinteraction.csv")

cluster1<- kmean_queries %>% 
  filter(k_means_12 == "12") %>% 
  rename(Bait = id )

test<- left_join(cluster1, master_table, by =  "Bait") %>% 
  filter(Prey != "NA")

write_csv(test, "cluster_12.1_mastertable")


head(master_table)
head(clean_long_master_table)

############### RBP disease association
RBPdiseases<-read_excel("/Users/katierothamel/Dropbox (VU Basic Sciences)/Ascano Lab/Manuscripts/DynamicPTGR_Review/NRG-20-113 Hentze Supplementary Table 2.xlsx")

RBPdiseases_reduced<- RBPdiseases %>% 
  rename(Prey = gene_name) %>%
  group_by(Prey,disease_name) %>% 
  summarise(Exp = unique(Therapeutic_area_name))


t_master_table<-read_excel("/Users/katierothamel/Dropbox (Personal)/Rothamel,Katherine/Projects/RBP_Networks/Trans_BaitPreyInteraction.xlsx")


disease_anno_KB_table<- full_join(master_table, RBPdiseases_reduced, by="Prey") %>% na.omit()

write_csv(disease_anno_KB_table, "disease_anno_KB_table")

ggplot(disease_anno_KB_table, aes(x=Therapeutic_area_name)) +
 geom_bar() + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text=element_text(size=rel(0.6)))

test<-disease_anno_KB_table %>% 
  group_by(Bait, Prey) %>% 
  summarise_()

reactome_table<- read.csv("/Users/katierothamel/Desktop/ReactomeTest.csv") 
reactome_table2<- reactome_table %>% 
  gather(Column, Gene, -'Pathway.name') 

write.csv(reactome_table2, "reactome_table")

clusters<- read.csv("/Users/katierothamel/Dropbox (Personal)/Rothamel,Katherine/Projects/RBP_Networks/ClustergeneList_ForPreyPreyHeatMap.csv")


long_clusters<- clusters %>% gather(Cluster, Prey, na.rm = TRUE) 

master_cluster<- full_join(long_clusters, master_table, by= "Prey")
master_cluster_test<- master_cluster[is.na(master_cluster$Bait) == FALSE , ]


write.csv(pot_bridges, "Prey_potential_bridges.csv")


