#COVID TwoBP - TA RNAseq
#October 11 2023 

#####Section 1 - makes figures 2, 3, and 4 
library(limma)
library(dplyr)
library(tidyverse)
library(vegan)
##Inputs:  
###Sample data post bacterial quality control from trach aspirates or nasal swabs 
###Timepoint samples arel a subset of total samples, with one sample selected taken at the closest point to 2BP for 2BP patients, 
#and one sample from each control patient at a similar time on ventilator 
TA_samples<-read.csv("TA_samples_postbacterialQC_longitudinal.csv")  
TA_samples_timepoint <- read.csv("TA_samples_postbacterialANDhostQC_timepoint.csv")
NS_samples<-read.csv("NS_samples_postbacterialQC_longitudinal.csv")  
NS_samples_timepoint<-read.csv("NS_samples_postbacterialQC_timepoint.csv")
####Patient metadata and culture results 
metadata<-read.csv("Patient_data_noPHI.csv")
culture_results <- read.csv("2BP_CultureResults.csv")
###Microbial results - read outputs from CZID, SARS-CoV-2 nt_RPM, and bacterial mass calculated
genus_microbe_reports <- read.csv("Microbe_reports_genus.csv")  
sarsrpm<-read.csv("BetaCoV_nt_rpM.csv")
bacterial_mass_filt <- read.csv("Analyzed_samples.csv")

### PART 1: Generating FIGURE 3 (descriptive) depicting rank at trach aspirates over time  
## Read in sample and patient metadata 
TA_TwoBP_samples <- TA_samples %>%
  filter(Final == "2BP")

#Read in all genus level data and filter out contaminants
contaminants <- c("Sphingomonas", "Bradyrhizobium", "Ralstonia", "Delftia", 
                  "Cutibacterium", "Methylobacterium", "Acidovorax", "Chryseobacterium", "Burkholderia")

genus_microbe_reports <- genus_microbe_reports %>%
  filter(nt_rpm > 0.1 & nr_rpm > 0.1) %>%
  filter(!name %in% contaminants)

#Filter genus-level data for TA samples from culture+ 2BP patients 
genus_TA_TwoBP <- genus_microbe_reports %>%
  filter(sample_name %in% TA_TwoBP_samples$dl_id) %>% 
  filter(category == "bacteria") 

ID<-TA_TwoBP_samples[,c("comet_id", "dl_id")]
genus_TA_TwoBP_1<-merge(genus_TA_TwoBP, ID, by.x ="sample_name", by.y ="dl_id", all.x=TRUE, all.y=FALSE)        

#Add rank # to genus (by rpm)
genus_TA_TwoBP <- genus_TA_TwoBP_1 %>%
  group_by(sample_name) %>%
  arrange(desc(nt_rpm)) %>%
  mutate(rank = row_number())

#Merge on COMET ID and genus combo, so we're only keeping rows from microbe reports that are the cultured pathogen 
genus_TA_TwoBP_pathogen <- genus_TA_TwoBP %>%
  inner_join(culture_results[, c(1,5)], by = c("comet_id" = "patient_id", "name" = "PNAOrg_genus"))

#Merge with TA sample list and add zeros when there was no result 
genus_TA_TwoBP_pathogen <- genus_TA_TwoBP_pathogen %>%
  right_join(TA_TwoBP_samples[, c(1, 2, 3, 4, 5)], by = c("sample_name" = "dl_id", "comet_id" = "comet_id")) %>%
  tidyr::replace_na(list(nt_rpm=0)) %>%
  tidyr::replace_na(list(rank=200)) 

#Change all rank > 10 to be 12 for Figure 2
genus_TA_TwoBP_pathogen$rank_new<-genus_TA_TwoBP_pathogen$rank
for (i in (1:nrow(genus_TA_TwoBP_pathogen))){
  if (genus_TA_TwoBP_pathogen$rank[i] > 10)
  {genus_TA_TwoBP_pathogen$rank_new[i] = 12}} 

##To visualize one pathogen/patient combination in figure 2, divide each patient by pathogens
###for 1145, want 3 = Klebsiella, Escherichia, and Staphylococcus
####for 1250, want 2 = Klebsiella and Pseudomonas
## for 1474 want 2 -staphylococcus and haemophilus
for (i in (1:nrow(genus_TA_TwoBP_pathogen))){
  if (genus_TA_TwoBP_pathogen$comet_id[i] =="1250" & genus_TA_TwoBP_pathogen$name[i] =="Pseudomonas"  )
  {genus_TA_TwoBP_pathogen$comet_id[i] = "1250_1"}} 

 for (i in (1:nrow(genus_TA_TwoBP_pathogen))){
  if (genus_TA_TwoBP_pathogen$comet_id[i] =="1250" & genus_TA_TwoBP_pathogen$name[i] =="Klebsiella"  )
  {genus_TA_TwoBP_pathogen$comet_id[i] = "1250_2"}} 

for (i in (1:nrow(genus_TA_TwoBP_pathogen))){
  if (genus_TA_TwoBP_pathogen$comet_id[i] =="1145" & genus_TA_TwoBP_pathogen$name[i] =="Escherichia"  )
  {genus_TA_TwoBP_pathogen$comet_id[i] = "1145_1"}} 

for (i in (1:nrow(genus_TA_TwoBP_pathogen))){
  if (genus_TA_TwoBP_pathogen$comet_id[i] =="1145" & genus_TA_TwoBP_pathogen$name[i] =="Klebsiella"  )
  {genus_TA_TwoBP_pathogen$comet_id[i] = "1145_2"}} 

for (i in (1:nrow(genus_TA_TwoBP_pathogen))){
  if (genus_TA_TwoBP_pathogen$comet_id[i] =="1145" & genus_TA_TwoBP_pathogen$name[i] =="Staphylococcus"  )
  {genus_TA_TwoBP_pathogen$comet_id[i] = "1145_3"}} 

for (i in (1:nrow(genus_TA_TwoBP_pathogen))){
  if (genus_TA_TwoBP_pathogen$comet_id[i] =="1474" & genus_TA_TwoBP_pathogen$name[i] =="Staphylococcus"  )
  {genus_TA_TwoBP_pathogen$comet_id[i] = "1474_1"}} 

for (i in (1:nrow(genus_TA_TwoBP_pathogen))){
  if (genus_TA_TwoBP_pathogen$comet_id[i] =="1474" & genus_TA_TwoBP_pathogen$name[i] =="Haemophilus"  )
  {genus_TA_TwoBP_pathogen$comet_id[i] = "1474_2"}} 

##FIgure 3###
#Plot genus rank over time in each patient - note that here we consider ranks >10 as >12 for visualization 
print(ggplot(genus_TA_TwoBP_pathogen, aes(x = sample_to_2BP, y = rank_new)) + geom_point(color = "#D41159") +
        geom_vline(xintercept=0, color = "black") +
        theme_bw() + facet_wrap(~ comet_id, nrow = 4) + labs(x="Days from 2BP", y="Bacterial genus rank") + 
        scale_y_continuous(trans = "reverse", breaks = c(11, 9, 7, 5, 3, 1), limits = c(12, 1)) + 
        theme(aspect.ratio = 1) + geom_hline(yintercept=10.5, color = "grey", linetype = "dashed") +
        theme(panel.grid.minor = element_blank(), panel.background = element_blank()) +
        scale_x_continuous(breaks = c(-6, -4, -2, 0, 2, 4, 6), limits = c(-7, 7)))


###### Generation of data that makes up Figure 2 = Microbial analysis of TA samples
#Bacterial mass analyses - 2A
mass_timepoint<-merge(TA_samples_timepoint, bacterial_mass_filt, by.x="dl_id", by.y="sample_name", all.x=TRUE, all.y=FALSE)
ggplot(mass_timepoint, aes(x=Final, y = log(mass))) + geom_boxplot() + geom_jitter(width = 0.25) + theme_bw()     
mass_timepoint$mass<-as.numeric(mass_timepoint$mass)
wilcox.test(as.numeric(mass) ~ Final, data = mass_timepoint)
write.csv(mass_timepoint, "TA_samples_timepoint_bacterial_mass.csv")

####Alpha diversity analyses - figures 2B-D
all_reports <- genus_microbe_reports %>%
  filter(sample_name %in% TA_samples$dl_id) %>% 
  filter(category == "bacteria") 

ID<-TA_samples[,c("comet_id", "dl_id")]
genus_TA_1<-merge(all_reports, ID, by.x ="sample_name", by.y ="dl_id", all.x=TRUE, all.y=FALSE)        

table <- genus_TA_1 %>%
  pivot_wider(values_from = nt_rpm, names_from = "sample_name", id_cols = "name", values_fill = 0) %>%
  column_to_rownames("name") %>% 
  as.matrix %>%
  t(.)

#Diversity calculations - calculate alpha diversity using package Vegan
alpha_diversity <- diversity(table, index = "shannon")
#Join alpha diversity to samples 
alpha_diversity <- as.data.frame(alpha_diversity) %>% rownames_to_column("sample_name")

TA_samples_alpha <- alpha_diversity %>%
  left_join(TA_samples, by = c("sample_name" = "dl_id"))

#All samples are part of longitudinal samples
TA_samples_alpha_longitudinal<- TA_samples_alpha 

#Pulling alpha diversity for timepoint analyses ##FIgure 2B
TA_samples_alpha_timepoint<- TA_samples_alpha %>% filter(sample_name %in% TA_samples_timepoint$dl_id)
ggplot(TA_samples_alpha_timepoint, aes(x=Final, y = alpha_diversity)) + geom_boxplot() + geom_jitter(width = 0.25) + theme_bw()     


#Write the alpha values into CSV files
write.csv(TA_samples_alpha_timepoint, "TA_samples_alpha_timepoint.csv")
write.csv(TA_samples_alpha_longitudinal, "TA_samples_alpha_longitudinal.csv")

###For Figure 2C - comparing alpha diversity between samples with cultured pathogen at highest rank vs controls
###Select the samples with the cultured pathogen at highest rank - note, that this INCLUDES
#the samples with multiple cultured pathogens as multiple pathogens 
best_TA_path <- genus_TA_TwoBP_pathogen[, c(1, 2, 5, 12, 13, 36, 37, 38, 39, 40, 41, 42)] %>%
  filter(sample_to_2BP <= 2 & sample_to_2BP >= -3) %>%
  group_by(comet_id) %>%
  dplyr::slice(which.min(rank)) 

#Remove the multiple pathogens from patients with multiple pathogens
best_TA_path_1<- best_TA_path %>% filter(comet_id != "1145_2")
best_TA_path_1<- best_TA_path_1 %>% filter(comet_id != "1145_3")
best_TA_path_1<- best_TA_path_1 %>% filter(comet_id != "1250_2")
best_TA_path_1<- best_TA_path_1 %>% filter(comet_id != "1474_2")          

##Add the control samples 
timepointcontrols<-TA_samples_timepoint %>% filter(Final == "No-BP")
best_TA_path_1$Final<- NA
#best_TA_path_1 <-best_TA_path_1 %>% rename( "dl_id"= "sample_name" )
best_TA_path_1$dl_id<-best_TA_path_1$sample_name
TA_samples_toprank_timepoint<-rbind(timepointcontrols[, c(1, 2,3,4,6)], best_TA_path_1[,c(7,13,9,10,11)])
TA_samples_alpha_topranked_timepoint<- merge(TA_samples_alpha, TA_samples_toprank_timepoint, by.x="sample_name", by.y="dl_id", all.x=FALSE, all.y=TRUE)
wilcox.test(alpha_diversity ~ Final.x, data = TA_samples_alpha_topranked_timepoint)
ggplot(TA_samples_alpha_topranked_timepoint, aes(x=Final.x, y = alpha_diversity)) + geom_boxplot() + geom_jitter(width = 0.25) + theme_bw()     
write.csv(TA_samples_alpha_topranked_timepoint, "TA_alpha_vs_rank_of_top_pathogen_at_toprankedsample.csv")

#mann-whitney test: TwoBP vs No TwoBP with timepoint samples 
wilcox.test(alpha_diversity ~ Final, data = TA_samples_alpha_timepoint)

#For Figure 3D, create a graph of alpha diversity over time wtih a fitted line - this is #3D
ggplot(TA_samples_alpha_longitudinal, aes(x = vent_to_sample, y = alpha_diversity, color = Final)) + geom_point() + geom_smooth(method = lm) +
  theme_bw() + labs(x="Days from Intubation", y="Shannon Diversity Index") + theme(aspect.ratio =1) + theme_bw() +scale_color_manual(values=c("#1A85FF","#D41159"))

#Figure 3D perform p values calculated using linear mixed-effects model and and likelihood ratio test package to look at effect of 2BP vs no-BP 
fit2 <- lme4::lmer(alpha_diversity ~ vent_to_sample + Final + (1|comet_id), data=TA_samples_alpha, REML=FALSE)
fit3 <- lme4::lmer(alpha_diversity ~ vent_to_sample + (1|comet_id), data=TA_samples_alpha, REML=FALSE)
anova(fit3, fit2) # P value of the co-efficent term, Final (2BP/No BP)

#match with their rank - genus_TA_TwoBP_pathogen has the rank of the cultured pathogen in each patient at each timepoint 
genus_TA_TwoBP_pathogen_timepoint <- genus_TA_TwoBP_pathogen %>% filter(sample_name %in% TA_samples_timepoint$dl_id)

#For Figure 3F - analyzing of how many patients, is there a topranked cultured pathogen in 7 days before and after?
#for graphing, show rank >10 as 10
genus_TA_TwoBP_pathogen <- genus_TA_TwoBP_pathogen %>%
  mutate(rank_new = case_when(
    rank > 10 ~ 10,
    TRUE ~ rank
  ))

#remove the multiple pathogens per patient
genus_TA_TwoBP_pathogen_onepathperpt<- genus_TA_TwoBP_pathogen %>% filter(comet_id != "1145_2")
genus_TA_TwoBP_pathogen_onepathperpt<- genus_TA_TwoBP_pathogen_onepathperpt %>% filter(comet_id != "1145_3")
genus_TA_TwoBP_pathogen_onepathperpt<- genus_TA_TwoBP_pathogen_onepathperpt %>% filter(comet_id != "1250_2")
genus_TA_TwoBP_pathogen_onepathperpt<- genus_TA_TwoBP_pathogen_onepathperpt %>% filter(comet_id != "1474_2")

##this table asks what the rank is of the cultured pathogen across all samples, with 10 indicating 10 or greater
table(genus_TA_TwoBP_pathogen_onepathperpt$rank_new)

#This table asks where the highest rank of the cultured pathogen is across all samples per each patient
best_TA_sample <- genus_TA_TwoBP_pathogen_onepathperpt %>%
  group_by(comet_id) %>%
  dplyr::slice(which.min(rank)) 

table(best_TA_sample$rank)

###FIGURE 3G: Look at top ranked bacteria by mNGS at time of TwoBP vs No-BP 
###cut out some of the unnecessary cols 
TA_samples_timepoint_short<-TA_samples_timepoint[,c("dl_id", "comet_id", "Final")]

#Filter genus-level data for TA samples from culture+ TwoBP patients 
genus_TA_all <- genus_microbe_reports %>%
  filter(sample_name %in% TA_samples_timepoint_short$dl_id) %>% 
  filter(category == "bacteria") %>%
  left_join(TA_samples_timepoint_short, by = c("sample_name" = "dl_id"))

#Add rank # to genus (by rpm)
genus_TA_all_rank <- genus_TA_all %>%
  group_by(sample_name) %>%
  arrange(desc(nt_rpm)) %>%
  mutate(rank = row_number())

y<-genus_TA_all_rank[,c("sample_name", "name", "nt_rpm", "rank")]

#Pull out only the microbes ranked top one
topone<-y %>% filter(rank<2)

#join it with timepoint data - this becomes 3G
#Manual adjudication of S. epi and S. aureus performed 
topone_timepoint<-merge(TA_samples_timepoint_short,topone, by.x ="dl_id", by.y="sample_name")
table(topone_timepoint$name, topone_timepoint$Final)
write.csv(topone_timepoint, "toppathogenrankedat_TA_timepoint_sample.csv") 

###SarsCOV2 RPM at timepoint samples - supplementary table
sars<-merge (TA_samples_timepoint, sarsrpm, by.x="dl_id", by.y="sample_name", all.x=TRUE, all.y=FALSE)
ggplot(sars, aes(x=Final, y = log(nt_rpm))) + geom_boxplot() + geom_jitter(width = 0.25) + theme_bw()     
wilcox.test(as.numeric(nt_rpm) ~ Final, data = sars)
write.csv(sars, "TS_samples_timepoint_sarscov2_ntrpm.csv")


### NASAL SWABS ### - FIGURE 4


###NS samples mass - FIgure 4A
mass_NS_timepoint<-merge(NS_samples_timepoint, bacterial_mass_filt, by.x="dl_id", by.y="sample_name", all.x=TRUE, all.y=FALSE)
ggplot(mass_NS_timepoint, aes(x=Final, y = log(mass))) + geom_boxplot() + geom_jitter(width = 0.25) + theme_bw()     
wilcox.test(mass ~ Final, data = mass_NS_timepoint)
write.csv(mass_NS_timepoint, "NS_samples_mass_timepoint.csv")

# Alpha diversity at time closest to TwoBP - this makes figure 4B
###Calculating NS Alpha diversity
# NS samples alpha diversity 
#NS_samples_timepoint <- read.csv("NS_samples_postbacterialQC_timepoint.csv")
NS_all_timepoints <- genus_microbe_reports %>%
  filter(sample_name %in% NS_samples_timepoint$dl_id) %>%
  filter(category == "bacteria") 

NS_table <- NS_all_timepoints %>%
  pivot_wider(values_from = nt_rpm, names_from = "sample_name", id_cols = "name", values_fill = 0) %>%
  column_to_rownames("name") %>% 
  as.matrix %>%
  t(.)
#Diversity calculations 
NS_alpha_diversity <- diversity(NS_table, index = "shannon")
NS_alpha_diversity <- as.data.frame(NS_alpha_diversity) %>% rownames_to_column("sample_name")

NS_samples_alpha <- NS_alpha_diversity %>%
  left_join(NS_samples_timepoint, by = c("sample_name" = "dl_id"))
#NS_samples_alpha <- read.csv()
ggplot(NS_samples_alpha, aes(x=Final, y = NS_alpha_diversity)) + geom_boxplot() + geom_jitter(width = 0.25) + theme_bw()     
#mann-whitney test: TwoBP vs No TwoBP at closest timepoint to TwoBP 
wilcox.test(NS_alpha_diversity ~ Final, data = NS_samples_alpha)
write.csv(NS_samples_alpha, "NS_samples_alpha_timepoint.csv")

##For FIgure 4D, create a timecourse file of this with rank of culture-confirmed pathogen over time (this will be figure 4D)
#Filter available nasal swab samples for TwoBP patients from samples that past host/bacterial QC 
#NS_TwoBP_samples_filt <- read.csv("NS_samples_postbacterialQC_longitudinal.csv") %>%  # 
#  filter(Final == "TwoBP")
  NS_TwoBP_samples_filt <- NS_samples %>%  # 
    filter(Final == "2BP")
#Filter genus-level data for NS samples from culture+ TwoBP patients 
#Only keep samples that passed QC --> left with 15 TwoBP patients
genus_NS_TwoBP <- genus_microbe_reports %>%
  filter(sample_name %in% NS_TwoBP_samples_filt$dl_id) %>% 
  filter(category == "bacteria") %>%
  left_join(NS_TwoBP_samples_filt[, c(1, 2, 3,4, 5,6)], by = c("sample_name" = "dl_id"))

#Add rank # to genus (by rpm)
genus_NS_TwoBP <- genus_NS_TwoBP %>%
  group_by(sample_name) %>%
  arrange(desc(nt_rpm)) %>%
  mutate(rank = row_number())

#For each patient, filter for results for cultured pathogen
#Merge on COMET ID and pathogen name combo
genus_NS_TwoBP$comet_id<-as.integer(genus_NS_TwoBP$comet_id)
genus_NS_TwoBP_pathogen <- genus_NS_TwoBP %>%
  inner_join(culture_results[, c(1,5)], by = c("comet_id" = "patient_id", "name" = "PNAOrg_genus"))
table(genus_NS_TwoBP_pathogen$rank)        
#1233 drops out because there is no serratia detected in it at all, indicated as lower than limit of detection 



#Figure 4C
#Compare rank in TA vs NS samples
#Change all rank > 10 to be>10 - this is for visualization in all figures 
genus_NS_TwoBP_pathogen <- genus_NS_TwoBP_pathogen %>%
  mutate(rank_new = case_when(
    rank > 10 ~ 10,
    TRUE ~ rank
  ))

genus_NS_TwoBP_pathogen_fortable<- genus_NS_TwoBP_pathogen %>% filter(comet_id != "1474_2")
genus_NS_TwoBP_pathogen_fortable<- genus_NS_TwoBP_pathogen_fortable %>% filter(comet_id != "1145_2")
genus_NS_TwoBP_pathogen_fortable<- genus_NS_TwoBP_pathogen_fortable %>% filter(comet_id != "1145_3")
table(genus_NS_TwoBP_pathogen_fortable$rank_new)

#Match all samples their rank - genus_TA_TwoBP_pathogen has the rank of the cultured pathogen in each patient at each timepoint 
genus_NS_TwoBP_pathogen_timepoint <- genus_NS_TwoBP_pathogen_fortable %>% filter(sample_name %in% NS_samples_timepoint$dl_id)
table(genus_NS_TwoBP_pathogen_timepoint$rank_new)


##Figure out of how many patients have the topranked pathogen in the 7 day window - figure 4C
##For 1233, the pathogen was never detected at all (hence >10)
best_NS_path <- genus_NS_TwoBP_pathogen_fortable %>%
  group_by(comet_id) %>%
  dplyr::slice(which.min(rank)) 
table(best_NS_path$rank_new)



##For visualization in figure 4D
##goal here is to divide patietns in to different 'patients' for each bug 
###for 1145, want 3 = Klebsiella, Escherichia, and Staphylococcus
####for 1250, want 2 = Klebsiella and Pseudomonas
## for 1474 want 2 -sNSphylococcus and haemophilus

for (i in (1:nrow(genus_NS_TwoBP_pathogen))){
  if (genus_NS_TwoBP_pathogen$comet_id[i] =="1145" & genus_NS_TwoBP_pathogen$name[i] =="Escherichia"  )
  {genus_NS_TwoBP_pathogen$comet_id[i] = "1145_1"}} 

for (i in (1:nrow(genus_NS_TwoBP_pathogen))){
  if (genus_NS_TwoBP_pathogen$comet_id[i] =="1145" & genus_NS_TwoBP_pathogen$name[i] =="Klebsiella"  )
  {genus_NS_TwoBP_pathogen$comet_id[i] = "1145_2"}} 

for (i in (1:nrow(genus_NS_TwoBP_pathogen))){
  if (genus_NS_TwoBP_pathogen$comet_id[i] =="1145" & genus_NS_TwoBP_pathogen$name[i] =="Staphylococcus"  )
  {genus_NS_TwoBP_pathogen$comet_id[i] = "1145_3"}} 

for (i in (1:nrow(genus_NS_TwoBP_pathogen))){
  if (genus_NS_TwoBP_pathogen$comet_id[i] =="1474" & genus_NS_TwoBP_pathogen$name[i] =="Staphylococcus"  )
  {genus_NS_TwoBP_pathogen$comet_id[i] = "1474_1"}} 

for (i in (1:nrow(genus_NS_TwoBP_pathogen))){
  if (genus_NS_TwoBP_pathogen$comet_id[i] =="1474" & genus_NS_TwoBP_pathogen$name[i] =="Haemophilus"  )
  {genus_NS_TwoBP_pathogen$comet_id[i] = "1474_2"}} 

#Merge with NS sample list and add zeros when there was no result 
genus_NS_TwoBP_pathogen$comet_id<-as.character(genus_NS_TwoBP_pathogen$comet_id)
NS_TwoBP_samples_filt$comet_id<-as.character(NS_TwoBP_samples_filt$comet_id)


#Change all rank > 10 to be 10
genus_NS_TwoBP_pathogen <- genus_NS_TwoBP_pathogen %>%
  mutate(rank_new = case_when(
    rank > 10 ~ 10,
    TRUE ~ rank
  ))

# TA and NS combined pathogen plots 
#Combine all samples
genus_TA_TwoBP_pathogen$comet_id<-as.character(genus_TA_TwoBP_pathogen$comet_id)
genus_TA_TwoBP_pathogen$sample_type <- 1

genus_NS_TwoBP_pathogen$comet_id<-as.character(genus_NS_TwoBP_pathogen$comet_id)
genus_NS_TwoBP_pathogen$sample_type<-0
all_samples <- rbind(genus_TA_TwoBP_pathogen[, c(1, 5, 12, 37, 42, 39,41, 38)], genus_NS_TwoBP_pathogen[, c(2, 5, 12, 37, 43,38,40, 42)])

#Only include comet IDs that have NS data (13) OR 1233, which has a NS sample in which the bug is not detected
all_samples <- all_samples %>%
  filter(comet_id %in% genus_NS_TwoBP_pathogen$comet_id | comet_id == "1233")

#Change all rank > 10 to be 10
all_samples <- all_samples %>%
  mutate(rank_new = case_when(
    rank > 10 ~ 12,
    TRUE ~ rank
  ))

all_samples$sample_type<-as.character(all_samples$sample_type)
##Graphing both TA and NS samples 
ggplot(all_samples, aes(x = sample_to_2BP, y = rank_new, color = sample_type)) + geom_point() +
  geom_vline(xintercept=0, color = "black") +
  theme_bw() + facet_wrap(~ comet_id, nrow = 2) + labs(x="Days from TwoBP", y="Bacterial genus rank") + 
  scale_y_continuous(trans = "reverse", breaks = c(11, 9, 7, 5, 3, 1), limits = c(12, 1)) + theme(aspect.ratio = 1) +
  geom_hline(yintercept=10.5, color = "darkgrey", linetype = "dashed") + 
  scale_color_manual(values = c("mediumpurple", "#D41159")) +
  theme(panel.grid.minor = element_blank(), panel.background = element_blank()) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2, 4, 6), limits = c(-7, 7))
