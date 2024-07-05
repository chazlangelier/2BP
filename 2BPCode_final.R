#COVID TwoBP - TA RNAseq
#June 24 2024

#####This code makes figures 2A-G, 3, and 4A-C 
library(limma)
library(dplyr)
library(tidyverse)
library(vegan)
library(ggplot2)
##Inputs:   

###Sample data post bacterial quality control from trach aspirates or nasal swabs 
###Timepoint samples are a subset of total samples, with one sample selected taken at the closest point to 2BP for 2BP patients, 
#and one sample from each control patient at a similar time on ventilator 

TA_samples<-read.csv("Inputs_Final/TA_samples_longitudinal.csv")  
TA_samples_timepoint <- read.csv("Inputs_Final/TA_samples_timepoint.csv")
NS_samples<-read.csv("Inputs_Final/NS_samples_longitudinal.csv")  
NS_samples_timepoint<-read.csv("Inputs_Final/NS_samples_timepoint.csv")

####Patient metadata and culture results 
metadata<-read.csv("Inputs_Final/Patient_data_noPHI.csv")
culture_results <- read.csv("Inputs_Final/2BP_CultureResults.csv")

###Microbial results - read outputs from CZID  and bacterial mass
genus_microbe_reports <- read.csv("Inputs_Final/Microbe_reports_genus.csv")  
bacterial_mass_filt <- read.csv("Inputs_Final/2BP_Mass_Read_Output.csv")

###Establish library means and ranges
mean(bacterial_mass_filt$total_reads)/10^7
range(bacterial_mass_filt$total_reads)/10^7
mean(bacterial_mass_filt$bacterial_reads) /10^6
range(bacterial_mass_filt$bacterial_reads)/10^6

### PART 1: Generating FIGURE 3 (descriptive) and all descriptive supplementary files depicting rank, mass, or normalized mass of cultured pathogens
#in at trach or nasal swab aspirates over time  
## Read in sample and patient metadata 
TA_TwoBP_samples <- TA_samples %>%
  filter(TwoBP_Status == "2BP")
NS_TwoBP_samples <- NS_samples %>%
  filter(TwoBP_Status == "2BP")
#Read in all genus level data and filter out contaminants
contaminants <- c("Sphingomonas", "Bradyrhizobium", "Ralstonia", "Delftia", 
                  "Cutibacterium", "Methylobacterium", "Acidovorax", "Chryseobacterium", "Burkholderia")

genus_microbe_reports <- genus_microbe_reports %>%
  filter(nt_rpm > 0.1 & nr_rpm > 0.1) %>%
  filter(!name %in% contaminants)

#Filter genus-level data for TA and NS samples from culture+ 2BP patients 
genus_TwoBP <- genus_microbe_reports %>%
  filter(sample_name %in% TA_TwoBP_samples$sample_name|sample_name %in% NS_TwoBP_samples$sample_name) %>% 
  filter(category == "bacteria") 
ID_TA<-TA_TwoBP_samples[,c("manuscript_id", "sample_name")]
ID_NS<-NS_TwoBP_samples[,c("manuscript_id", "sample_name")]
ID<-rbind(ID_TA, ID_NS)
genus_TwoBP_1<-merge(genus_TwoBP, ID, by.x ="sample_name", by.y ="sample_name", all.x=TRUE, all.y=FALSE)        
####divide into TA and NS samples 
genus_TA_TwoBP<-genus_TwoBP_1%>% filter(sample_name %in% TA_TwoBP_samples$sample_name)
genus_NS_TwoBP<-genus_TwoBP_1%>% filter(sample_name %in% NS_TwoBP_samples$sample_name)

#Add rank # to TA and NA (by rpm) - keeping these separate as ranks should be separate 
genus_TA_TwoBP <- genus_TA_TwoBP %>%
  group_by(sample_name) %>%
  arrange(desc(nt_rpm)) %>%
  mutate(rank = row_number())

genus_NS_TwoBP <- genus_NS_TwoBP %>%
  group_by(sample_name) %>%
  arrange(desc(nt_rpm)) %>%
  mutate(rank = row_number())

#Merge on patient ID and genus combo, so we're only keeping rows from microbe reports that are the cultured pathogen 
genus_TA_TwoBP_pathogen <- genus_TA_TwoBP %>%
  inner_join(culture_results[, c(1,7)], by = c("manuscript_id" = "manuscript_id", "name" = "PNAOrg_genus"))

genus_NS_TwoBP_pathogen <- genus_NS_TwoBP %>%
  inner_join(culture_results[, c(1,7)], by = c("manuscript_id" = "manuscript_id", "name" = "PNAOrg_genus"))

#Merge with sample list and add zeros when there was no result 
genus_TA_TwoBP_pathogen <- genus_TA_TwoBP_pathogen %>%
  right_join(TA_TwoBP_samples[, c(1, 2, 3, 4, 5,6,7)], by = c("sample_name" = "sample_name", "manuscript_id" = "manuscript_id")) %>%
  tidyr::replace_na(list(nt_rpm=0)) %>%
  tidyr::replace_na(list(rank=200)) 

genus_NS_TwoBP_pathogen <- genus_NS_TwoBP_pathogen %>%
  right_join(NS_TwoBP_samples[, c(1, 2, 3, 4, 5,6,7,8)], by = c("sample_name" = "sample_name", "manuscript_id" = "manuscript_id")) %>%
  tidyr::replace_na(list(nt_rpm=0)) %>%
  tidyr::replace_na(list(rank=200)) 

genus_NS_TwoBP_pathogen$sample_type = 0
genus_TA_TwoBP_pathogen$sample_type = 1

######Recombine nasal swabs and trach aspirate samples 
genus_TwoBP_pathogen_all<-rbind(genus_TA_TwoBP_pathogen,genus_NS_TwoBP_pathogen)

#For visualization only, change all rank > 10 to be 11 for Figures
genus_TwoBP_pathogen_all$rank_new<-genus_TwoBP_pathogen_all$rank
for (i in (1:nrow(genus_TwoBP_pathogen_all))){
  if (genus_TwoBP_pathogen_all$rank[i] > 10)
  {genus_TwoBP_pathogen_all$rank_new[i] = 11}
} 

##Add mass data
pathogen_merge<-merge(genus_TwoBP_pathogen_all, bacterial_mass_filt, by.x="sample_name", by.y="sample_name", all.x=TRUE, all.y=FALSE)
###simplify columns
pathogen_merge<-pathogen_merge[,c(36,1,5,6,12,37,38,39,40, 41,42,44,45,46,47,48,49)]


#Calculate total mass
pathogen_merge$totalmass
for (i in (1:nrow(pathogen_merge))){
  pathogen_merge$totalmass[i] = ((pathogen_merge$total_reads[i]/pathogen_merge$total_ercc_reads[i])*pathogen_merge$ERCC_mass[i])*10^(-12)
} 

#Calculate (pathogen nt_rpm)
pathogen_merge$logpathogenrpm<-NA
for (i in (1:nrow(pathogen_merge))){
  pathogen_merge$pathogenrpm[i] = (pathogen_merge$nt_rpm[i])
} 

#Calculate (pathogen mass)
pathogen_merge$pathogenmass<-NA
for (i in (1:nrow(pathogen_merge))){
  pathogen_merge$pathogenmass[i] = ((pathogen_merge$nt_rpm[i]/pathogen_merge$total_ercc_reads[i])*pathogen_merge$ERCC_mass[i])*10^(-12)
} 

#Calculate (pathogen mass/totalmass)
pathogen_merge$normalizedpathmass<-NA
for (i in (1:nrow(pathogen_merge))){
  pathogen_merge$normalizedpathmass[i] = (pathogen_merge$pathogenmass[i]/pathogen_merge$totalmass[i])
} 

####For visualization only, create pathogen mass cutoffs so that very low reads still appear on the graph
pathogen_merge$pathogenmass_cutoff<-NA
for (i in (1:nrow(pathogen_merge))){
 (pathogen_merge$pathogenmass_cutoff[i] = (pathogen_merge$pathogenmass[i]*10^12))
} 

for (i in (1:nrow(pathogen_merge))){
  if (pathogen_merge$pathogenmass_cutoff[i] < 0.0001)
  {pathogen_merge$pathogenmass_cutoff_new[i] = 0.0001} else {pathogen_merge$pathogenmass_cutoff_new[i] = pathogen_merge$pathogenmass_cutoff[i]}
} 

####For visualization only, create pathogen mass ratio cutoffs so that very low reads still appear on the graphs
pathogen_merge$massratiocutoff<-NA
for (i in (1:nrow(pathogen_merge))){
  if ((pathogen_merge$normalizedpathmass[i]<10^-8))
  {pathogen_merge$massratiocutoff[i] = 10^-8} else {pathogen_merge$massratiocutoff[i] = (pathogen_merge$normalizedpathmass[i])}
} 


##To visualize one pathogen/patient combination for these graphs, divide each patient by pathogens
###for patient 1145/12, want 3 = Klebsiella, Escherichia, and Staphylococcus
####for patient 1250/29, want 2 = Klebsiella and Pseudomonas
## for patient 1474/45 want 2 - Staphylococcus and Haemophilus
for (i in (1:nrow(pathogen_merge))){
  if (pathogen_merge$manuscript_id[i] =="29" & pathogen_merge$name[i] =="Pseudomonas"  )
  {pathogen_merge$manuscript_id[i] = "29_1"}
  } 

for (i in (1:nrow(pathogen_merge))){
  if (pathogen_merge$manuscript_id[i] =="29" & pathogen_merge$name[i] =="Klebsiella"  )
  {pathogen_merge$manuscript_id[i] = "29_2"}
  } 

for (i in (1:nrow(pathogen_merge))){
  if (pathogen_merge$manuscript_id[i] =="12" & pathogen_merge$name[i] =="Escherichia"  )
  {pathogen_merge$manuscript_id[i] = "12_1"}
  } 

for (i in (1:nrow(pathogen_merge))){
  if (pathogen_merge$manuscript_id[i] =="12" & pathogen_merge$name[i] =="Klebsiella"  )
  {pathogen_merge$manuscript_id[i] = "12_2"}
  } 

for (i in (1:nrow(pathogen_merge))){
  if (pathogen_merge$manuscript_id[i] =="12" & pathogen_merge$name[i] =="Staphylococcus"  )
  {pathogen_merge$manuscript_id[i] = "12_3"}
  } 

for (i in (1:nrow(pathogen_merge))){
  if (pathogen_merge$manuscript_id[i] =="45" & pathogen_merge$name[i] =="Staphylococcus"  )
  {pathogen_merge$manuscript_id[i] = "45_1"}
  } 

for (i in (1:nrow(pathogen_merge))){
  if (pathogen_merge$manuscript_id[i] =="45" & pathogen_merge$name[i] =="Haemophilus"  )
  {pathogen_merge$manuscript_id[i] = "45_2"}
  } 

#Replace ptID 3 and 4 with 03 and 04 for ordering
for (i in (1:nrow(pathogen_merge))){
  if (pathogen_merge$manuscript_id[i] =="3")
  {pathogen_merge$manuscript_id[i] = "03"}
  } 

for (i in (1:nrow(pathogen_merge))){
  if (pathogen_merge$manuscript_id[i] =="4")
  {pathogen_merge$manuscript_id[i] = "04"}
}  

######For trach aspirate figures only, select out TA samples
TA_pathogen<-pathogen_merge %>% filter(sample_type==1)
##Supplemental figures with trach aspirate 
#Plot genus rank over time in each patient - note that here we consider ranks >10 as >12 for visualization 
print(p1<-ggplot(TA_pathogen, aes(x = sample_to_2BP, y = rank_new)) + geom_point(color = "#D41159", size=2) +
        geom_vline(xintercept=0, color = "black") +
        theme_bw() + facet_wrap(~ manuscript_id, nrow = 4) + labs(x="Days from 2BP", y="Bacterial genus rank") + 
        scale_y_continuous(trans = "reverse", breaks = c(11, 9, 7, 5, 3, 1), limits = c(12, 0)) + 
        theme(aspect.ratio = 1) +
        theme(panel.grid.minor = element_blank(), panel.background = element_blank()) +
        scale_x_continuous(breaks = c(-6, -4, -2, 0, 2, 4, 6), limits = c(-7, 7)))

###SupplementaryFigure8
ggsave(
  "Outputs_Final/TA_Rank_Plots_All.pdf",
  plot=p1,
)

####Total mass plot
print(p1<-ggplot(TA_pathogen, aes(x = sample_to_2BP, y = pathogenmass_cutoff_new)) + geom_point(color = "#D41159", size=2) +
        geom_vline(xintercept=0, color = "black") +
        theme_bw() + facet_wrap(~ manuscript_id, nrow=4) + labs(x="Days from 2BP", y="Pathogen mass") + 
        scale_y_log10(breaks = c(0.0001,0.01, 1,100, 10000), limits = c(0.000001, 10000)) +
        theme(aspect.ratio = 1) + 
        theme(panel.grid.minor = element_blank(), panel.background = element_blank()) +
        scale_x_continuous(breaks = c(-6, -4, -2, 0, 2, 4, 6), limits = c(-7, 7)))

###SupplementaryFigure9
ggsave(
  "Outputs_Final/TA_Mass_Plots_All.pdf",
  plot=p1,
)


  ####Normalized mass plot
print(p1<-ggplot(TA_pathogen, aes(x = sample_to_2BP, y =massratiocutoff )) + geom_point(color = "#D41159",size=2) +
        geom_vline(xintercept=0, color = "black") +
        theme_bw() + facet_wrap(~ manuscript_id, nrow=4) + labs(x="Days from 2BP", y="Pathogen mass/Total RNA Mass") + 
        scale_y_log10(limits = c(0.000000001, 0.01)) + 
        theme(aspect.ratio = 1) +
        theme(panel.grid.minor = element_blank(), panel.background = element_blank()) +
        scale_x_continuous(breaks = c(-6, -4, -2, 0, 2, 4, 6), limits = c(-7, 7)))

#SupplementaryFigure10
ggsave(
  "Outputs_Final/TA_NormalizedMass_Plots_All.pdf",
  plot=p1,
)


# TA and NS combined pathogen plots 
#Only include patientIDs that have NS data (13) OR 1233/24, which has a NS sample in which the bug is not detected
pathogen_merge_NS <- pathogen_merge %>%
  filter(manuscript_id %in% genus_NS_TwoBP_pathogen$manuscript_id | manuscript_id == "24"|manuscript_id == "12_1"|manuscript_id == "12_2"|manuscript_id == "12_3"
         |manuscript_id == "45_1"|manuscript_id == "45_2")

pathogen_merge_NS$sample_type<-as.factor(pathogen_merge_NS$sample_type)
##Plot genus rank over time in each patient - note that here we consider ranks >10 as >12 for visualization 
print(p1<-ggplot(pathogen_merge_NS, aes(x = sample_to_2BP, y = rank_new, color = sample_type)) + geom_point(size=2) +
        scale_color_manual(values = c("mediumpurple","#D41159")) +
        geom_vline(xintercept=0, color = "black") +
        theme_bw() + facet_wrap(~ manuscript_id, nrow = 3) + labs(x="Days from 2BP", y="Bacterial genus rank") + 
        scale_y_continuous(trans = "reverse", breaks = c(11, 9, 7, 5, 3, 1), limits = c(12, 0)) + 
        theme(aspect.ratio = 1) + 
        theme(panel.grid.minor = element_blank(), panel.background = element_blank()) +
        scale_x_continuous(breaks = c(-6, -4, -2, 0, 2, 4, 6), limits = c(-7, 7)))
##SupplementaryFigure12
ggsave(
  "Outputs_Final/NSandTA_Rank_Plots_All.pdf",
  plot=p1,
)


print(p1<-ggplot(pathogen_merge_NS, aes(x = sample_to_2BP, y = pathogenmass_cutoff_new, color = sample_type)) + geom_point(size=2) +
        scale_color_manual(values = c("mediumpurple","#D41159")) +
        geom_vline(xintercept=0, color = "black") +
        theme_bw() + facet_wrap(~ manuscript_id, nrow=3) + labs(x="Days from 2BP", y="Pathogen mass") + 
        scale_y_log10(breaks = c(0.0001,0.01, 1,100, 10000), limits = c(0.000001, 10000)) +
        theme(aspect.ratio = 1) + 
        theme(panel.grid.minor = element_blank(), panel.background = element_blank()) +
        scale_x_continuous(breaks = c(-6, -4, -2, 0, 2, 4, 6), limits = c(-7, 7)))
####SupplementaryFigure13
ggsave(
  "Outputs_Final/NSandTA_Mass_Plots_All.pdf",
  plot=p1,
)

####Normalized mass plot
print(p1<-ggplot(pathogen_merge_NS, aes(x = sample_to_2BP, y =massratiocutoff,color = sample_type)) + geom_point(size=2) +
        scale_color_manual(values = c("mediumpurple","#D41159")) +
        geom_vline(xintercept=0, color = "black") +
        theme_bw() + facet_wrap(~ manuscript_id, nrow=3) + labs(x="Days from 2BP", y="Pathogen mass/Total RNA Mass") + 
        scale_y_log10(limits = c(0.000000001, 0.01)) + 
        theme(aspect.ratio = 1) +
        theme(panel.grid.minor = element_blank(), panel.background = element_blank()) +
        scale_x_continuous(breaks = c(-6, -4, -2, 0, 2, 4, 6), limits = c(-7, 7)))
####SupplementaryFigure13
ggsave(
  "Outputs_Final/NSandTA_NormalizedMass_Plots_All.pdf",
  plot=p1,
)


#######FIGURE 3 - Example plots of tracheal aspirates 
####Limited to 4 patients - 1 increasing, 1 decreasing, 1 persistor, and 1 resistor 
path_limited<-TA_pathogen %>% filter (manuscript_id=="16"|manuscript_id=="23"|manuscript_id=="18"|manuscript_id=="13")
####Reorder for appearance
for (i in (1:nrow(path_limited))){
  if (path_limited$manuscript_id[i] =="16")
  {path_limited$manuscript_id[i] = "_16"}
} 

for (i in (1:nrow(path_limited))){
  if (path_limited$manuscript_id[i] =="23")
  {path_limited$manuscript_id[i] = "__23"}
} 

for (i in (1:nrow(path_limited))){
  if (path_limited$manuscript_id[i] =="18")
  {path_limited$manuscript_id[i] = "___18"}
} 

###Plot Rank
print(p1<-ggplot(path_limited, aes(x = sample_to_2BP, y = rank_new)) + geom_point(color = "#D41159", size=5) +
        geom_vline(xintercept=0, color = "black") +
        theme_bw() + facet_wrap(~ manuscript_id, nrow = 1) + labs(x="Days from 2BP", y="Bacterial genus rank") + 
        scale_y_continuous(trans = "reverse", breaks = c(11, 9, 7, 5, 3, 1), limits = c(12, 0)) + 
        theme(aspect.ratio = 1) +
        theme(panel.grid.minor = element_blank(), panel.background = element_blank()) +
        scale_x_continuous(breaks = c(-6, -4, -2, 0, 2, 4, 6), limits = c(-7, 7)))
####Figure 3 - Part 1
ggsave(
  "Outputs_Final/Fig3_part1.pdf",
  plot=p1,
)
print(p1<-ggplot(path_limited, aes(x = sample_to_2BP, y = pathogenmass_cutoff_new)) + geom_point(color = "#D41159", size=5) +
        geom_vline(xintercept=0, color = "black") +
        theme_bw() + facet_wrap(~ manuscript_id, nrow=1) + labs(x="Days from 2BP", y="Pathogen mass") + 
        scale_y_log10(limits = c(0.00001, 100)) + 
        theme(aspect.ratio = 1) + 
        theme(panel.grid.minor = element_blank(), panel.background = element_blank()) +
        scale_x_continuous(breaks = c(-6, -4, -2, 0, 2, 4, 6), limits = c(-7, 7)))
####Figure 3 - Part 2
ggsave(
  "Outputs_Final/Fig3_part2.pdf",
  plot=p1,
)
print(p1<-ggplot(path_limited, aes(x = sample_to_2BP, y =massratiocutoff)) + geom_point(color = "#D41159",size=5) +
        geom_vline(xintercept=0, color = "black") +
        theme_bw() + facet_wrap(~ manuscript_id, nrow=1) + labs(x="Days from 2BP", y="Pathogen mass/Total RNA Mass") + 
        scale_y_log10(limits = c(0.000000001, 0.01)) + 
        theme(aspect.ratio = 1) + 
        theme(panel.grid.minor = element_blank(), panel.background = element_blank()) +
        scale_x_continuous(breaks = c(-6, -4, -2, 0, 2, 4, 6), limits = c(-7, 7)))
####Figure 3 - Part 3
ggsave(
  "Outputs_Final/Fig3_part3.pdf",
  plot=p1,
)

#match with their rank - genus_TA_TwoBP_pathogen has the rank of the cultured pathogen in each patient at each timepoint 
genus_TA_TwoBP_pathogen_timepoint <- genus_TA_TwoBP_pathogen %>% filter(sample_name %in% TA_samples_timepoint$sample_name)

#####Figure 2 analyses comparing the tracheal aspirate microbiome of patients with 2BP with patients with NoBP 
##COVID analyses - Supp fig
genus_microbe_reports_covid <- genus_microbe_reports %>%
  filter(nt_rpm > 0.1 & nr_rpm > 0.1) %>%
  filter(name == "Betacoronavirus")

covid_timepoint<-merge(TA_samples_timepoint, genus_microbe_reports_covid, by.x="sample_name", by.y="sample_name", all.x=TRUE, all.y=FALSE)
p1<-ggplot(covid_timepoint, aes(x=TwoBP_Status, y = log(nt_rpm))) + geom_boxplot() + geom_jitter(width = 0.25) + theme_bw()   
print(p1)
covid_timepoint_short<-covid_timepoint[,c(1,2,4,5,18)]
colnames(covid_timepoint_short)[5]<-"BetaCoV_nt_rpm"
##saved to inputs, as this is used in host analyses
write.csv(Outputs_Final/covid_timepoint_short, "TA_SARS2_rpm.csv")



#Bacterial mass analyses - Figure 2A
mass_timepoint<-merge(TA_samples_timepoint, bacterial_mass_filt, by.x="sample_name", by.y="sample_name", all.x=TRUE, all.y=FALSE)
ggplot(mass_timepoint, aes(x=TwoBP_Status, y = log(bacterial_mass))) + geom_boxplot() + geom_jitter(width = 0.25) + theme_bw()     
mass_timepoint$bacterial_mass<-as.numeric(mass_timepoint$bacterial_mass)
wilcox.test(as.numeric(bacterial_mass) ~ TwoBP_Status, data = mass_timepoint)

#take cols of interest
mass_timepoint<-mass_timepoint[, c("manuscript_id", "TwoBP_Status", "bacterial_mass")]
write.csv(mass_timepoint, "Outputs_Final/TA_samples_timepoint_bacterial_mass.csv")

####Alpha diversity analyses - figures 2B-C
all_reports <- genus_microbe_reports %>%
  filter(sample_name %in% TA_samples$sample_name) %>% 
  filter(category == "bacteria") 

ID<-TA_samples[,c("manuscript_id", "sample_name")]
genus_TA_1<-merge(all_reports, ID, by.x ="sample_name", by.y ="sample_name", all.x=TRUE, all.y=FALSE)        

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
  left_join(TA_samples, by = c("sample_name" = "sample_name"))

##Output alpha diversity for supp figure of longitudinal alpha diversity over time - supplementary figure 1
write.csv(TA_samples_alpha, "Outputs_Final/Supp_Fig1_TA_samples_alpha_allsamples.csv")

#Pulling alpha diversity for timepoint analyses ##Figure 2B
TA_samples_alpha_timepoint<- TA_samples_alpha %>% filter(sample_name %in% TA_samples_timepoint$sample_name)
ggplot(TA_samples_alpha_timepoint, aes(x=TwoBP_Status, y = alpha_diversity)) + geom_boxplot() + geom_jitter(width = 0.25) + theme_bw()     
TA_samples_alpha_timepoint_forcsv<-TA_samples_alpha_timepoint[, c("manuscript_id", "TwoBP_Status", "alpha_diversity")]
#Write the alpha values into CSV files for figure 2B 
write.csv(TA_samples_alpha_timepoint_forcsv, "Outputs_Final/TA_samples_alpha_timepoint.csv")

###For Figure 2C - comparing alpha diversity between samples with cultured pathogen at highest rank vs controls
###Select the samples with the cultured pathogen at highest rank - note, that this INCLUDES
#the samples with multiple cultured pathogens as multiple pathogens 
best_TA_path <- genus_TA_TwoBP_pathogen[, c(1, 2, 5, 12, 13, 36, 37, 38, 39, 40, 41, 42)] %>%
  filter(sample_to_2BP <= 2 & sample_to_2BP >= -3) %>%
  group_by(manuscript_id) %>%
  dplyr::slice(which.min(rank)) 

#Remove the multiple pathogens from patients with multiple pathogens
best_TA_path_1<- best_TA_path %>% filter(manuscript_id != "12_2")
best_TA_path_1<- best_TA_path_1 %>% filter(manuscript_id != "12_3")
best_TA_path_1<- best_TA_path_1 %>% filter(manuscript_id != "29_2")
best_TA_path_1<- best_TA_path_1 %>% filter(manuscript_id != "45_2")          

##Add the control samples 
timepointcontrols<-TA_samples_timepoint %>% filter(TwoBP_Status == "No-BP")
best_TA_path_1$TwoBP_Status<- NA
#best_TA_path_1 <-best_TA_path_1 %>% rename( "dl_id"= "sample_name" )
best_TA_path_1$sample_name<-best_TA_path_1$sample_name

###Here we create a comparison of topranked 2BP samples with control samples 
##Figure 2C
TA_samples_toprank_timepoint<-rbind(timepointcontrols[, c(1,3,7,5)], best_TA_path_1[,c(6,1,11,10)])
TA_samples_alpha_topranked_timepoint<- merge(TA_samples_alpha, TA_samples_toprank_timepoint, by.x="sample_name", by.y="sample_name", all.x=FALSE, all.y=TRUE)
wilcox.test(alpha_diversity ~ TwoBP_Status.x, data = TA_samples_alpha_topranked_timepoint)
ggplot(TA_samples_alpha_topranked_timepoint, aes(x=TwoBP_Status.x, y = alpha_diversity)) + geom_boxplot() + geom_jitter(width = 0.25) + theme_bw()     
TA_samples_alpha_topranked_timepoint<-TA_samples_alpha_topranked_timepoint[, c("manuscript_id.x", "TwoBP_Status.x", "alpha_diversity")]
write.csv(TA_samples_alpha_topranked_timepoint, "Outputs_Final/TA_alpha_vs_rank_of_top_pathogen_at_toprankedsample.csv")


#####For supplementary figure 3 we compare alpha diversity by days on ventilator, days on steroids, SARS-CoV-2 RPM, and bacterial mass. 
#Join covid reads per million, bacterial mass, and pertinent metadata (days on vent and days on steroids before sample) to alpha diversity 
covid<-genus_microbe_reports_covid[,c(1,5,12)]
TA_samples_alpha_supp<-merge(TA_samples_alpha_timepoint, covid, by="sample_name", all.x=TRUE, all.y=FALSE)
TA_samples_alpha_supp<-merge(TA_samples_alpha_supp, bacterial_mass_filt, by="sample_name", all.x=TRUE, all.y=FALSE)
metadata_short<-metadata[,c(1,26)]
TA_samples_alpha_supp<-merge(TA_samples_alpha_supp, metadata_short, by="manuscript_id", all.x=TRUE, all.y=FALSE)
TA_samples_alpha_supp_2BP<-TA_samples_alpha_supp %>% filter (TwoBP_Status=="2BP")
TA_samples_alpha_supp_NoBP<-TA_samples_alpha_supp %>% filter (TwoBP_Status=="No-BP")

###divide days of steroids to sample into high-low
median(TA_samples_alpha_supp_2BP$Steroid_days_before_sample_collection)
for (i in (1:nrow(TA_samples_alpha_supp_2BP))){
  if (TA_samples_alpha_supp_2BP$Steroid_days_before_sample_collection[i] > 7)
  {TA_samples_alpha_supp_2BP$steroidshighlow[i] = "high"} else  {TA_samples_alpha_supp_2BP$steroidshighlow[i] = "low"}
} 

median(TA_samples_alpha_supp_NoBP$Steroid_days_before_sample_collection)
for (i in (1:nrow(TA_samples_alpha_supp_NoBP))){
  if (TA_samples_alpha_supp_NoBP$Steroid_days_before_sample_collection[i] > 2)
  {TA_samples_alpha_supp_NoBP$steroidshighlow[i] = "high"} else  {TA_samples_alpha_supp_NoBP$steroidshighlow[i] = "low"}
} 
wilcox.test(alpha_diversity ~ steroidshighlow, data = TA_samples_alpha_supp_2BP)
wilcox.test(alpha_diversity ~ steroidshighlow, data = TA_samples_alpha_supp_NoBP)

###divide days from intubation to sample into high-low
median(TA_samples_alpha_supp_2BP$vent_to_sample)
for (i in (1:nrow(TA_samples_alpha_supp_2BP))){
  if (TA_samples_alpha_supp_2BP$vent_to_sample[i] > 6)
  {TA_samples_alpha_supp_2BP$daysonventhighlow[i] = "high"} else  {TA_samples_alpha_supp_2BP$daysonventhighlow[i] = "low"}
} 

median(TA_samples_alpha_supp_NoBP$vent_to_sample)
for (i in (1:nrow(TA_samples_alpha_supp_NoBP))){
  if (TA_samples_alpha_supp_NoBP$vent_to_sample[i] > 6)
  {TA_samples_alpha_supp_NoBP$daysonventhighlow[i] = "high"} else  {TA_samples_alpha_supp_NoBP$daysonventhighlow[i] = "low"}
} 

wilcox.test(alpha_diversity ~ daysonventhighlow, data = TA_samples_alpha_supp_2BP)
wilcox.test(alpha_diversity ~ daysonventhighlow, data = TA_samples_alpha_supp_NoBP)

###divide covid into high-low
###replace NA wtih 0 (not detectable)
for (i in (1:nrow(TA_samples_alpha_supp_2BP))){
  if (is.na(TA_samples_alpha_supp_2BP$nt_rpm[i]) == TRUE)
  {TA_samples_alpha_supp_2BP$covidrpm[i] = 0} else  
    {TA_samples_alpha_supp_2BP$covidrpm[i] = (TA_samples_alpha_supp_2BP$nt_rpm[i])}
} 

for (i in (1:nrow(TA_samples_alpha_supp_NoBP))){
  if (is.na(TA_samples_alpha_supp_NoBP$nt_rpm[i]) == TRUE)
  {TA_samples_alpha_supp_NoBP$covidrpm[i] = 0} else  
  {TA_samples_alpha_supp_NoBP$covidrpm[i] = (TA_samples_alpha_supp_NoBP$nt_rpm[i])}
} 

###now divide into high low
median(TA_samples_alpha_supp_2BP$covidrpm)
for (i in (1:nrow(TA_samples_alpha_supp_2BP))){
  if (TA_samples_alpha_supp_2BP$covidrpm[i] > 15.91242)
  {TA_samples_alpha_supp_2BP$covidhighlow[i] = "high"} else  {TA_samples_alpha_supp_2BP$covidhighlow[i] = "low"}
} 

median(TA_samples_alpha_supp_NoBP$covidrpm)
for (i in (1:nrow(TA_samples_alpha_supp_NoBP))){
  if (TA_samples_alpha_supp_NoBP$covidrpm[i] > 93.20993)
  {TA_samples_alpha_supp_NoBP$covidhighlow[i] = "high"} else  {TA_samples_alpha_supp_NoBP$covidhighlow[i] = "low"}
} 
wilcox.test(alpha_diversity ~ covidhighlow, data = TA_samples_alpha_supp_2BP)
wilcox.test(alpha_diversity ~ covidhighlow, data = TA_samples_alpha_supp_NoBP)


####divide mass into high-low
median(TA_samples_alpha_supp_2BP$bacterial_mass)
for (i in (1:nrow(TA_samples_alpha_supp_2BP))){
  if (TA_samples_alpha_supp_2BP$bacterial_mass[i] > 90.93706)
  {TA_samples_alpha_supp_2BP$masshighlow[i] = "high"} else  {TA_samples_alpha_supp_2BP$masshighlow[i] = "low"}
} 
wilcox.test(alpha_diversity ~ masshighlow, data = TA_samples_alpha_supp_2BP)

median(TA_samples_alpha_supp_NoBP$bacterial_mass)
for (i in (1:nrow(TA_samples_alpha_supp_NoBP))){
  if (TA_samples_alpha_supp_NoBP$bacterial_mass[i] > 5.139369)
  {TA_samples_alpha_supp_NoBP$masshighlow[i] = "high"} else  {TA_samples_alpha_supp_NoBP$masshighlow[i] = "low"}
} 
wilcox.test(alpha_diversity ~ masshighlow, data = TA_samples_alpha_supp_NoBP)

####write raw data into supplementaries 
write.csv(TA_samples_alpha_supp_2BP, "Outputs_Final/SupplementaryFigAlphaDiversity_2BP_060324.csv")
write.csv(TA_samples_alpha_supp_NoBP, "Outputs_Final/SupplementaryFigAlphaDiversity_NoBP_060324.csv")



###Continuing tracheal aspirate analyses for Figure 2#####

###Figure 2E - generate two exemplary abundance plots
#Select two patients/samples from microbe data and create excel files with the name of detected microbes
#and nt_rpm
exemplary_12<-genus_microbe_reports %>% filter(sample_name == "MVIR1-HS145-D1ETA1-MNGS1")
exemplary_12<-exemplary_12[,c("name", "nt_rpm")]
write.csv(exemplary_12, "Outputs_Final/exemplary_12.csv")
exemplary_21<-genus_microbe_reports %>% filter(sample_name == "MVIR1-HS204-D0ETA1-MNGS1")
exemplary_21<-exemplary_21[,c("name", "nt_rpm")]
write.csv(exemplary_21, "Outputs_Final/exemplary_21.csv")

#For Figure 2F - analyzing of how many patients, is there a top-ranked cultured pathogen in 7 days before and after?
#for graphing, show rank >10 as 10
genus_TA_TwoBP_pathogen <- genus_TA_TwoBP_pathogen %>%
  mutate(rank_new = case_when(
    rank > 10 ~ 10,
    TRUE ~ rank
  ))

#remove the multiple pathogens per patient
genus_TA_TwoBP_pathogen_onepathperpt<- genus_TA_TwoBP_pathogen %>% filter(manuscript_id != "12_2")
genus_TA_TwoBP_pathogen_onepathperpt<- genus_TA_TwoBP_pathogen_onepathperpt %>% filter(manuscript_id != "12_3")
genus_TA_TwoBP_pathogen_onepathperpt<- genus_TA_TwoBP_pathogen_onepathperpt %>% filter(manuscript_id != "29_2")
genus_TA_TwoBP_pathogen_onepathperpt<- genus_TA_TwoBP_pathogen_onepathperpt %>% filter(manuscript_id != "45_2")

##this table asks what the rank is of the cultured pathogen across all samples, with 10 indicating 10 or greater
table(genus_TA_TwoBP_pathogen_onepathperpt$rank_new)

#This table asks where the highest rank of the cultured pathogen is across all samples per each patient
best_TA_sample <- genus_TA_TwoBP_pathogen_onepathperpt %>%
  group_by(manuscript_id) %>%
  dplyr::slice(which.min(rank)) 

table(best_TA_sample$rank)
write.csv(table(best_TA_sample$rank), "Outputs_Final/TA_toprankofculturedpathogen.csv")


###FIGURE 2G: Look at top ranked bacteria by mNGS at time of TwoBP vs No-BP 
###cut out some of the unnecessary cols 
TA_samples_timepoint_short<-TA_samples_timepoint[,c("sample_name", "manuscript_id", "TwoBP_Status")]

#Filter genus-level data for TA samples from culture+ TwoBP patients 
genus_TA_all <- genus_microbe_reports %>%
  filter(sample_name %in% TA_samples_timepoint_short$sample_name) %>% 
  filter(category == "bacteria") %>%
  left_join(TA_samples_timepoint_short, by = c("sample_name" = "sample_name"))

#Add rank # to genus (by rpm)
genus_TA_all_rank <- genus_TA_all %>%
  group_by(sample_name) %>%
  arrange(desc(nt_rpm)) %>%
  mutate(rank = row_number())

y<-genus_TA_all_rank[,c("sample_name", "name", "nt_rpm", "rank")]

#Pull out only the microbes ranked top one
topone<-y %>% filter(rank<2)

#join it with timepoint data - this becomes 2G
#Manual adjudication of S. epi and S. aureus performed 
topone_timepoint<-merge(TA_samples_timepoint_short,topone, by.x ="sample_name", by.y="sample_name")
table(topone_timepoint$name, topone_timepoint$TwoBP_Status)
topone_timepoint<-topone_timepoint[,c("manuscript_id", "TwoBP_Status", "name")]
write.csv(topone_timepoint, "Outputs_Final/toppathogenrankedat_TA_timepoint_sample.csv") 



### NASAL SWABS ### - FIGURE 4
###NS samples mass - FIgure 4A
mass_NS_timepoint<-merge(NS_samples_timepoint, bacterial_mass_filt, by.x="sample_name", by.y="sample_name", all.x=TRUE, all.y=FALSE)
ggplot(mass_NS_timepoint, aes(x=TwoBP_Status, y = log(bacterial_mass))) + geom_boxplot() + geom_jitter(width = 0.25) + theme_bw()     
wilcox.test(bacterial_mass ~ TwoBP_Status, data = mass_NS_timepoint)
mass_NS_timepoint<-mass_NS_timepoint[,c("manuscript_id", "TwoBP_Status", "bacterial_mass")]
write.csv(mass_NS_timepoint, "Outputs_Final/NS_samples_mass_timepoint.csv")

# Alpha diversity at time closest to TwoBP - this makes figure 4B
###Calculating NS Alpha diversity
# NS samples alpha diversity 
NS_all_timepoints <- genus_microbe_reports %>%
  filter(sample_name %in% NS_samples_timepoint$sample_name) %>%
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
  left_join(NS_samples_timepoint, by = c("sample_name" = "sample_name"))
#NS_samples_alpha <- read.csv()
ggplot(NS_samples_alpha, aes(x=TwoBP_Status, y = NS_alpha_diversity)) + geom_boxplot() + geom_jitter(width = 0.25) + theme_bw()     
#mann-whitney test: TwoBP vs No TwoBP at closest timepoint to TwoBP 
wilcox.test(NS_alpha_diversity ~ TwoBP_Status, data = NS_samples_alpha)
NS_samples_alpha<-NS_samples_alpha[,c("manuscript_id", "TwoBP_Status", "NS_alpha_diversity")]
write.csv(NS_samples_alpha, "Outputs_Final/NS_samples_alpha_timepoint.csv")

#Filter available nasal swab samples for TwoBP patients for pathogen analyses (Figure 4C)
NS_TwoBP_samples_filt <- NS_samples %>%  # 
  filter(TwoBP_Status == "2BP")
#Filter genus-level data for NS samples from culture+ TwoBP patients 
#Only keep samples that passed QC --> left with 15 TwoBP patients
genus_NS_TwoBP <- genus_microbe_reports %>%
  filter(sample_name %in% NS_TwoBP_samples_filt$sample_name) %>% 
  filter(category == "bacteria") %>%
  left_join(NS_TwoBP_samples_filt[, c(1, 2, 3,4, 5,6,8)], by = c("sample_name" = "sample_name"))

genus_NS_TwoBP <- genus_microbe_reports %>%
  filter(sample_name %in% NS_TwoBP_samples_filt$sample_name) %>% 
  filter(category == "bacteria") %>%
  left_join(NS_TwoBP_samples_filt, by = c("sample_name" = "sample_name"))

#Add rank # to genus (by rpm)
genus_NS_TwoBP <- genus_NS_TwoBP %>%
  group_by(sample_name) %>%
  arrange(desc(nt_rpm)) %>%
  mutate(rank = row_number())

#For each patient, filter for results for cultured pathogen
#Merge on COMET ID and pathogen name combo
genus_NS_TwoBP$manuscript_id<-as.integer(genus_NS_TwoBP$manuscript_id)
genus_NS_TwoBP_pathogen <- genus_NS_TwoBP %>%
  inner_join(culture_results, by = c("manuscript_id" = "manuscript_id", "name" = "PNAOrg_genus"))
table(genus_NS_TwoBP_pathogen$rank)        
#1233/22 drops out because there is no serratia detected in it at all, indicated as lower than limit of detection 


#Figure 4C
#Change all rank > 10 to show rank as 10
genus_NS_TwoBP_pathogen <- genus_NS_TwoBP_pathogen %>%
  mutate(rank_new = case_when(
    rank > 10 ~ 10,
    TRUE ~ rank
  ))
genus_NS_TwoBP_pathogen_fortable<- genus_NS_TwoBP_pathogen %>% filter(manuscript_id != "45_2")
genus_NS_TwoBP_pathogen_fortable<- genus_NS_TwoBP_pathogen_fortable %>% filter(manuscript_id != "12_2")
genus_NS_TwoBP_pathogen_fortable<- genus_NS_TwoBP_pathogen_fortable %>% filter(manuscript_id != "12_3")
table(genus_NS_TwoBP_pathogen_fortable$rank_new)

#Match all samples their rank - genus_TA_TwoBP_pathogen has the rank of the cultured pathogen in each patient at each timepoint 
genus_NS_TwoBP_pathogen_timepoint <- genus_NS_TwoBP_pathogen_fortable %>% filter(sample_name %in% NS_samples_timepoint$sample_name)
table(genus_NS_TwoBP_pathogen_timepoint$rank_new)

##Figure out of how many patients have the topranked pathogen in the 7 day window - figure 4C
##For 1233/22, the pathogen was never detected at all (hence >10)
best_NS_path <- genus_NS_TwoBP_pathogen_fortable %>%
  group_by(manuscript_id) %>%
  dplyr::slice(which.min(rank)) 
table(best_NS_path$rank_new)
write.csv(table(best_NS_path$rank_new), "Outputs_Final/NS_toprankofculturedpathogen.csv")


