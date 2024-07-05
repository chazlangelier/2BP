# Author: Hoang Van Phan
# This script is for analyzing host gene expression from 2BP and No-BP patients

# Load packages -----
library(tidyverse)
library(ggplot2)
library(patchwork)

library(limma)

library(fgsea)

sessionInfo()

# ggplot theme
my.theme <- theme_classic() + 
  theme(
    plot.title = element_text(hjust = 0.5, size=12, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.3,1,0.7,0), "cm")
  )

# Load data -----
## Metadata -----
metadata <- read.csv(
  "./Inputs_Final/TA_samples_timepoint.csv",
  row.names = NULL
)
metadata$manuscript_id <- as.character(metadata$manuscript_id)
metadata$patient_id <- as.character(metadata$patient_id)
# Convert TwoBP_Status column to factor, keep No-BP as the reference level
metadata$TwoBP_Status <- factor(metadata$TwoBP_Status, levels=c("No-BP","2BP"))

# Check if there's any patient with multiple samples
stopifnot(!duplicated(metadata$patient_id))

# Number of patients per group
print(table(metadata$TwoBP_Status))

# Days of ventilation per group
print(metadata %>%
        group_by(TwoBP_Status) %>%
        summarise(vent=median(vent_to_sample))
      )

# Plot distribution of ventilation time. They are a bit different
ggplot(data=metadata) +
  geom_histogram(aes(vent_to_sample)) +
  facet_wrap(~TwoBP_Status, ncol=1)

# Viral load
sc2 <- read.csv(
  "./Inputs_Final/TA_SARS2_rpm.csv"
)
# Add to metadata
metadata <- merge(
  metadata, sc2[,c("sample_name","BetaCoV_nt_rpm")],
  by="sample_name",
  all.x=TRUE, all.y=FALSE
)
# 3 samples don't have SARS-CoV-2, so they have NA nt_rpm
metadata[is.na(metadata$BetaCoV_nt_rpm),"BetaCoV_nt_rpm"] <- 0

## Count data -----
# Mapping from ensembl id to gene names
gene.names <- read.csv(
  "./Inputs_Final/host_gene_name_mapping.csv",
  row.names=1)

# Host counts
host.counts <- read.csv(
  "./Inputs_Final/host_gene_counts.csv",
  row.names=1)
stopifnot(metadata$sample_name == colnames(host.counts))

# Filter for genes with >=10 counts in >=20% of samples
keep <- rowSums(host.counts >= 10) >= 0.2*ncol(host.counts)

## Load pathways for GSEA -----
# Hallmark gene sets
hallmark <- msigdbr::msigdbr(species="Homo sapiens",
                             category="H")
hallmark.list <- split(x=hallmark$ensembl_gene, f=hallmark$gs_name)
hallmark.list2 <- hallmark %>% select(gs_name, ensembl_gene)

# DE testing: 2BP vs no-BP -----
## Compare viral load -----
ggplot(data=metadata, aes(x=TwoBP_Status, y=log10(BetaCoV_nt_rpm+1))) +
  geom_boxplot() +
  geom_jitter(width=0.2) +
  labs(x="", y="log10(rpM+1)", title="SARS-CoV-2 viral load vs 2BP")

## Covariate: viral load -----
design <- model.matrix(~ TwoBP_Status + log10(BetaCoV_nt_rpm+1), data = metadata)
print(colnames(design))
colnames(design)[2] <- "status.2BP"

# limma-voom
vwts <- voom(host.counts[keep,], 
             design = design,
             normalize.method = "quantile",
             plot = T) 
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)
VAP.DE.sc2 <- topTable(
  vfit,
  coef = "status.2BP", sort.by = "none", 
  number = Inf)
print("Effects of 2BP on gene expression:")
print(sprintf("Number of genes with FDR < 0.1: %d", sum(VAP.DE.sc2$adj.P.Val<0.1)))

# Add gene names
VAP.DE.sc2$gene_name <- gene.names[rownames(VAP.DE.sc2), "gene_name"]

# Add a new column for coloring
VAP.DE.sc2.tmp <- VAP.DE.sc2 %>%
  mutate(tmp=ifelse(
    adj.P.Val>=0.1, "ns",
    ifelse(logFC>0, "up", "down")
  ))

# Volcano plot
p1 <- ggplot(data=VAP.DE.sc2.tmp, 
             aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(aes(color=as.factor(tmp)), size=0.5) +
  coord_cartesian(xlim=c(-4,4), ylim=c(-0.07,2), expand=FALSE) +
  scale_y_continuous(breaks=c(0,1,2)) +
  scale_color_manual(
    values=c("ns"="grey","up"="#D41159","down"="#1A85FF"), guide="none") +
  labs(x="", y="") + # add x and y labels manually
  my.theme
ggsave(
  "bulk_host_volcano_covariate-viralload.png",
  plot=p1, dpi=600, width=3.1, height=3.1, units="in"
)

# Export DE results
write.csv(
  VAP.DE.sc2,
  "DE_2BP-vs-noBP_covariate-viralload.csv"
)

### GSEA: Hallmark -----
gene.ranks <- VAP.DE.sc2$t
names(gene.ranks) <- rownames(VAP.DE.sc2)
gene.ranks <- sort(gene.ranks, decreasing = TRUE)

set.seed(1)
hallmark.gsea.sc2 <- fgseaMultilevel(
  pathways = hallmark.list,
  stats = gene.ranks,
  minSize = 15,
  maxSize = 500,
  nproc = 1
)

p2 <- ggplot(data=hallmark.gsea.sc2,
             aes(x=NES, 
                 y=reorder(pathway,NES))) +
  geom_col(aes(color=as.factor(NES>0)), fill="white", width=0.7) +
  geom_col(data=. %>% subset(padj<0.05), aes(fill=as.factor(NES>0)), width=0.7) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "black") +
  scale_x_continuous(limits=c(-2.5,2.5), expand=c(0,0)) +
  labs(x="NES", y="Hallmark pathways",
       title="Pathways associated with 2BP") +
  theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, size=12, face="plain"),
    axis.text.x = element_text(size=11, color="black"),
    axis.text.y = element_text(size=11, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0,1,1,0), "cm"),
    panel.grid.minor=element_blank(),
    # panel.grid.major.y = element_blank(),
    panel.spacing = unit(0.5, "cm"),
    legend.position = "none"
  )

# Select pathways
my.pathways <- c(
  "HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)
tmp <- hallmark.gsea.sc2 %>%
  subset(pathway %in% my.pathways)
# Fix pathway names
tmp$pathway <- gsub(
  "HALLMARK_","", tmp$pathway, fixed=TRUE)
tmp$pathway <- gsub(
  "_"," ", tmp$pathway, fixed=TRUE)
tmp$pathway <- stringr::str_to_sentence(tmp$pathway)
my.dict <- list("Il"="IL", "jak"="JAK", "stat"="STAT", "Tgf"="TGF",
                "Tnfa"="TNFa", "nfkb"="NF-kB", "Interferon"="IFN")
for (i in names(my.dict)) {
  tmp$pathway <- gsub(i, my.dict[[i]], tmp$pathway, fixed=TRUE)
}

p2 <- ggplot(data=tmp,
             aes(x=NES, 
                 y=reorder(pathway,NES))) +
  geom_col(aes(color=as.factor(NES>0)), fill="white", width=0.7) +
  geom_col(data=. %>% subset(padj<0.05), aes(fill=as.factor(NES>0)), width=0.7) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "black") +
  geom_text(
    data=. %>% subset(NES>0),
    aes(x=0,label=sprintf("P = %.2g", padj)), 
    hjust=1.07, size=3.4) +
  geom_text(
    data=. %>% subset(NES<0),
    aes(x=0, label=sprintf("P = %.2g", padj)),
    hjust=-0.07, size=3.4) +
  scale_x_continuous(limits=c(-2,2), expand=c(0,0)) +
  scale_color_manual(values=c("FALSE"="#1A85FF","TRUE"="#D41159")) +
  scale_fill_manual(values=c("FALSE"="#1A85FF","TRUE"="#D41159")) +
  labs(x="GSEA Normalized\nEnrichment Score", y="",
       title="Pathways associated with 2BP") +
  theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, size=12, face="plain"),
    axis.text.x = element_text(size=11, color="black"),
    axis.text.y = element_text(size=11, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.5,1,0.5,0), "cm"),
    panel.grid.minor=element_blank(),
    # panel.grid.major.y = element_blank(),
    panel.spacing = unit(0.5, "cm"),
    legend.position = "none"
  )
ggsave(
  "bulk_2BP-vs-noBP_gsea-hallmark.svg",
  plot=p2,
  width=4.5, height=3.5, units="in"
)

# Export GSEA results
data.table::fwrite(
  hallmark.gsea.sc2,
  file="gsea-hallmark_2BP-vs-noBP_covariate-viralload.tsv",
  sep = "\t", sep2 = c("", " ", "")
)

## Covariate: viral load (steroid patients only) -----
# Load steroid data
# Steroid_days_before_sample_collection: total days a patient was treated with steroids before sample collection
metadata.steroid <- read.csv(
  "./Inputs_Final/Patient_data_noPHI.csv"
) %>%
  select(manuscript_id,patient_id,TwoBP_Status,Steroid_days_before_sample_collection)

# Merge steroid info
metadata <- merge(
  metadata, metadata.steroid[,c("patient_id","Steroid_days_before_sample_collection")],
  by="patient_id",
  all.x=TRUE, all.y=FALSE, sort=FALSE
)
# All patients should have Steroid days
stopifnot(!is.na(metadata$Steroid_days_before_sample_collection))

# Extract samples who received steroid
metadata.steroid.yes <- metadata %>%
  subset(Steroid_days_before_sample_collection>0)

design <- model.matrix(~ TwoBP_Status + log10(BetaCoV_nt_rpm+1),
                       data = metadata.steroid.yes)
print(colnames(design))
colnames(design)[2] <- "status.2BP"

# limma-voom
vwts <- voom(host.counts[keep,metadata.steroid.yes$sample_name], 
             design = design,
             normalize.method = "quantile",
             plot = T) 
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)
VAP.DE.sc2.steroid <- topTable(
  vfit,
  coef = "status.2BP", sort.by = "none", 
  number = Inf)
print("Effects of 2BP on gene expression:")
print(sprintf("Number of genes with FDR < 0.1: %d",
              sum(VAP.DE.sc2.steroid$adj.P.Val<0.1)))

# Add gene names
VAP.DE.sc2.steroid$gene_name <- gene.names[rownames(VAP.DE.sc2.steroid), "gene_name"]

# Export DE results
write.csv(
  VAP.DE.sc2.steroid,
  "DE_2BP-vs-noBP_steroid_covariate-viralload.csv"
)

### GSEA: Hallmark -----
gene.ranks <- VAP.DE.sc2.steroid$t
names(gene.ranks) <- rownames(VAP.DE.sc2.steroid)
gene.ranks <- sort(gene.ranks, decreasing = TRUE)

set.seed(1)
hallmark.gsea.sc2.steroid <- fgseaMultilevel(
  pathways = hallmark.list,
  stats = gene.ranks,
  minSize = 15,
  maxSize = 500,
  nproc = 1
)

# Select pathways
my.pathways <- c(
  "HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)
tmp <- hallmark.gsea.sc2.steroid %>%
  subset(pathway %in% my.pathways)
# Fix pathway names
tmp$pathway <- gsub(
  "HALLMARK_","", tmp$pathway, fixed=TRUE
)
tmp$pathway <- gsub(
  "_"," ", tmp$pathway, fixed=TRUE
)
tmp$pathway <- stringr::str_to_sentence(tmp$pathway)
my.dict <- list("Il"="IL", "jak"="JAK", "stat"="STAT", "Tgf"="TGF",
                "Tnfa"="TNFa", "nfkb"="NF-kB", "Interferon"="IFN")
for (i in names(my.dict)) {
  tmp$pathway <- gsub(i, my.dict[[i]], tmp$pathway, fixed=TRUE)
}

p3 <- ggplot(data=tmp,
             aes(x=NES, 
                 y=reorder(pathway,NES))) +
  geom_col(aes(color=as.factor(NES>0)), fill="white", width=0.7) +
  geom_col(data=. %>% subset(padj<0.05), aes(fill=as.factor(NES>0)), width=0.7) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "black") +
  geom_text(
    data=. %>% subset(NES>0),
    aes(x=0,label=sprintf("P = %.2g", padj)), 
    hjust=1.07, size=3.4) +
  geom_text(
    data=. %>% subset(NES<0),
    aes(x=0, label=sprintf("P = %.2g", padj)),
    hjust=-0.07, size=3.4) +
  scale_x_continuous(limits=c(-2,2), expand=c(0,0)) +
  scale_color_manual(values=c("FALSE"="#1A85FF","TRUE"="#D41159")) +
  labs(x="GSEA Normalized\nEnrichment Score", y="",
       title="Pathways associated with 2BP\n(steroid recipients)") +
  theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, size=12, face="plain"),
    axis.text.x = element_text(size=11, color="black"),
    axis.text.y = element_text(size=11, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.5,1,0.5,0), "cm"),
    panel.grid.minor=element_blank(),
    # panel.grid.major.y = element_blank(),
    panel.spacing = unit(0.5, "cm"),
    legend.position = "none"
  )
ggsave(
  "bulk_2BP-vs-noBP_steroid_gsea-hallmark.svg",
  plot=p3,
  width=4.5, height=3.7, units="in"
)

# Export GSEA results
data.table::fwrite(
  hallmark.gsea.sc2.steroid,
  file="gsea-hallmark_2BP-vs-noBP_steroid_covariate-viralload.tsv",
  sep = "\t", sep2 = c("", " ", "")
)

# DE testing: continuous steroid day number -----
## In No-BP patients -----
metadata.noBP <- metadata %>%
  subset(TwoBP_Status=="No-BP")

design <- model.matrix(
  ~ Steroid_days_before_sample_collection,
  data = metadata.noBP)
print(colnames(design))

# limma-voom
vwts <- voom(host.counts[keep,metadata.noBP$sample_name], 
             design = design,
             normalize.method = "quantile",
             plot = T) 
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)
VAP.DE.steroid.days.a <- topTable(
  vfit,
  coef = "Steroid_days_before_sample_collection",
  sort.by = "none",
  number = Inf)
print("Effects of steroid days on gene expression:")
print(sprintf("Number of genes with FDR < 0.1: %d", sum(VAP.DE.steroid.days.a$adj.P.Val<0.1)))

# Add gene names
VAP.DE.steroid.days.a$gene_name <- gene.names[rownames(VAP.DE.steroid.days.a), "gene_name"]
# Export DE results
write.csv(
  VAP.DE.steroid.days.a,
  "DE_noBP_steroid-days.csv"
)

### GSEA: Hallmark -----
gene.ranks <- VAP.DE.steroid.days.a$t
names(gene.ranks) <- rownames(VAP.DE.steroid.days.a)
gene.ranks <- sort(gene.ranks, decreasing = TRUE)

set.seed(1)
hallmark.gsea.steroid.days <- fgseaMultilevel(
  pathways = hallmark.list,
  stats = gene.ranks,
  minSize = 15,
  maxSize = 500,
  nproc = 1
)

# Export GSEA results
data.table::fwrite(
  hallmark.gsea.steroid.days,
  file="gsea-hallmark_noBP_steroid-days.tsv",
  sep = "\t", sep2 = c("", " ", "")
)

# Select pathways
my.pathways <- c(
  "HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)
tmp <- hallmark.gsea.steroid.days %>%
  subset(pathway %in% my.pathways)
# Fix pathway names
tmp$pathway <- gsub(
  "HALLMARK_","", tmp$pathway, fixed=TRUE)
tmp$pathway <- gsub(
  "_"," ", tmp$pathway, fixed=TRUE)
tmp$pathway <- stringr::str_to_sentence(tmp$pathway)
my.dict <- list("Il"="IL", "jak"="JAK", "stat"="STAT", "Tgf"="TGF",
                "Tnfa"="TNFa", "nfkb"="NF-kB", "Interferon"="IFN")
for (i in names(my.dict)) {
  tmp$pathway <- gsub(i, my.dict[[i]], tmp$pathway, fixed=TRUE)
}

p.a <- ggplot(data=tmp,
             aes(x=NES, 
                 y=reorder(pathway,NES))) +
  geom_col(aes(color=as.factor(NES>0)), fill="white", width=0.6) +
  geom_col(data=. %>% subset(padj<0.05), aes(fill=as.factor(NES>0)), width=0.6) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "black") +
  geom_text(
    data=. %>% subset(NES>0),
    aes(x=0,label=sprintf("P = %.2g", padj)), 
    hjust=1.07, size=3.4) +
  geom_text(
    data=. %>% subset(NES<0),
    aes(x=0, label=sprintf("P = %.2g", padj)),
    hjust=-0.07, size=3.4) +
  scale_x_continuous(limits=c(-4,4), expand=c(0,0)) +
  scale_color_manual(values=c("FALSE"="#5e3c99","TRUE"="#e66101")) +
  scale_fill_manual(values=c("FALSE"="#5e3c99","TRUE"="#e66101")) +
  labs(x="GSEA Normalized\nEnrichment Score", y="",
       title="Pathways associated with steroid days\nin No-BP patients") +
  theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, size=12, face="plain"),
    axis.text.x = element_text(size=11, color="black"),
    axis.text.y = element_text(size=11, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.5,1,0.5,0), "cm"),
    panel.grid.minor=element_blank(),
    # panel.grid.major.y = element_blank(),
    panel.spacing = unit(0.5, "cm"),
    legend.position = "none"
  )

## In 2BP patients -----
metadata.2BP <- metadata %>%
  subset(TwoBP_Status=="2BP")

design <- model.matrix(
  ~ Steroid_days_before_sample_collection,
  data = metadata.2BP)
print(colnames(design))

# limma-voom
vwts <- voom(host.counts[keep,metadata.2BP$sample_name], 
             design = design,
             normalize.method = "quantile",
             plot = T) 
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)
VAP.DE.steroid.days.b <- topTable(
  vfit,
  coef = "Steroid_days_before_sample_collection",
  sort.by = "none", 
  number = Inf)
print("Effects of steroid days on gene expression:")
print(sprintf("Number of genes with FDR < 0.1: %d", sum(VAP.DE.steroid.days.b$adj.P.Val<0.1)))

# Add gene names
VAP.DE.steroid.days.b$gene_name <- gene.names[rownames(VAP.DE.steroid.days.b), "gene_name"]
# Export DE results
write.csv(
  VAP.DE.steroid.days.b,
  "DE_2BP_steroid-days.csv"
)

### GSEA: Hallmark -----
gene.ranks <- VAP.DE.steroid.days.b$t
names(gene.ranks) <- rownames(VAP.DE.steroid.days.b)
gene.ranks <- sort(gene.ranks, decreasing = TRUE)

set.seed(1)
hallmark.gsea.steroid.days <- fgseaMultilevel(
  pathways = hallmark.list,
  stats = gene.ranks,
  minSize = 15,
  maxSize = 500,
  nproc = 1
)

# Export GSEA results
data.table::fwrite(
  hallmark.gsea.steroid.days,
  file="gsea-hallmark_2BP_steroid-days.tsv",
  sep = "\t", sep2 = c("", " ", "")
)

# Select pathways
my.pathways <- c(
  "HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)
tmp <- hallmark.gsea.steroid.days %>%
  subset(pathway %in% my.pathways)
# Fix pathway names
tmp$pathway <- gsub(
  "HALLMARK_","", tmp$pathway, fixed=TRUE
)
tmp$pathway <- gsub(
  "_"," ", tmp$pathway, fixed=TRUE
)
tmp$pathway <- stringr::str_to_sentence(tmp$pathway)
my.dict <- list("Il"="IL", "jak"="JAK", "stat"="STAT", "Tgf"="TGF",
                "Tnfa"="TNFa", "nfkb"="NF-kB", "Interferon"="IFN")
for (i in names(my.dict)) {
  tmp$pathway <- gsub(i, my.dict[[i]], tmp$pathway, fixed=TRUE)
}

p.b <- ggplot(data=tmp,
              aes(x=NES, 
                  y=reorder(pathway,NES))) +
  geom_col(aes(color=as.factor(NES>0)), fill="white", width=0.6) +
  geom_col(data=. %>% subset(padj<0.05), aes(fill=as.factor(NES>0)), width=0.6) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "black") +
  geom_text(
    data=. %>% subset(NES>0),
    aes(x=0,label=sprintf("P = %.2g", padj)), 
    hjust=1.07, size=3.4) +
  geom_text(
    data=. %>% subset(NES<0),
    aes(x=0, label=sprintf("P = %.2g", padj)),
    hjust=-0.07, size=3.4) +
  scale_x_continuous(limits=c(-4,4), expand=c(0,0)) +
  scale_color_manual(values=c("FALSE"="#5e3c99","TRUE"="#e66101")) +
  scale_fill_manual(values=c("FALSE"="#5e3c99","TRUE"="#e66101")) +
  labs(x="GSEA Normalized\nEnrichment Score", y="",
       title="Pathways associated with steroid days\nin 2BP patients") +
  theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, size=12, face="plain"),
    axis.text.x = element_text(size=11, color="black"),
    axis.text.y = element_text(size=11, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.5,1,0.5,0), "cm"),
    panel.grid.minor=element_blank(),
    # panel.grid.major.y = element_blank(),
    panel.spacing = unit(0.5, "cm"),
    legend.position = "none"
  )

ggsave(
  "bulk_2BP_noBP_steroid-days_gsea-hallmark.svg",
  plot=p.b + p.a,
  width=9.2, height=3.7, units="in"
)

# Compare gene expression vs bacterial mass -----
bact.mass <- read.csv(
  "./Inputs_Final/TA_samples_timepoint_bacterial_mass.csv",
  header=TRUE
)
bact.mass$manuscript_id <- as.character(bact.mass$manuscript_id)

# Add bacterial mass data
metadata.2 <- merge(
  metadata, bact.mass[,c("manuscript_id","bacterial_mass")],
  by="manuscript_id",
  all.x=TRUE, all.y=FALSE,
  sort=FALSE
)
print(any(is.na(metadata.2$bacterial_mass)))
boxplot(log10(bacterial_mass) ~ TwoBP_Status, data=metadata.2)
wilcox.test(log10(bacterial_mass) ~ TwoBP_Status, data=metadata.2)

## In No-BP patients -----
metadata.2a <- metadata.2 %>%
  subset(TwoBP_Status=="No-BP")

design <- model.matrix(~ log10(bacterial_mass), data = metadata.2a)
print(colnames(design))

# limma-voom
vwts <- voom(host.counts[keep,metadata.2a$sample_name], 
             design = design,
             normalize.method = "quantile",
             plot = T) 
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)
VAP.DE.mass.a <- topTable(
  vfit,
  coef = "log10(bacterial_mass)",
  sort.by = "none",
  number = Inf)
print("Effects of bacterial mass on gene expression:")
print(sprintf("Number of genes with FDR < 0.1: %d", sum(VAP.DE.mass.a$adj.P.Val<0.1)))

# Add gene names
VAP.DE.mass.a$gene_name <- gene.names[rownames(VAP.DE.mass.a), "gene_name"]

# Add a new column for coloring
VAP.DE.mass.a.tmp <- VAP.DE.mass.a %>%
  mutate(tmp=ifelse(
    adj.P.Val>=0.1, "ns",
    ifelse(logFC>0, "up", "down")
  ))

# Volcano plot
p1 <- ggplot(data=VAP.DE.mass.a.tmp, 
             aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(aes(color=as.factor(tmp)), size=0.5) +
  coord_cartesian(xlim=c(-3,3), ylim=c(-0.1,3), expand=FALSE, clip="off") +
  scale_color_manual(
    values=c("ns"="grey","up"="#5ab4ac","down"="#d8b365"), guide="none") +
  labs(x="", y="") + # add x and y labels manually
  my.theme
ggsave(
  "bulk_noBP_mass_volcano.png",
  plot=p1, dpi=600, width=2.9, height=3.1, units="in"
)

# Export DE results
write.csv(
  VAP.DE.mass.a,
  "DE_noBP_mass.csv"
)

# GSEA Hallmark
gene.ranks <- VAP.DE.mass.a$t
names(gene.ranks) <- rownames(VAP.DE.mass.a)
gene.ranks <- sort(gene.ranks, decreasing = TRUE)

set.seed(1)
hallmark.gsea.mass <- fgseaMultilevel(
  pathways = hallmark.list,
  stats = gene.ranks,
  minSize = 15,
  maxSize = 500,
  nproc = 1
)

# Export GSEA results
data.table::fwrite(
  hallmark.gsea.mass,
  file="gsea-hallmark_noBP_mass.tsv",
  sep = "\t", sep2 = c("", " ", "")
)

p3 <- ggplot(data=hallmark.gsea.mass,
             aes(x=NES, 
                 y=reorder(pathway,NES))) +
  geom_col(aes(color=as.factor(NES>0)), fill="white", width=0.6) +
  geom_col(data=. %>% subset(padj<0.05), aes(fill=as.factor(NES>0)), width=0.6) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "black") +
  scale_x_continuous(limits=c(-4,4), expand=c(0,0)) +
  labs(x="NES", y="Hallmark pathways",
       title="Pathways associated with bacterial mass\nin No-BP patients") +
  theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, size=12, face="plain"),
    axis.text.x = element_text(size=11, color="black"),
    axis.text.y = element_text(size=11, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0,1,1,0), "cm"),
    panel.grid.minor=element_blank(),
    # panel.grid.major.y = element_blank(),
    panel.spacing = unit(0.5, "cm"),
    legend.position = "none"
  )

# Select pathways
my.pathways <- c(
  "HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)
tmp <- hallmark.gsea.mass %>%
  subset(pathway %in% my.pathways)
# Fix pathway names
tmp$pathway <- gsub(
  "HALLMARK_","", tmp$pathway, fixed=TRUE)
tmp$pathway <- gsub(
  "_"," ", tmp$pathway, fixed=TRUE)
tmp$pathway <- stringr::str_to_sentence(tmp$pathway)
my.dict <- list("Il"="IL", "jak"="JAK", "stat"="STAT", "Tgf"="TGF",
                "Tnfa"="TNFa", "nfkb"="NF-kB", "Interferon"="IFN")
for (i in names(my.dict)) {
  tmp$pathway <- gsub(i, my.dict[[i]], tmp$pathway, fixed=TRUE)
}

p4 <- ggplot(data=tmp,
             aes(x=NES, 
                 y=reorder(pathway,NES))) +
  geom_col(aes(color=as.factor(NES>0)), fill="white", width=0.7) +
  geom_col(data=. %>% subset(padj<0.05), aes(fill=as.factor(NES>0)), width=0.7) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "black") +
  geom_text(
    data=. %>% subset(NES>0),
    aes(x=0,label=sprintf("P = %.2g", padj)), 
    hjust=1.07, size=3.4) +
  geom_text(
    data=. %>% subset(NES<0),
    aes(x=0, label=sprintf("P = %.2g", padj)),
    hjust=-0.07, size=3.4) +
  scale_x_continuous(limits=c(-4,4), expand=c(0,0)) +
  scale_color_manual(values=c("FALSE"="#d8b365","TRUE"="#5ab4ac")) +
  scale_fill_manual(values=c("FALSE"="#d8b365","TRUE"="#5ab4ac")) +
  labs(x="GSEA Normalized\nEnrichment Score", y="",
       title="Pathways associated with bacterial mass\nin No-BP patients") +
  theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, size=12, face="plain"),
    axis.text.x = element_text(size=11, color="black"),
    axis.text.y = element_text(size=11, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.5,1,0.5,0), "cm"),
    panel.grid.minor=element_blank(),
    # panel.grid.major.y = element_blank(),
    panel.spacing = unit(0.5, "cm"),
    legend.position = "none"
  )
ggsave(
  "bulk_noBP_mass_gsea-hallmark.svg",
  plot=p4,
  width=4.3, height=3.6, units="in"
)

## In 2BP patients -----
metadata.2b <- metadata.2 %>%
  subset(TwoBP_Status=="2BP")

# Compare SARS-CoV-2 viral load vs bacterial mass
plot(log10(BetaCoV_nt_rpm+1) ~ log10(bacterial_mass), data=metadata.2b,
     ylab="log10(SARS2 rpM + 1)")
print(summary(
  lm(log10(BetaCoV_nt_rpm+1) ~ log10(bacterial_mass), data=metadata.2b)
))
# viral load and bacterial mass is not correlated among 2BP patients

# DE analysis
design <- model.matrix(~ log10(bacterial_mass), data = metadata.2b)
print(colnames(design))

# limma-voom
vwts <- voom(host.counts[keep,metadata.2b$sample_name], 
             design = design,
             normalize.method = "quantile",
             plot = T) 
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)
VAP.DE.mass.b <- topTable(
  vfit,
  coef = "log10(bacterial_mass)",
  sort.by = "none",
  number = Inf)
print("Effects of bacterial mass on gene expression:")
print(sprintf("Number of genes with FDR < 0.1: %d", sum(VAP.DE.mass.b$adj.P.Val<0.1)))

# Add gene names
VAP.DE.mass.b$gene_name <- gene.names[rownames(VAP.DE.mass.b), "gene_name"]

# Add a new column for coloring
VAP.DE.mass.b.tmp <- VAP.DE.mass.b %>%
  mutate(tmp=ifelse(
    adj.P.Val>=0.1, "ns",
    ifelse(logFC>0, "up", "down")
  ))

# Volcano plot
p1 <- ggplot(data=VAP.DE.mass.b.tmp, 
             aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(aes(color=as.factor(tmp)), size=0.5) +
  coord_cartesian(xlim=c(-3,3), ylim=c(-0.1,3), expand=FALSE, clip="off") +
  scale_color_manual(
    values=c("ns"="grey","up"="#5ab4ac","down"="#d8b365"), guide="none") +
  labs(x="", y="") + # add x and y labels manually
  my.theme
ggsave(
  "bulk_2BP_mass_volcano.png",
  plot=p1, dpi=600, width=2.9, height=3.1, units="in"
)

# Plot certain genes vs bacterial mass
tmp <- vwts$E[rownames(gene.names)[gene.names$gene_name %in% c("HLA-DRA","C1QC")],]
rownames(tmp) <- gene.names[rownames(tmp),"gene_name"]
tmp <- as.data.frame(t(tmp))
tmp$bacterial_mass <- metadata.2b$bacterial_mass

pp <- list()
for (i in c("HLA-DRA","C1QC")) {
  pp[[i]] <- ggplot(data=tmp,
                    aes(x=bacterial_mass, y=!!sym(i))) +
    geom_point(color="#D41159") +
    scale_x_log10(
      # limits=c(1e0,1e6), expand=c(0,0),
      breaks=c(1e0,1e1,1e3,1e5),
      labels=scales::trans_format("log10", scales::math_format(10^.x))) +
    geom_smooth(method='lm', formula= y ~ x,
                color="black", level=0.95) +
    labs(x="Bacterial mass (pg)", y="Normalized gene expression", title=i) +
    my.theme
}
pp[[1]] + pp[[2]]
ggsave(
  "bulk_2BP_mass_genes.svg",
  plot=pp[[1]]+pp[[2]],
  width=5.8, height=3.5, units="in"
)

# GSEA Hallmark
gene.ranks <- VAP.DE.mass.b$t
names(gene.ranks) <- rownames(VAP.DE.mass.b)
gene.ranks <- sort(gene.ranks, decreasing = TRUE)

set.seed(1)
hallmark.gsea.mass <- fgseaMultilevel(
  pathways = hallmark.list,
  stats = gene.ranks,
  minSize = 15,
  maxSize = 500,
  nproc = 1
)

# Select pathways
my.pathways <- c(
  "HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)
tmp <- hallmark.gsea.mass %>%
  subset(pathway %in% my.pathways)
# Fix pathway names
tmp$pathway <- gsub(
  "HALLMARK_","", tmp$pathway, fixed=TRUE)
tmp$pathway <- gsub(
  "_"," ", tmp$pathway, fixed=TRUE)
tmp$pathway <- stringr::str_to_sentence(tmp$pathway)
my.dict <- list("Il"="IL", "jak"="JAK", "stat"="STAT", "Tgf"="TGF",
                "Tnfa"="TNFa", "nfkb"="NF-kB", "Interferon"="IFN")
for (i in names(my.dict)) {
  tmp$pathway <- gsub(i, my.dict[[i]], tmp$pathway, fixed=TRUE)
}

p4 <- ggplot(data=tmp,
             aes(x=NES, 
                 y=reorder(pathway,NES))) +
  geom_col(aes(color=as.factor(NES>0)), fill="white", width=0.6) +
  geom_col(data=. %>% subset(padj<0.05), aes(fill=as.factor(NES>0)), width=0.6) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "black") +
  geom_text(
    data=. %>% subset(NES>0),
    aes(x=0,label=sprintf("P = %.2g", padj)), 
    hjust=1.07, size=3.4) +
  geom_text(
    data=. %>% subset(NES<0),
    aes(x=0, label=sprintf("P = %.2g", padj)),
    hjust=-0.07, size=3.4) +
  scale_x_continuous(limits=c(-4,4), expand=c(0,0)) +
  scale_color_manual(values=c("FALSE"="#d8b365","TRUE"="#5ab4ac")) +
  scale_fill_manual(values=c("FALSE"="#d8b365","TRUE"="#5ab4ac")) +
  labs(x="GSEA Normalized\nEnrichment Score", y="",
       title="Pathways associated with bacterial mass\nin 2BP patients") +
  theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, size=12, face="plain"),
    axis.text.x = element_text(size=11, color="black"),
    axis.text.y = element_text(size=11, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.5,1,0.5,0), "cm"),
    panel.grid.minor=element_blank(),
    # panel.grid.major.y = element_blank(),
    panel.spacing = unit(0.5, "cm"),
    legend.position = "none"
  )
ggsave(
  "bulk_2BP_mass_gsea-hallmark.svg",
  plot=p4,
  width=4.3, height=3.6, units="in"
)

# Export DE results
write.csv(
  VAP.DE.mass.b,
  "DE_2BP_mass.csv"
)

# Export GSEA results
data.table::fwrite(
  hallmark.gsea.mass,
  file="gsea-hallmark_2BP_mass.tsv",
  sep = "\t", sep2 = c("", " ", "")
)

## In No-BP steroid-receiving patients -----
# Plot bacterial mass vs steroid days
p5 <- ggplot(
  data=metadata.2,
  aes(x=Steroid_days_before_sample_collection,
      y=log10(bacterial_mass))) +
  geom_point(aes(color=TwoBP_Status), size=0.8) +
  geom_smooth(method='lm', formula= y~x,
              se=TRUE, level=0.95,
              color="black") +
  labs(x="Number of days on steroid\nbefore sample collection",
       y="log10(bacterial mass)") +
  scale_color_manual(
    values=c("2BP"="#D41159","No-BP"="#1A85FF"), guide="none") +
  facet_wrap(~fct_relevel(TwoBP_Status,"2BP","No-BP")) +
  theme_bw() +
  theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, size=12, face="plain"),
    axis.text = element_text(size=12, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.5,1,0.5,0), "cm"),
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    # panel.grid.major.y = element_blank(),
    panel.spacing = unit(0.5, "cm"),
    legend.position = "none"
  )

metadata.2c <- metadata.2 %>%
  subset((TwoBP_Status=="No-BP") & (Steroid_days_before_sample_collection>0))

design <- model.matrix(~ log10(bacterial_mass), data = metadata.2c)
print(colnames(design))

# limma-voom
vwts <- voom(host.counts[keep,metadata.2c$sample_name], 
             design = design,
             normalize.method = "quantile",
             plot = T) 
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)
VAP.DE.mass.c <- topTable(
  vfit,
  coef = "log10(bacterial_mass)",
  sort.by = "none",
  number = Inf)
print("Effects of bacterial mass on gene expression:")
print(sprintf("Number of genes with FDR < 0.1: %d", sum(VAP.DE.mass.c$adj.P.Val<0.1)))

# Add gene names
VAP.DE.mass.c$gene_name <- gene.names[rownames(VAP.DE.mass.c), "gene_name"]

# Export DE results
write.csv(
  VAP.DE.mass.c,
  "DE_noBP_mass_steroid.csv"
)

# GSEA Hallmark
gene.ranks <- VAP.DE.mass.c$t
names(gene.ranks) <- rownames(VAP.DE.mass.c)
gene.ranks <- sort(gene.ranks, decreasing = TRUE)

set.seed(1)
hallmark.gsea.mass <- fgseaMultilevel(
  pathways = hallmark.list,
  stats = gene.ranks,
  minSize = 15,
  maxSize = 500,
  nproc = 1
)

# Export GSEA results
data.table::fwrite(
  hallmark.gsea.mass,
  file="gsea-hallmark_noBP_mass_steroid.tsv",
  sep = "\t", sep2 = c("", " ", "")
)

# Select pathways
my.pathways <- c(
  "HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)
tmp <- hallmark.gsea.mass %>%
  subset(pathway %in% my.pathways)
# Fix pathway names
tmp$pathway <- gsub(
  "HALLMARK_","", tmp$pathway, fixed=TRUE
)
tmp$pathway <- gsub(
  "_"," ", tmp$pathway, fixed=TRUE
)
tmp$pathway <- stringr::str_to_sentence(tmp$pathway)
my.dict <- list("Il"="IL", "jak"="JAK", "stat"="STAT", "Tgf"="TGF",
                "Tnfa"="TNFa", "nfkb"="NF-kB", "Interferon"="IFN")
for (i in names(my.dict)) {
  tmp$pathway <- gsub(i, my.dict[[i]], tmp$pathway, fixed=TRUE)
}

p4 <- ggplot(data=tmp,
             aes(x=NES, 
                 y=reorder(pathway,NES))) +
  geom_col(aes(color=as.factor(NES>0)), fill="white", width=0.6) +
  geom_col(data=. %>% subset(padj<0.05), aes(fill=as.factor(NES>0)), width=0.6) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "black") +
  geom_text(
    data=. %>% subset(NES>0),
    aes(x=0,label=sprintf("P = %.2g", padj)), 
    hjust=1.07, size=3.4) +
  geom_text(
    data=. %>% subset(NES<0),
    aes(x=0, label=sprintf("P = %.2g", padj)),
    hjust=-0.07, size=3.4) +
  scale_x_continuous(limits=c(-4,4), expand=c(0,0)) +
  scale_color_manual(values=c("FALSE"="#d8b365","TRUE"="#5ab4ac")) +
  scale_fill_manual(values=c("FALSE"="#d8b365","TRUE"="#5ab4ac")) +
  labs(x="GSEA Normalized\nEnrichment Score", y="",
       title="Pathways associated with bacterial mass\nin No-BP steroid recipients") +
  theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, size=12, face="plain"),
    axis.text.x = element_text(size=11, color="black"),
    axis.text.y = element_text(size=11, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.5,1,0.5,0), "cm"),
    panel.grid.minor=element_blank(),
    # panel.grid.major.y = element_blank(),
    panel.spacing = unit(0.5, "cm"),
    legend.position = "none"
  )
ggsave(
  "bulk_noBP_steroid_mass_gsea-hallmark.svg",
  plot=p4,
  width=4.3, height=3.6, units="in"
)

## In 2BP steroid-receiving patients -----
# Plot bacterial mass vs steroid days
metadata.2d <- metadata.2 %>%
  subset((TwoBP_Status=="2BP") & (Steroid_days_before_sample_collection>0))

design <- model.matrix(~ log10(bacterial_mass), data = metadata.2d)
print(colnames(design))

# limma-voom
vwts <- voom(host.counts[keep,metadata.2d$sample_name], 
             design = design,
             normalize.method = "quantile",
             plot = T) 
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)
VAP.DE.mass.d <- topTable(vfit,
                          coef = "log10(bacterial_mass)", sort.by = "none", 
                          number = Inf)
print("Effects of bacterial mass on gene expression:")
print(sprintf("Number of genes with FDR < 0.1: %d", sum(VAP.DE.mass.d$adj.P.Val<0.1)))

# Add gene names
VAP.DE.mass.d$gene_name <- gene.names[rownames(VAP.DE.mass.d), "gene_name"]

# Add a new column for coloring
VAP.DE.mass.d.tmp <- VAP.DE.mass.d %>%
  mutate(tmp=ifelse(
    adj.P.Val>=0.1, "ns",
    ifelse(logFC>0, "up", "down")
  ))

# Export DE results
write.csv(
  VAP.DE.mass.d,
  "DE_2BP_mass_steroid.csv"
)

# GSEA Hallmark
gene.ranks <- VAP.DE.mass.d$t
names(gene.ranks) <- rownames(VAP.DE.mass.d)
gene.ranks <- sort(gene.ranks, decreasing = TRUE)

set.seed(1)
hallmark.gsea.mass <- fgseaMultilevel(
  pathways = hallmark.list,
  stats = gene.ranks,
  minSize = 15,
  maxSize = 500,
  nproc = 1
)

# Export GSEA results
data.table::fwrite(
  hallmark.gsea.mass,
  file="gsea-hallmark_2BP_mass_steroid.tsv",
  sep = "\t", sep2 = c("", " ", "")
)

# Select pathways
my.pathways <- c(
  "HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)
tmp <- hallmark.gsea.mass %>%
  subset(pathway %in% my.pathways)
# Fix pathway names
tmp$pathway <- gsub(
  "HALLMARK_","", tmp$pathway, fixed=TRUE
)
tmp$pathway <- gsub(
  "_"," ", tmp$pathway, fixed=TRUE
)
tmp$pathway <- stringr::str_to_sentence(tmp$pathway)
my.dict <- list("Il"="IL", "jak"="JAK", "stat"="STAT", "Tgf"="TGF",
                "Tnfa"="TNFa", "nfkb"="NF-kB", "Interferon"="IFN")
for (i in names(my.dict)) {
  tmp$pathway <- gsub(i, my.dict[[i]], tmp$pathway, fixed=TRUE)
}

p4 <- ggplot(data=tmp,
             aes(x=NES, 
                 y=reorder(pathway,NES))) +
  geom_col(aes(color=as.factor(NES>0)), fill="white", width=0.6) +
  geom_col(data=. %>% subset(padj<0.05), aes(fill=as.factor(NES>0)), width=0.6) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "black") +
  geom_text(
    data=. %>% subset(NES>0),
    aes(x=0,label=sprintf("P = %.2g", padj)), 
    hjust=1.07, size=3.4) +
  geom_text(
    data=. %>% subset(NES<0),
    aes(x=0, label=sprintf("P = %.2g", padj)),
    hjust=-0.07, size=3.4) +
  scale_x_continuous(limits=c(-4,4), expand=c(0,0)) +
  scale_color_manual(values=c("FALSE"="#d8b365","TRUE"="#5ab4ac")) +
  scale_fill_manual(values=c("FALSE"="#d8b365","TRUE"="#5ab4ac")) +
  labs(x="GSEA Normalized\nEnrichment Score", y="",
       title="Pathways associated with bacterial mass\nin 2BP steroid recipients") +
  theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, size=12, face="plain"),
    axis.text.x = element_text(size=11, color="black"),
    axis.text.y = element_text(size=11, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.5,1,0.5,0), "cm"),
    panel.grid.minor=element_blank(),
    # panel.grid.major.y = element_blank(),
    panel.spacing = unit(0.5, "cm"),
    legend.position = "none"
  )
ggsave(
  "bulk_2BP_steroid_mass_gsea-hallmark.svg",
  plot=p4,
  width=4.3, height=3.6, units="in"
)

## In No-BP non-steroid patients -----
metadata.2e <- metadata.2 %>%
  subset((TwoBP_Status=="No-BP") & (Steroid_days_before_sample_collection==0))

design <- model.matrix(~ log10(bacterial_mass), data = metadata.2e)
print(colnames(design))

# limma-voom
vwts <- voom(host.counts[keep,metadata.2e$sample_name], 
             design = design,
             normalize.method = "quantile",
             plot = T) 
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)
VAP.DE.mass.e <- topTable(
  vfit,
  coef = "log10(bacterial_mass)",
  sort.by = "none",
  number = Inf)
print("Effects of bacterial mass on gene expression:")
print(sprintf("Number of genes with FDR < 0.1: %d", sum(VAP.DE.mass.e$adj.P.Val<0.1)))

# Add gene names
VAP.DE.mass.e$gene_name <- gene.names[rownames(VAP.DE.mass.e), "gene_name"]

# Export DE results
write.csv(
  VAP.DE.mass.e,
  "DE_noBP_mass_non-steroid.csv"
)

# GSEA Hallmark
gene.ranks <- VAP.DE.mass.e$t
names(gene.ranks) <- rownames(VAP.DE.mass.e)
gene.ranks <- sort(gene.ranks, decreasing = TRUE)

set.seed(1)
hallmark.gsea.mass <- fgseaMultilevel(
  pathways = hallmark.list,
  stats = gene.ranks,
  minSize = 15,
  maxSize = 500,
  nproc = 1
)

# Export GSEA results
data.table::fwrite(
  hallmark.gsea.mass,
  file="gsea-hallmark_noBP_mass_non-steroid.tsv",
  sep = "\t", sep2 = c("", " ", "")
)

# Select pathways
my.pathways <- c(
  "HALLMARK_IL2_STAT5_SIGNALING","HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)
tmp <- hallmark.gsea.mass %>%
  subset(pathway %in% my.pathways)
# Fix pathway names
tmp$pathway <- gsub(
  "HALLMARK_","", tmp$pathway, fixed=TRUE
)
tmp$pathway <- gsub(
  "_"," ", tmp$pathway, fixed=TRUE
)
tmp$pathway <- stringr::str_to_sentence(tmp$pathway)
my.dict <- list("Il"="IL", "jak"="JAK", "stat"="STAT", "Tgf"="TGF",
                "Tnfa"="TNFa", "nfkb"="NF-kB", "Interferon"="IFN")
for (i in names(my.dict)) {
  tmp$pathway <- gsub(i, my.dict[[i]], tmp$pathway, fixed=TRUE)
}

p6 <- ggplot(data=tmp,
             aes(x=NES, 
                 y=reorder(pathway,NES))) +
  geom_col(aes(color=as.factor(NES>0)), fill="white", width=0.6) +
  geom_col(data=. %>% subset(padj<0.05), aes(fill=as.factor(NES>0)), width=0.6) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "black") +
  geom_text(
    data=. %>% subset(NES>0),
    aes(x=0,label=sprintf("P = %.2g", padj)), 
    hjust=1.07, size=3.4) +
  geom_text(
    data=. %>% subset(NES<0),
    aes(x=0, label=sprintf("P = %.2g", padj)),
    hjust=-0.07, size=3.4) +
  scale_x_continuous(limits=c(-4,4), expand=c(0,0)) +
  scale_color_manual(values=c("FALSE"="#d8b365","TRUE"="#5ab4ac")) +
  scale_fill_manual(values=c("FALSE"="#d8b365","TRUE"="#5ab4ac")) +
  labs(x="GSEA Normalized\nEnrichment Score", y="",
       title="Pathways associated with bacterial mass\nin No-BP non-steroid recipients") +
  theme_bw() + theme(
    plot.title = element_text(hjust = 0.5, size=12, face="plain"),
    axis.text.x = element_text(size=11, color="black"),
    axis.text.y = element_text(size=11, color="black"),
    text = element_text(size=12, family="Arial"),
    plot.margin = unit(c(0.5,1,0.5,0), "cm"),
    panel.grid.minor=element_blank(),
    # panel.grid.major.y = element_blank(),
    panel.spacing = unit(0.5, "cm"),
    legend.position = "none"
  )
ggsave(
  "bulk_noBP_non-steroid_mass_gsea-hallmark.svg",
  plot=p6,
  width=4.3, height=3.6, units="in"
)
