# This code is used to analyze host genes for the COMET 2BP study

# Load packages -----
library(tidyverse)
library(ggplot2)
library(patchwork)

library(limma)
# library(edgeR)

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
  "Inputs/TA_samples_postbacterialANDhostQC_timepoint.csv",
  row.names = 1
)
# Convert some columns to factor
for (i in c("comet_id")) {
  metadata[,i] <- as.factor(metadata[,i])
}

# Check if there's any patient with multiple samples
stopifnot(!any(duplicated(metadata$comet_id)))

# Number of patients per group
print(table(metadata$Final))

# Days of ventilation per group
print(metadata %>%
        group_by(Final) %>%
        summarise(vent=median(vent_to_sample))
      )

# Plot distribution of ventilation time. They are a bit different
ggplot(data=metadata) +
  geom_histogram(aes(vent_to_sample)) +
  facet_wrap(~Final, ncol=1)

# Viral load
sc2 <- read.csv(
  "Inputs/BetaCoV_nt_rpM.csv"
)
# Add to metadata
metadata <- merge(
  metadata, sc2[,c("sample_name","nt_rpm")],
  by.x="dl_id", by.y="sample_name",
  all.x=TRUE, all.y=FALSE
)
# 3 samples don't have SARS-CoV-2, so they have NA nt_rpm
metadata[is.na(metadata$nt_rpm),"nt_rpm"] <- 0

## Count data -----
# Get the mapping from ensembl id to gene names
gene.names <- read.csv(
  "host_gene_name_mapping.csv",
  row.names=1
)

# Keep selected patients
host.counts <- read.csv(
  "Inputs/host_gene_counts.csv",
  row.names=1
)

# Filter for genes with >=10 counts in >=20% of samples
keep <- rowSums(host.counts >= 10) >= 0.2*ncol(host.counts)

## Load pathways for GSEA -----
# Hallmark gene sets
hallmark <- msigdbr::msigdbr(species="Homo sapiens",
                             category="H")
hallmark.list <- split(x=hallmark$ensembl_gene, f=hallmark$gs_name)
hallmark.list2 <- hallmark %>% select(gs_name, ensembl_gene)

# DE testing: host vs no-host -----
## Compare viral load -----
ggplot(data=metadata, aes(x=Final, y=log10(nt_rpm+1))) +
  geom_boxplot() +
  geom_jitter(width=0.2) +
  labs(x="", y="log10(rpM+1)", title="SARS-CoV-2 viral load vs host")

## Covariate: viral load -----
design <- model.matrix(~ Final + log10(nt_rpm+1), data = metadata)
print(colnames(design))
colnames(design)[2] <- "host.pos"

# limma-voom
vwts <- voom(host.counts[keep,], 
             design = design,
             normalize.method = "quantile",
             plot = T) 
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)
host.DE.sc2 <- topTable(vfit,
                   coef = "host.pos", sort.by = "none", 
                   number = Inf)
print("Effects of host on gene expression:")
print(sprintf("Number of genes with FDR < 0.1: %d", sum(host.DE.sc2$adj.P.Val<0.1)))

# Add gene names
host.DE.sc2$gene_name <- gene.names[rownames(host.DE.sc2), "gene_name"]

# Add a new column for coloring
host.DE.sc2.tmp <- host.DE.sc2 %>%
  mutate(tmp=ifelse(
    adj.P.Val>=0.1, "ns",
    ifelse(logFC>0, "up", "down")
  ))

# Volcano plot
p1 <- ggplot(data=host.DE.sc2.tmp, 
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
) # Fig 6A

### Box plot of highlighted genes -----
print(
  host.DE.sc2 %>% subset(gene_name %in% c("CST7","MUC16"))
)
# my.cpm <- edgeR::cpm(vwts, log=FALSE)
my.cpm <- vwts$E
temp <- my.cpm[rownames(gene.names)[gene.names$gene_name %in% c("CST7","MUC16")],]
rownames(temp) <- gene.names[rownames(temp),"gene_name"]
temp <- as.data.frame(t(temp))
temp <- merge(
  temp, metadata[,c("dl_id","comet_id","Final","nt_rpm")],
  by.x="row.names", by.y="dl_id",
  all.x=TRUE, all.y=FALSE, sort=FALSE
)

pp <- list()
for (i in c("CST7","MUC16")) {
  pp[[i]] <- ggplot(data=temp,
                    aes(x=Final, y=!!sym(i))) +
    geom_boxplot(aes(color=Final), outlier.shape=NA, width=0.5) +
    ggbeeswarm::geom_quasirandom(aes(color=Final), size=0.7) +
    scale_color_manual(values=c("#1A85FF","#D41159"), guide="none") +
    labs(x="", y="Normalized gene expression", title=i) +
    my.theme
}
ggsave(
  "bulk_host_boxplot_genes.svg",
  plot=pp[[1]]+pp[[2]],
  width=6.1, height=3.5, units="in"
) # Fig 6B

# Export results
write.csv(
  host.DE.sc2,
  "DE_host_vs_nohost_covariate-viralload.csv"
)

### GSEA: Hallmark -----
gene.ranks <- host.DE.sc2$t
names(gene.ranks) <- rownames(host.DE.sc2)
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
  "bulk_host-vs-nohost_gsea-hallmark.svg",
  plot=p2,
  width=4.5, height=3.5, units="in"
) # Fig 6C

# Export GSEA results
data.table::fwrite(
  hallmark.gsea.sc2,
  file="gsea-hallmark_host_vs_nohost_covariate-viralload.tsv",
  sep = "\t", sep2 = c("", " ", "")
)

# Compare gene expression vs bacterial mass (2BP patients only) -----
bact.mass <- read.csv(
  "Inputs/TA_samples_timepoint_bacterial_mass.csv",
  header=TRUE
)

# Add bacterial mass data
metadata.2 <- merge(
  metadata, bact.mass[,c("dl_id","mass")],
  by="dl_id",
  all.x=TRUE, all.y=FALSE
)
print(any(is.na(metadata.2$mass)))
boxplot(log10(mass)~Final, data=metadata.2)
wilcox.test(log10(mass)~Final, data=metadata.2)

# Isolate 2BP patients
metadata.2b <- metadata.2 %>%
  subset(Final=="2BP-culture pos")

# Compare SARS-CoV-2 viral load vs bacterial mass
plot(log10(nt_rpm+1) ~ log10(mass), data=metadata.2b,
     ylab="log10(SARS2 rpM + 1)")
print(summary(
  lm(log10(nt_rpm+1) ~ log10(mass), data=metadata.2b)
))
# viral load and bacterial mass is not correlated among 2BP patients

# DE analysis
design <- model.matrix(~ log10(mass), data = metadata.2b)
print(colnames(design))
# colnames(design)[2] <- "log10.mass"

# limma-voom
vwts <- voom(host.counts[keep,metadata.2b$dl_id], 
             design = design,
             normalize.method = "quantile",
             plot = T) 
vfit <- lmFit(vwts)
vfit <- eBayes(vfit)
host.DE.mass.b <- topTable(vfit,
                       coef = "log10(mass)", sort.by = "none", 
                       number = Inf)
print("Effects of bacterial mass on gene expression:")
print(sprintf("Number of genes with FDR < 0.1: %d", sum(host.DE.mass.b$adj.P.Val<0.1)))

# Add gene names
host.DE.mass.b$gene_name <- gene.names[rownames(host.DE.mass.b), "gene_name"]

# Add a new column for coloring
host.DE.mass.b.tmp <- host.DE.mass.b %>%
  mutate(tmp=ifelse(
    adj.P.Val>=0.1, "ns",
    ifelse(logFC>0, "up", "down")
  ))

# Volcano plot
p1 <- ggplot(data=host.DE.mass.b.tmp, 
             aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(aes(color=as.factor(tmp)), size=0.5) +
  coord_cartesian(xlim=c(-3,3), ylim=c(-0.1,3), expand=FALSE, clip="off") +
  scale_color_manual(
    values=c("ns"="grey","up"="#5ab4ac","down"="#d8b365"), guide="none") +
  labs(x="", y="") + # add x and y labels manually
  my.theme
ggsave(
  "bulk_host_mass_volcano.png",
  plot=p1, dpi=600, width=3.1, height=3.1, units="in"
) # Fig 6D

# Plot certain genes vs bacterial mass
tmp <- vwts$E[rownames(gene.names)[gene.names$gene_name %in% c("HLA-DRA","C1QC")],]
rownames(tmp) <- gene.names[rownames(tmp),"gene_name"]
tmp <- as.data.frame(t(tmp))
tmp$mass <- metadata.2b$mass
tmp$Final <- metadata.2b$Final

pp <- list()
for (i in c("HLA-DRA","C1QC")) {
  pp[[i]] <- ggplot(data=tmp,
                    aes(x=mass, y=!!sym(i))) +
    geom_point(color="#D41159") +
    scale_x_log10(
      limits=c(1e0,1e6), expand=c(0,0),
      labels=scales::trans_format("log10", scales::math_format(10^.x))) +
    geom_smooth(method='lm', formula= y ~ x,
                color="black", level=0.95) +
    labs(x="Bacterial mass (pg)", y="Normalized gene expression", title=i) +
    my.theme
}
pp[[1]] + pp[[2]]
ggsave(
  "bulk_host_mass_genes.svg",
  plot=pp[[1]]+pp[[2]],
  width=5.8, height=3.5, units="in"
) # Fig 6E

# GSEA Hallmark
gene.ranks <- host.DE.mass.b$t
names(gene.ranks) <- rownames(host.DE.mass.b)
gene.ranks <- sort(gene.ranks, decreasing = TRUE)

set.seed(1)
hallmark.gsea.mass <- fgseaMultilevel(
  pathways = hallmark.list,
  stats = gene.ranks,
  minSize = 15,
  maxSize = 500,
  nproc = 1
)

p4 <- ggplot(data=hallmark.gsea.mass,
             aes(x=NES, 
                 y=reorder(pathway,NES))) +
  geom_col(aes(color=as.factor(NES>0)), fill="white", width=0.7) +
  geom_col(data=. %>% subset(padj<0.05), aes(fill=as.factor(NES>0)), width=0.7) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "black") +
  scale_x_continuous(limits=c(-4,4), expand=c(0,0)) +
  labs(x="NES", y="Hallmark pathways",
       title="Pathways associated with bacterial mass\nin 2BP patients") +
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
  "bulk_host_mass_gsea-hallmark.svg",
  plot=p4,
  width=4.5, height=3.7, units="in"
) # Fig 6F

# Export DE results
write.csv(
  host.DE.mass,
  "DE_host_mass.csv"
)

# Export GSEA results
data.table::fwrite(
  hallmark.gsea.mass,
  file="gsea-hallmark_host_mass.tsv",
  sep = "\t", sep2 = c("", " ", "")
)

# Visually check gene expression vs bacterial mass
tmp <- vwts$E
plot(log10(metadata.2$mass), tmp["ENSG00000188536",])
