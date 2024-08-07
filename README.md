# 2BP: Secondary Bacteria Pneumonia in COVID-19

This repository contains the code for the paper **Microbial Dynamics and Pulmonary Immune Responses in COVID-19 Secondary Bacterial Pneumonia**.

# Table of contents

1. [Code](#code)
2. [Code usage](#code-usage)
3. [Figures](#figures)
    1. [Main Figures](#main-figures)
    2. [Supplementary Figures](#supplementary-figures)
4. [Supplementary Tables](#supplementary-tables)
5. [Data files](#data-files)
6. [Required hardware and software](#required-hardware-and-software)

## Code

[2BPCode_final.R](2BPCode_final.R) is used to analyze data at the microbe level and produce manuscript Figures 2A-C, 2E-G, Figure 3, and  Figure 4A-4C, as well as Supplementary Figures 1, 3, 6, and 8-14. This code uses the following input files: 
* [TA_samples_longitudinal.csv](Inputs_Final/TA_samples_longitudinal.csv): a list of tracheal aspirate samples for all patients at all timepoints available, after bacterial QC.
* [TA_samples_timepoint.csv](Inputs_Final/TA_samples_timepoint.csv): a list of tracheal aspirate samples at timepoint closest to 2BP diagnosis, or, for No-BP controls, at matched days of mechanical ventilation, after host and bacterial QC.
* [NS_samples_longitudinal.csv](Inputs_Final/NS_samples_longitudinal.csv): a list of nasal swab samples for all patients at all timepoints available, after bacterial QC.
* [NS_samples_timepoint.csv](Inputs_Final/NS_samples_timepoint.csv): a list of nasal swab samples at timepoint closest to 2BP diagnosis, or, for No-BP controls, at matched days of mechanical ventilation, after host and bacterial QC.
* [Patient_data_noPHI.csv](Inputs_Final/Patient_data_noPHI.csv): patient metadata, deidentified.
* [2BP_CultureResults.csv](Inputs_Final/2BP_CultureResults.csv): culture results from conventional microbiologic assays for patients with 2BP. 
* [Microbe_reports_genus.csv](Inputs_Final/Microbe_reports_genus.csv): outputs from CZID pipeline.
* [2BP_Mass_Read_Output.csv](Inputs_Final/2BP_Mass_Read_Output.csv): bacterial mass for each sample.

[2BP_AMRCode.Rmd](2BP_AMRCode.Rmd) is used to make manuscript Figure 2D, 2H, 4D, 5, and Supplementary Figures 2, 4, 5, and 15-17. This code uses the following additional input files: 
* [AMR_genes_of_interest.csv](Inputs_Final/AMR_genes_of_interest.csv): AMR genes of public health significance evaluted in Supplementary Figure 16.
* [ARGANNOT_genelist.csv](Inputs_Final/ARGANNOT_genelist.csv): AMR genes in the ARGANNOT database.
* [idseq_amr_heatmap_values_final.csv](Inputs_Final/amr_heatmap_values.csv): output of CZ ID AMR pipeline (which uses SRST2 and ARGANNOT database) for this dataset.

[2BP_host_gene_analysis.R](2BP_host_gene_analysis.R) is used to analyze data at the human gene expression data and make Figures 6 and 7, as well as Supplementary Data Files 3-14. This code uses the following files:
* [TA_samples_timepoint.csv](Inputs_Final/TA_samples_timepoint.csv): a list of tracheal aspirate samples at timepoint closest to 2BP diagnosis, or, for No-BP controls, at matched days of mechanical ventilation, after host and bacterial QC.
* [TA_SARS2_rpm.csv](Inputs_Final/TA_SARS2_rpm.csv): SARS-CoV-2 viral abundance, in reads per million output from the CZ ID pipeline. Generated by 2BPCode_final.R 
* [host_gene_name_mapping.csv](Inputs_Final/host_gene_name_mapping.csv): for conversion between gene symbol and gene ensembl ID.
* [host_gene_counts.csv](Input_Finals/host_gene_counts.csv): host gene count table.
* [Patient_data_noPHI.csv](Inputs_Final/Patient_data_noPHI.csv): steroid information.
* [TA_samples_timepoint_bacterial_mass.csv](Inputs_Final/TA_samples_timepoint_bacterial_mass.csv): Table of bacterial mass.

[2BP_microbial_pathways_analysis.Rmd](2BP_microbial_pathways_analysis.Rmd) is used to analyze microbial functional pathways and make Supplementary Figure 7. This code uses the following files:
* [TA_samples_timepoint.csv](Inputs_Final/TA_samples_timepoint.csv): a list of tracheal aspirate samples at timepoint closest to 2BP diagnosis, or, for No-BP controls, at matched days of mechanical ventilation, after host and bacterial QC. 
* [humann_pathabundance_relab.tsv](Inputs_Final/humann_pathabundance_relab.tsv). This file was produced by renormalizing and joining HUMAnN's outputs with the commands `humann_renorm_table()` and `humann_join_tables()`.

## Code usage

The full data sets are available in the folder [Inputs_Final](Inputs_Final). To reproduce the results in the paper, the uploaded codes can be run directly. The code outputs are located in the folder [Outputs_Final](Outputs_Final).

## Figures

### Main Figures

Figure 1. Figure 1a created in Biorender. Figure 1b made in Prism and source data is available in Supplementary Data 1. 

Figure 2: csv output files for panels A, B, C, E, G were generated by 2BPCode_final.R, and then final graphs created in Prism. Figure 2D was created in R using 2BP_AMRCode.Rmd and final panel arrangement performed in Prism (standardizing font size).  Figure 2F csv output file was generated by 2BPCode_final and plot was created in Prism. Figure 2H csv output files was generated by 2BP_AMRCode.Rmd and plot was created in Prism. 
* A: [TA_samples_timepoint_bacterial_mass.csv](Outputs_Final/TA_samples_timepoint_bacterial_mass.csv)
* B: [TA_samples_alpha_timepoint.csv](Outputs_Final/TA_samples_alpha_timepoint.csv)
* C: [TA_alpha_vs_rank_of_top_pathogen_at_toprankedsample.csv](Outputs_Final/TA_alpha_vs_rank_of_top_pathogen_at_toprankedsample.csv)
* D: [Figure2D.pdf](Outputs_Final/Figure2D.pdf)
* E: [exemplary_12.csv](Outputs_Final/exemplary_12.csv) and [exemplary_21.csv](Outputs_Final/exemplary_21.csv)
* F: [toppathogenrankedat_TA_timepoint_sample.csv](Outputs_Final/toppathogenrankedat_TA_timepoint_sample.csv) 
* G: [TA_toprankofculturedpathogen.csv](Outputs_Final/TA_toprankofculturedpathogen.csv)
* H: [TA_DA data - TA2BPvsNoBP.csv](<Outputs_Final/TA_DA data - TA2BPvsNoBP.csv>)

Figure 3: Graph generated by 2BPCode_final.R and with final edits performed in Illustrator (adding pathogen names, adding empty circles for the samples in which pathogen was not detected, adding antibiotic data)
* [Fig3_part1.pdf](Outputs_Final/Fig3_part1.pdf), [Fig3_part2.pdf](Outputs_Final/Fig3_part2.pdf) and [Fig3_part3.pdf](Outputs_Final/Fig3_part3.pdf)

Figure 4: csv output files for panels A, B, C were generated by 2BPCode_final.R, and then final graphs created in Prism. Figure 4D was created as a graph in R using code 2BP_AMRCode.Rmd and final edits performed in Illustrator (combining the 2BP and no BP plots, adding statistical test results performed in R using code 2BP_AMRCode.Rmd, adjusting font). 
* A: [NS_samples_mass_timepoint.csv](Outputs_Final/NS_samples_mass_timepoint.csv)
* B: [NS_samples_alpha_timepoint.csv](Outputs_Final/NS_samples_alpha_timepoint.csv)
* C: [NS_toprankofculturedpathogen.csv](Outputs_Final/NS_toprankofculturedpathogen.csv)
* D: [Fig4d_barplot.pdf](Outputs_Final/Fig4d_barplot.pdf)

Figure 5: Created as a graph in R using code 2BP_AMRCode.Rmd and finalized in Illustrator (standardizing font sizes) 
* [Fig5-AMR_heatmap_v2.pdf](Outputs_Final/Fig5-AMR_heatmap_v2.pdf)

Figure 6: generated by the code 2BP_host_gene_analysis.R. The panels were exported under the following names, and edited in Inkscape (adding axis labels, removing white spaces, editing panel titles, simplyfing small P-values):
* A: [bulk_host_volcano_covariate-viralload.png](Outputs_Final/bulk_host_volcano_covariate-viralload.png)
* B: [bulk_2BP-vs-noBP_gsea-hallmark.svg](Outputs_Final/bulk_2BP-vs-noBP_gsea-hallmark.svg)
* C: [bulk_2BP-vs-noBP_steroid_gsea-hallmark.svg](Outputs_Final/bulk_2BP-vs-noBP_steroid_gsea-hallmark.svg)
* D, E: [bulk_2BP_noBP_steroid-days_gsea-hallmark.svg](Outputs_Final/bulk_2BP_noBP_steroid-days_gsea-hallmark.svg)

Figure 7: generated by the code 2BP_host_gene_analysis.R. The panels were exported under the following names, and edited in Inkscape (adding axis labels, removing white spaces, editing panel titles, simplifying small P-values):
* A: [bulk_2BP_mass_volcano.png](Outputs_Final/bulk_2BP_mass_volcano.png)
* B: [bulk_2BP_mass_genes.svg](Outputs_Final/bulk_2BP_mass_genes.svg)
* C: [bulk_noBP_mass_volcano.png](Outputs_Final/bulk_noBP_mass_volcano.png)
* D: [bulk_2BP_mass_gsea-hallmark.svg](Outputs_Final/bulk_2BP_mass_gsea-hallmark.svg)
* E: [bulk_noBP_mass_gsea-hallmark.svg](Outputs_Final/bulk_noBP_mass_gsea-hallmark.svg)
* F: [bulk_2BP_steroid_mass_gsea-hallmark.svg](Outputs_Final/bulk_2BP_steroid_mass_gsea-hallmark.svg)
* G: [bulk_noBP_steroid_mass_gsea-hallmark.svg](Outputs_Final/bulk_noBP_steroid_mass_gsea-hallmark.svg)
* H: [bulk_noBP_non-steroid_mass_gsea-hallmark.svg](Outputs_Final/bulk_noBP_non-steroid_mass_gsea-hallmark.svg)

### Supplementary Figures

Supp Figure 1: generated by the code 2BPCode_final.R. The csv files were generated by 2BPCode_final.R and then final graph created in Prism.
* [Supp_Fig1_TA_samples_alpha_allsamples.csv](Outputs_Final/Supp_Fig1_TA_samples_alpha_allsamples.csv)

Supp Figure 2: generated by the code 2BP_AMRCode.R and finalized in Illustrator (standardizing font sizes)
* A: [SuppFig-Jaccard beta div.pdf](<Outputs_Final/SuppFig-Jaccard beta div.pdf>)

Supp Figure 3:  generated by the code 2BPCode_final.R. The csv files were generated by 2BPCode_final.R and then final graph created in Prism.
* [SupplementaryFigAlphaDiversity_2BP.csv](Outputs_Final/SupplementaryFigAlphaDiversity_2BP.csv) and [SupplementaryFigAlphaDiversity_NoBP.csv](Outputs_Final/SupplementaryFigAlphaDiversity_NoBP.csv)

Supp Figure 4: generated by the code 2BP_AMRCode.R and finalized in Illustrator (standardizing font sizes)
* A: [SuppFig - betadiv contsteroids.pdf](<Outputs_Final/SuppFig - betadiv contsteroids.pdf>)
* B: [SuppFig - betadiv vent day.pdf](<Outputs_Final/SuppFig - betadiv vent day.pdf>)
* C: [SuppFig - betadiv adjSC2.pdf](<Outputs_Final/SuppFig - betadiv adjSC2.pdf>)
* D: [SuppFig - betadiv adjbactmass.pdf](<Outputs_Final/SuppFig - betadiv adjbactmass.pdf>)

Supp Figure 5: generated by the code 2BP_AMRCode.R. The csv files were generated by 2BP_AMRCode.R and plotted in Prism.
* [TA_DA data - TA2BPvsNoBP_steroids only.csv](<Inputs_Final/TA_DA data - TA2BPvsNoBP_steroids only.csv>)

Supp Figure 6: generated by the code 2BPCode_final.R. The csv files were generated by 2BPCode_final.R and then final plot created in Prism. 
* [TA_SARS2_rpm.csv](Inputs_Final/TA_SARS2_rpm.csv)

Supp Figure 7: generated by the code 2BP_microbial_pathways_analysis.Rmd.
* [microbial_pathway_fdr0.25.svg](Outputs_Final/microbial_pathway_fdr0.25.svg)

Supp Figure 8: generated by the code 2BPCode_final.R and finalized in Illustrator (standardizing font sizes, adding patient numbers and pathogen ID, adding antibiotic data)
* [TA_Rank_Plots_All.pdf](Outputs_Final/TA_Rank_Plots_All.pdf)

Supp Figure 9: generated by the code 2BPCode_final.R and finalized in Illustrator (standardizing font sizes, adding patient numbers and pathogen ID, adding antibiotic data)
* [TA_Mass_Plots_All.pdf](Outputs_Final/TA_Mass_Plots_All.pdf)

Supp Figure 10: generated by code 2BPCode_final.R and finalized in Illustrator (standardizing font sizes, adding patient numbers and pathogen ID, adding antibiotic data)
* [TA_NormalizedMass_Plots_All.pdf](Outputs_Final/TA_NormalizedMass_Plots_All.pdf)

Supp Figure 11: plot produced in Prism.
* [Persistor_notpersistor.csv](Outputs_Final/Persistor_notpersistor.csv)

Supp Figure 12: generated by the code 2BPCode_final.R and finalized in Illustrator (standardizing font sizes, adding patient numbers and pathogen ID, adding antibiotic data)
* [NSandTA_Rank_Plots_All.pdf](Outputs_Final/NSandTA_Rank_Plots_All.pdf)

Supp Figure 13: generated by the code 2BPCode_final.R and finalized in Illustrator (standardizing font sizes, adding patient numbers and pathogen ID, adding antibiotic data)
* [NSandTA_Mass_Plots_All.pdf](Outputs_Final/NSandTA_Mass_Plots_All.pdf)

Supp Figure 14: generated by the code 2BPCode_final.R and finalized in Illustrator (standardizing font sizes, adding patient numbers and pathogen ID, adding antibiotic data)
* [NSandTA_NormalizedMass_Plots_All.pdf](Outputs_Final/NSandTA_NormalizedMass_Plots_All.pdf)

Supp Figure 15: generated by the code 2BP_AMRCode.Rmd and finalized in Illustrator (standardizing font sizes)
* A: [SuppFig - betadiv 2BP 2BP - NS vs TA.pdf](<Outputs_Final/SuppFig - betadiv 2BP - NS vs TA.pdf>)
* A: [SuppFig - betadiv 2BP NoBP - NS vs TA.pdf](<Outputs_Final/SuppFig - betadiv NoBP - NS vs TA.pdf>)

Supp Figure 16: generated by the code 2BP_AMRCode.Rmd and finalized in Illustrator (standardizing font sizes)
* A: [AMR_plot_GNR.pdf](Outputs_Final/AMR_plot_GNR.pdf)
* B: [AMR_plot_GPC.pdf](Outputs_Final/AMR_plot_GPC.pdf)
* C: [AMR_plot_both.pdf](Outputs_Final/AMR_plot_both.pdf)

Supp Figure 17: generated by the code code 2BP_AMRCode.Rmd and finalized in Illustrator (standardizing font sizes)
* [Pseudomonas_v2.pdf](Outputs_Final/Pseudomonas_v2.pdf)

## Supplementary Tables

Table S1: Clinical data, gathered from clinical record and deidentified.

Table S2: AMR/phenotype data generated in Excel from AMR data per patient (Supp Figure 16) and antimicrobial susceptibilities (Supp Data 1)


## Data files

Supp Data 1: Antimicrobial susceptibilities from standard microbiological test, gathered from clinical record and deidentified.

Supp Data 2: generated by the code 2BP_microbial_pathways_analysis.Rmd, saved as SuppData2_MicrobialPathways_MaAsLin2_2BP-vs-noBP.

Supp Data 3: generated by the code 2BP_host_gene_analysis.R, saved as SuppData3_DE_2BP-vs-noBP_covariate-viralload.

Supp Data 4: generated by the code 2BP_host_gene_analysis.R, saved as SuppData4_gsea-hallmark_2BP-vs-noBP_covariate-viralload.csv.

Supp Data 5: generated by the code 2BP_host_gene_analysis.R, saved as SuppData5_DE_2BP-vs-noBP_steroid_covariate-viralload.csv.

Supp Data 6: generated by the code 2BP_host_gene_analysis.R, saved as SuppData6_gsea-hallmark_2BP-vs-noBP_steroid_covariate-viralload.csv.

Supp Data 7: generated by the code 2BP_host_gene_analysis.R, saved as SuppData7_gsea-hallmark_2BP_steroid-days.csv.

Supp Data 8: generated by the code 2BP_host_gene_analysis.R, saved as SuppData8_gsea-hallmark_noBP_steroid-days.csv.

Supp Data 9: generated by the code 2BP_host_gene_analysis.R, saved as SuppData9_DE_2BP_mass.csv.

Supp Data 10: generated by the code 2BP_host_gene_analysis.R, saved as SuppData10_DE_noBP_mass.csv.

Supp Data 11: generated by the code 2BP_host_gene_analysis.R, saved as SuppData11_gsea-hallmark_2BP_mass.csv.

Supp Data 12: generated by the code 2BP_host_gene_analysis.R, saved as SuppData12_gsea-hallmark_noBP_mass.csv.

Supp Data 13: generated by the code 2BP_host_gene_analysis.R, saved as SuppData13_gsea-hallmark_2BP_mass_steroid.csv.

Supp Data 14: generated by the code 2BP_host_gene_analysis.R, saved as SuppData14_gsea-hallmark_noBP_mass_steroid.csv.


## Required hardware and software

The codes were run on a Mac laptop. The required R packages (see below for more details) can be installed with the commands `install.package()` and `BiocManager::install()`. Each package can take up to 1 minute to install. Please refer to each package's website for more information on the installation.

### 2BPCode_final.R

```
R version 4.3.0 (2023-04-21)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Ventura 13.3.1

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Los_Angeles
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] vegan_2.6-4     lattice_0.22-5  permute_0.9-7   lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1  
 [7] purrr_1.0.2     readr_2.1.5     tidyr_1.3.0     tibble_3.2.1    ggplot2_3.4.4   tidyverse_2.0.0
[13] dplyr_1.1.4     limma_3.56.2   

loaded via a namespace (and not attached):
 [1] utf8_1.2.4        generics_0.1.3    stringi_1.8.3     lme4_1.1-35.1     hms_1.1.3        
 [6] digest_0.6.34     magrittr_2.0.3    timechange_0.3.0  evaluate_0.23     grid_4.3.0       
[11] fastmap_1.1.1     Matrix_1.6-5      mgcv_1.9-1        fansi_1.0.6       scales_1.3.0     
[16] cli_3.6.2         rlang_1.1.3       munsell_0.5.0     splines_4.3.0     withr_3.0.0      
[21] yaml_2.3.8        parallel_4.3.0    tools_4.3.0       tzdb_0.4.0        nloptr_2.0.3     
[26] minqa_1.2.6       colorspace_2.1-0  boot_1.3-28.1     vctrs_0.6.5       R6_2.5.1         
[31] lifecycle_1.0.4   MASS_7.3-60.0.1   cluster_2.1.6     pkgconfig_2.0.3   pillar_1.9.0     
[36] gtable_0.3.4      glue_1.7.0        Rcpp_1.0.12       xfun_0.41         tidyselect_1.2.0 
[41] rstudioapi_0.15.0 knitr_1.45        farver_2.1.1      htmltools_0.5.7   nlme_3.1-164     
[46] labeling_0.4.3    rmarkdown_2.25    compiler_4.3.0   
```

### 2BP_AMRCode.Rmd

```
R version 4.2.1 (2022-06-23)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Ventura 13.2

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] rstatix_0.7.2               lubridate_1.9.3             forcats_1.0.0               stringr_1.5.1              
 [5] dplyr_1.1.4                 purrr_1.0.2                 readr_2.1.4                 tidyr_1.3.0                
 [9] tibble_3.2.1                ggplot2_3.5.1               tidyverse_2.0.0             DESeq2_1.38.3              
[13] SummarizedExperiment_1.28.0 Biobase_2.58.0              MatrixGenerics_1.10.0       matrixStats_1.2.0          
[17] GenomicRanges_1.50.2        GenomeInfoDb_1.34.9         IRanges_2.32.0              S4Vectors_0.36.2           
[21] BiocGenerics_0.44.0         vegan_2.6-4                 lattice_0.22-5              permute_0.9-7              

loaded via a namespace (and not attached):
 [1] nlme_3.1-164           bitops_1.0-7           bit64_4.0.5            RColorBrewer_1.1-3     httr_1.4.7            
 [6] tools_4.2.1            backports_1.4.1        utf8_1.2.4             R6_2.5.1               DBI_1.2.3             
[11] mgcv_1.9-1             colorspace_2.1-0       withr_3.0.0            tidyselect_1.2.0       bit_4.0.5             
[16] compiler_4.2.1         cli_3.6.2              DelayedArray_0.24.0    scales_1.3.0           digest_0.6.33         
[21] rmarkdown_2.18         XVector_0.38.0         pkgconfig_2.0.3        htmltools_0.5.7        fastmap_1.1.1         
[26] rlang_1.1.2            rstudioapi_0.16.0      RSQLite_2.3.4          generics_0.1.3         BiocParallel_1.32.6   
[31] car_3.1-2              RCurl_1.98-1.13        magrittr_2.0.3         GenomeInfoDbData_1.2.9 Matrix_1.6-4          
[36] Rcpp_1.0.11            munsell_0.5.1          fansi_1.0.6            abind_1.4-5            lifecycle_1.0.4       
[41] stringi_1.8.3          yaml_2.3.8             carData_3.0-5          MASS_7.3-58.1          zlibbioc_1.44.0       
[46] grid_4.2.1             blob_1.2.4             parallel_4.2.1         crayon_1.5.2           Biostrings_2.66.0     
[51] splines_4.2.1          annotate_1.76.0        hms_1.1.3              KEGGREST_1.38.0        locfit_1.5-9.8        
[56] knitr_1.41             pillar_1.9.0           geneplotter_1.76.0     codetools_0.2-20       XML_3.99-0.16         
[61] glue_1.6.2             evaluate_0.24.0        png_0.1-8              vctrs_0.6.5            tzdb_0.4.0            
[66] gtable_0.3.5           cachem_1.0.8           xfun_0.41              xtable_1.8-4           broom_1.0.6           
[71] AnnotationDbi_1.60.2   memoise_2.0.1          cluster_2.1.6          timechange_0.2.0  

```

### 2BP_host_gene_analysis.R

```
R version 4.3.2 (2023-10-31)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Sonoma 14.5

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: US/Pacific
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] fgsea_1.28.0    limma_3.58.1    patchwork_1.2.0 lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2    
 [9] readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0

loaded via a namespace (and not attached):
 [1] Matrix_1.6-5        gtable_0.3.5        babelgene_22.9      compiler_4.3.2      Rcpp_1.0.12         tidyselect_1.2.1    parallel_4.3.2     
 [8] scales_1.3.0        BiocParallel_1.36.0 statmod_1.5.0       lattice_0.22-6      R6_2.5.1            generics_0.1.3      munsell_0.5.1      
[15] pillar_1.9.0        tzdb_0.4.0          rlang_1.1.4         fastmatch_1.1-4     utf8_1.2.4          stringi_1.8.4       msigdbr_7.5.1      
[22] timechange_0.3.0    cli_3.6.2           withr_3.0.0         magrittr_2.0.3      grid_4.3.2          rstudioapi_0.16.0   hms_1.1.3          
[29] cowplot_1.1.3       lifecycle_1.0.4     vctrs_0.6.5         glue_1.7.0          data.table_1.15.4   codetools_0.2-20    fansi_1.0.6        
[36] colorspace_2.1-0    tools_4.3.2         pkgconfig_2.0.3    
```

### 2BP_microbial_pathways_analysis.Rmd

```
R version 4.3.2 (2023-10-31)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Sonoma 14.5

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: US/Pacific
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] Maaslin2_1.16.0 patchwork_1.2.0 lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2     readr_2.1.5    
 [9] tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0

loaded via a namespace (and not attached):
 [1] gtable_0.3.5      compiler_4.3.2    tidyselect_1.2.1  optparse_1.7.5    scales_1.3.0      R6_2.5.1          generics_0.1.3   
 [8] pcaPP_2.0-4       biglm_0.9-2.1     robustbase_0.99-2 knitr_1.47        munsell_0.5.1     DBI_1.2.3         pillar_1.9.0     
[15] tzdb_0.4.0        rlang_1.1.4       getopt_1.20.4     utf8_1.2.4        stringi_1.8.4     xfun_0.44         timechange_0.3.0 
[22] cli_3.6.2         withr_3.0.0       magrittr_2.0.3    grid_4.3.2        mvtnorm_1.2-5     rstudioapi_0.16.0 hms_1.1.3        
[29] lifecycle_1.0.4   DEoptimR_1.1-3    vctrs_0.6.5       glue_1.7.0        data.table_1.15.4 fansi_1.0.6       colorspace_2.1-0 
[36] tools_4.3.2       pkgconfig_2.0.3  
```
