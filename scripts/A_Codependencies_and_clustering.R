# Manuscript:   A reference map of CDK-cyclin dependencies across human cancers
# Project:      Codependencies
# Part:         A. Codependencies and clustering
# Author:       Agustín Sánchez-Belmonte
# Date:         03/03/23
# Last update:  01/12/25


# ---------------- INDEX ---------------- 

# 0. SET UP
# 1. LOAD DATA
#   a. Metadata
#   b. Expression
#   c. Dependency
#   d. Genomic landscape
#   e. CNV
#   f. Mutations
#   g. Pathways
#   h. Merge 
# 2. DATA EXPLORATION
# 3. CORRELATION ANALYSIS
#   a. Lineage
#   b. Primary
#   c. Plotting specific correlations
# 4. CLUSTERING
#   a. Determine n of clusters
#   b. Clustering CDKs and Cyclins
#   c. Most variable genes
#   d. Determine n of clusters selected genes
#   e. Clustering selected genes
# 5. MUTATION ANALYSIS
#   a. TMB
#   b. Plotting CIN, Aneuploidy, LoH, MSI
# 6. CLUSTER CHARACTERIZATION
#   a. Proportions by clusters
#   b. Clustering

#############################################################################################################
# 0. SET UP
#############################################################################################################

# ---------------- PATHS ---------------- 

setwd("/CODEPENDENCIES_CDKS_CYCLINS/")
PLOTS <- 'CODEPENDENCIES_CDKS_CYCLINS/25Q3_RESULTS/'
Desktop <- 'Desktop/'

# ---------------- PACKAGES ---------------- 

library(corrplot)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(readxl)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(NMF)
library(RColorBrewer)
library(pheatmap)
library(factoextra)
library(dplyr)
library(stringr)
library(limma)
library(UpSetR)
library(fmsb)
library(scales)
library(tibble)

# ---------------- FUNCTIONS ---------------- 

# Function to create a table with CDKs and Cyclins genes correlation and pvalue by tumour type
correlation_asb <- function(df, 
                               mode = c("Pancancer", "Primary", "Lineage","Cluster"),
                               type = NULL,
                               exc = NULL,
                               path = NULL,
                               plot = TRUE,
                               save = TRUE) {
  
  # ---- Packages ----
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(RColorBrewer)
  library(cowplot)
  
  mode <- match.arg(mode)
  
  # --- Validations ---
  if (mode != "Pancancer" && is.null(type)) {
    stop("Must specify 'type' when is Primary o Lineage.")
  }
  
  if (!is.null(exc) && any(!exc %in% colnames(df))) {
    stop("Some columns in 'exc' are not in the dataframe.")
  }
  
  # --- Mode selection ---
  data_type <- switch(
    mode,
    "Pancancer" = df,
    "Primary"   = df[df$OncotreePrimaryDisease == type, ],
    "Lineage"   = df[df$OncotreeLineage == type, ],
    "Cluster"   = df[df$cluster == type, ],
    stop("mode must be 'Pancancer', 'Primary','Lineage' or 'Cluster'")
  )
  t_type <- if (mode == "Pancancer") {"Pancancer"} else {type}
  # --- Column filtering ---
  numeric_df <- data_type %>% 
    select(where(is.numeric), -all_of(exc))
  
  # --- Corr and pval matrix ---
  corr_mat <- cor(numeric_df, use = "pairwise.complete.obs")
  testRes  <- cor.mtest(numeric_df, conf.level = 0.95)
  p_mat    <- testRes$p
  
  # --- Long format ---
  corr_long <- corr_mat %>%
    as.data.frame() %>%
    rownames_to_column("Partner1") %>%
    pivot_longer(-Partner1, names_to = "Partner2", values_to = "Corr")
  
  pval_long <- p_mat %>%
    as.data.frame() %>%
    rownames_to_column("Partner1") %>%
    pivot_longer(-Partner1, names_to = "Partner2", values_to = "Pvalue")
  
  # --- Merge ---
  res <- corr_long %>%
    left_join(pval_long, by = c("Partner1","Partner2")) %>%
    mutate(
      Padj = p.adjust(Pvalue, method = "bonferroni"),
      Pair = paste0(Partner1, "_", Partner2)
    )
  
  # --- Type of tumor ---
  res <- res %>% mutate(Tumor_type = t_type)
  
  # --- Mean dependency ---
  dep_means <- numeric_df %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    pivot_longer(everything(), names_to = "Partner", values_to = "Dep")
  
  # --- Merge dependencies ---
  res <- res %>%
    left_join(dep_means, by = c("Partner1" = "Partner")) %>%
    rename(Dep_P1 = Dep) %>%
    left_join(dep_means, by = c("Partner2" = "Partner")) %>%
    rename(Dep_P2 = Dep)
  
  # --- Round ---
  res <- res %>%
    mutate(
      Dep_P1 = round(Dep_P1,3),
      Dep_P2 = round(Dep_P2,3)
    )
  
  # --- Score and mix dependencies---
  res <- res %>%
    mutate(
      Score = Corr * -log10(Pvalue),
      Deps_mix = Dep_P1 + Dep_P2,
      Score = round(Score,3)
    )
  
  
  
  # --- Correlation Plot ---
  if (plot == TRUE) {
    
    if (is.null(path)) {stop("You must provide a 'path' to save the plot.")}
    
    res$Pval <- ifelse(res$Padj < 0.05 , '< 0.05','No sign')
    res$Pval <- factor(res$Pval,levels = c('< 0.05','No sign'))
    ALL <- c('CDK1','CDK2','CDK3','CDK4','CDK5','CDK6','CDK7','CDK8','CDK9','CDK10','CDK12','CDK13','CDK14','CDK15','CDK16','CDK17','CDK18','CDK19','CDK20',"CCNA1", "CCNA2", "CCNB1", "CCNB2", "CCNB3", "CCNC", "CCND1", "CCND2", "CCND3", "CCNE1", "CCNE2", "CCNF","CCNG1", "CCNG2", "CCNH", "CCNI", "CCNI2", "CCNJ", "CCNJL", "CCNK", "CCNL1", "CCNL2", "CCNO", "CCNP","CCNT1", "CCNT2", "CCNY", "CCNYL1","CNPPD1", "CNTD1", "CABLES1", "CABLES2","CDK5R1", "CDK5R2","SPDYA", "SPDYE1","SPDYE4")#  "CCNQ", "SPDYE3"
    
    res$Partner1 <- factor(res$Partner1,
                                levels = ALL)
    res$Partner2 <- factor(res$Partner2,
                           levels = rev(ALL))
    
    max <- max(res$Deps_mix)
    min <- min(res$Deps_mix)
    
    # Whole correlation plot
    png(filename = file.path(path, paste0(t_type,'_','whole_corr.png')), width = 3400, height = 3200, res = 305)
    
    p <- ggplot(res, aes(Partner1, Partner2)) +
      geom_point(aes(colour = Pval,
                     size=Deps_mix*-1,          
                     shape = Pval,
                     fill = Corr)) +
      scale_size('Deps_mix',limits = c(max*-1,min*-1)) +
      theme_bw() +
      scale_color_manual(values=c("black","#E0E0E0")) +    
      scale_fill_gradientn(colours = rev(brewer.pal(name = 'RdBu', n=11)),
                           limits = c(-1,1),
                           oob = scales::squish,
                           name = 'Correlation') +
      scale_x_discrete(position = "top") +
      cowplot::theme_cowplot() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_shape_manual(values=c(21,21))
  
    print(p)
    dev.off()
    
    # Codependencies correlation plot

    CDKS <- c('CDK1','CDK2','CDK3','CDK4','CDK5','CDK6','CDK7','CDK8','CDK9','CDK10','CDK12','CDK13','CDK14','CDK15','CDK16','CDK17','CDK18','CDK19','CDK20') # ['CDK11', "CDK11A", "CDK11B"] not in index
    CYCLINS <- c("CCNA1", "CCNA2", "CCNB1", "CCNB2", "CCNB3", "CCNC", "CCND1", "CCND2", "CCND3", "CCNE1", "CCNE2", "CCNF","CCNG1", "CCNG2", "CCNH", "CCNI", "CCNI2", "CCNJ", "CCNJL", "CCNK", "CCNL1", "CCNL2", "CCNO", "CCNP","CCNT1", "CCNT2", "CCNY", "CCNYL1","CNPPD1", "CNTD1", "CABLES1", "CABLES2","CDK5R1", "CDK5R2","SPDYA", "SPDYE1","SPDYE4") ## [''CCNYL2'] is a
    
    res_filt <- res[res$Partner2 %in% CDKS,]
    res_filt <- res_filt[res_filt$Partner1 %in% CYCLINS,]
    
    res_filt$Partner1 <- factor(res_filt$Partner1,
                          levels = CYCLINS)
    res_filt$Partner2 <- factor(res_filt$Partner2,
                          levels = rev(CDKS))
    
    png(filename = file.path(path, paste0(t_type,'_','codep_corr.png')), width = 2500, height = 1500, res = 305)
    
    p <- ggplot(res_filt, aes(Partner1, Partner2)) +
      geom_point(aes(colour = Pval,
                     size=Deps_mix*-1,          # Corregido: size no negativo
                     shape = Pval,
                     fill = Corr)) +
      scale_size('Deps_mix',limits = c(max*-1,min*-1)) +
      theme_bw() +
      scale_color_manual(values = c("< 0.05" = "black",  
                                    "No sign" = "#E0E0E0")) +    
      scale_fill_gradientn(colours = rev(brewer.pal(name = 'RdBu', n=11)),
                           limits = c(-1,1),
                           oob = scales::squish,
                           name = 'Correlation') +
      scale_x_discrete(position = "top") +
      cowplot::theme_cowplot() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_shape_manual(values=c("< 0.05" = 21,  
                                  "No sign" = 21))
    print(p)
    dev.off()
    
  
  }
  
  
  # --- Remove self-corr ---
  res <- res %>% filter(Partner1 != Partner2)
  
  # --- Remove duplicates ---
  res <- res %>%
    mutate(pair_id = paste(
      pmin(as.character(Partner1), as.character(Partner2)),
      pmax(as.character(Partner1), as.character(Partner2)),
      sep = "_")) %>%
    distinct(pair_id, .keep_all = TRUE) %>%
    select(-pair_id)
  
  # --- Ranking ---
  res <- res %>%
    mutate(Rank = min_rank(desc(Score)))
  
  # --- Save table corr ---
  if (save == TRUE) {
    write.table(res, file.path(path, paste0(t_type,'_','table_corr.tsv')), sep = "\t", row.names = FALSE, quote = FALSE)
  }
  # --- Volcan Plot ---
  if (plot == TRUE) {
    
    
    # Volcan plot
    res$delabel <- ifelse((res$Corr > 0.2 & res$Padj < 0.05), 'Positive', ifelse((res$Corr < -0.2 & res$Padj < 0.05), 'Negative','NA'))
    res$gene_lab <- ifelse((abs(res$Corr) > 0.2 & res$Padj < 0.05), res$Pair, '')
    res$Target <-  ifelse(( abs(res$Corr) > 0.2 & res$Padj < 0.05), 'Yes', 'No')
    
    max <- max(res$Corr)
    min <- min(res$Corr)
    val <- ifelse(abs(min) > abs(max), abs(min),abs(max))
    
    png(file.path(path, paste0(t_type,'_','volcan.png')),width=1200, height=1200, res = 300)
    p <- ggplot(data=res, aes(x=Corr, y=-log10(Padj), label=gene_lab)) + 
      geom_point(aes(shape=delabel,fill=delabel, size=delabel,color=Target)) +
      theme_minimal() +
      scale_shape_manual(values=c(20,21,21))+
      #geom_text_repel(box.padding = 0.2,  
      #                segment.size = 0.2,  
      #                segment.color = "black",
      #                max.overlaps = 100,
      #                size = 2) + 
      scale_fill_manual(values = c("Positive" = "#FF6666",  
                                   "Negative" = "#1874CD"), 
                        na.value = "#E0E0E0" ) +
      scale_color_manual(values=c("#E0E0E0","black"))+
      scale_size_manual(values=c(0.8,1.8,1.8))+
      geom_vline(xintercept=c(0),col="#C0C0C0", size=0.2) +
      xlim(c(-val,val)) +
      xlab("FC") + 
      ylab("-log10(Padj)") +
      theme(panel.grid = element_blank(),
            plot.margin = margin(3, 3, 3, 3, "pt"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(color = "black"))+
      guides(color = FALSE, shape = FALSE, fill=FALSE, size=FALSE) 
    print(p)
    dev.off()
  }
  
  # --- Clean table ---
  res <- res %>% select(-Pval, -delabel, -gene_lab, -Target)
  
  # --- Return corr table ---
  return(res)
}

source('/Users/asanchezb/Desktop/PROJECTS/CDK_PREDICT_CDK4_6_BRCA/OLD_ANALYSIS/R/theme.R')
figure_font_size=12

#############################################################################################################
# 1. LOAD DATA
#############################################################################################################

# Data were downloaded from DepMap 25Q3 (DepMap, 2025).

# ---------------- METADATA ---------------- 

metadata <- read.csv('25Q3_data/Model.csv', check.names = F)

# ---------------- EXPRESSION ---------------- 

expression_data <- read.csv('25Q3_data/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv', check.names = F)

names <- c()
for(i in colnames(expression_data)){
  new <- strsplit(i, " ")[[1]][1]
  names <- c(names,new)
}
expression_data <- expression_data[expression_data$IsDefaultEntryForModel == 'Yes',]
colnames(expression_data) <- names

names(expression_data)[1] <- 'Test'
expression_data <- expression_data %>%
  select(-Test, -SequencingID, -IsDefaultEntryForMC, -IsDefaultEntryForModel, -ModelConditionID)

cell_names <- expression_data$ModelID
expression_data <- expression_data[,-c(1)]
row.names(expression_data) <- cell_names

# ---------------- DEPENDENCY ---------------- 

dependency_data <- read.csv('25Q3_data/CRISPRGeneEffect.csv', check.names = F)

names <- c()
for(i in colnames(dependency_data)){
  new <- strsplit(i, " ")[[1]][1]
  names <- c(names,new)
}

colnames(dependency_data) <- names
names(dependency_data)[1] <- 'ModelID'

# ---------------- GENOMIC LANDSCAPE ---------------- 

genomic_profile <- read.csv('25Q3_data/OmicsGlobalSignatures.csv', check.names = F)
genomic_profile <- genomic_profile[genomic_profile$IsDefaultEntryForModel == 'Yes',]
genomic_profile <-  genomic_profile[,c("ModelID","MSIScore","Ploidy","WGD","LoHFraction","CIN","Aneuploidy")]

metadata <- merge(metadata, genomic_profile, by = "ModelID", all.x=TRUE)


# ---------------- CNV ---------------- 

cnv_data <- read.csv('25Q3_data/OmicsCNGeneWGS.csv', check.names = F)

names <- c()
for(i in colnames(cnv_data)){
  new <- strsplit(i, " ")[[1]][1]
  names <- c(names,new)
}

colnames(cnv_data) <- names
cnv_data <- cnv_data[cnv_data$IsDefaultEntryForModel == 'Yes',]
cnv_data <- cnv_data[, -c(1, 2, 4, 5, 6)]

# Add CNV to metadata

cnv_relevants <- merge(cnv_data,metadata,by.x='ModelID',by.y='ModelID')

gene_cols <- colnames(cnv_data)[-1]   # quitando ModelID si era la primera

cnv_long <- cnv_relevants %>%
  pivot_longer(
    cols = all_of(gene_cols),
    names_to = "Gene",
    values_to = "CNV"
  )

cnv_long <- cnv_long[,c('ModelID','Ploidy','Gene','CNV')]

cnv_long$CNV_rel <- log2(cnv_long$CNV / cnv_long$Ploidy)
cnv_long$CNV_rel_2 <- log2(cnv_long$CNV / cnv_long$Ploidy + 1)



# add ampl/loss

genes_amp <- c("CCND1", "CCND2","CCND3", "CDK4", "CDK6", "CDK2","CCNE1", "CCNE2", "MYC",'MYCN',"MYCL",'ERBB2')

for (g in genes_amp) {
  
  metadata[[paste0(g, "_amp")]] <- ifelse(
    metadata$ModelID %in%
      cnv_long$ModelID[
        cnv_long$Gene == g &
          cnv_long$CNV_rel_2 > 0.55
        ],
    "amp",
    "WT"
  )
}

genes_loss <- c("CDKN2A", "CDKN2B","CDKN2C", "CDKN2D","RB1",'PTEN')

for (g in genes_loss) {
  
  metadata[[paste0(g, "_loss")]] <- ifelse(
    metadata$ModelID %in%
      cnv_long$ModelID[
        cnv_long$Gene == g &
          cnv_long$CNV_rel_2 < 0.3
        ],
    "loss",
    "WT"
  )
}


# ---------------- MUTATIONS ---------------- 

mutations_data <- read.csv('25Q3_data/OmicsSomaticMutations.csv', check.names = F, row.names = 1)

# Relevants types of mutations
relevant_mut <- c(
  "missense_variant",
  "frameshift_variant",
  "stop_gained",
  "stop_lost",
  "start_lost",
  "inframe_deletion",
  "inframe_insertion",
  "protein_altering_variant",
  "splice_donor_variant",
  "splice_acceptor_variant",
  "splice_region_variant"
)


# Filter relevant mutations
mutaciones_relevants <- mutations_data %>%
  filter(str_detect(VariantInfo, paste(relevant_mut, collapse = "|")))

# Driver genes
genes <- c("RB1","TP53","PIK3CA","AKT1","CDKN2A","CDKN2B","CDKN2C","CDKN2D","PTEN","EGFR","KRAS","RBL1","RBL2",
           "CCND1","CCND2","CCND3","CDK4","CDK6","CTNNB1","APC","AXIN1","AXIN2","ALK","MET","ROS1","BRAF","RET",'CDK12')

# Status genes
for(g in genes){
  modelos_mutados <- mutaciones_relevants$ModelID[mutaciones_relevants$HugoSymbol == g]
  metadata[[paste0(g, "_status")]] <- ifelse(metadata$ModelID %in% modelos_mutados, "mut", "WT")
}

# ---------------- PATHWAYS ---------------- 
metadata$PI3K_pathway <- ifelse(rowSums(metadata[, c("PIK3CA_status","AKT1_status","PTEN_status")] == "mut") > 0 |metadata$PTEN_loss == "loss","mut","WT")
metadata$RBL1_2_status  <- ifelse(rowSums(metadata[, c("RBL1_status","RBL2_status")] == "mut") > 0, "mut", "WT")
metadata$RB1_act <- ifelse((metadata$RB1_status == 'mut')|(metadata$RB1_loss == 'loss'),'mut','WT')
metadata$WNT_pathway    <- ifelse(rowSums(metadata[, c("CTNNB1_status","APC_status","AXIN1_status","AXIN2_status")] == "mut") > 0, "mut", "WT")

metadata$Rb_pathway <- ifelse(
  (metadata$CCND1_amp == 'amp') |
    (metadata$CCND2_amp == 'amp') |
    (metadata$CCND3_amp == 'amp') |
    (metadata$CDK4_amp  == 'amp') |
    (metadata$CDK6_amp  == 'amp') |
    (metadata$CCND1_status == 'mut') |
    (metadata$CCND2_status == 'mut') |
    (metadata$CCND3_status == 'mut') |
    (metadata$CDK4_status  == 'mut') |
    (metadata$CDK6_status  == 'mut') |
    (metadata$CDKN2A_loss == 'loss') |
    (metadata$CDKN2B_loss == 'loss') |
    (metadata$CDKN2C_loss == 'loss') |
    (metadata$CDKN2D_loss == 'loss') |
    (metadata$CDKN2A_status == 'mut') |
    (metadata$CDKN2B_status == 'mut') |
    (metadata$CDKN2C_status == 'mut') |
    (metadata$CDKN2D_status == 'mut'),
  'mut',
  'WT'
)

# Group by RB
metadata$Group <- ifelse(metadata$Rb_pathway == 'mut', 'RBupstream','WT')
metadata$Group <- ifelse(metadata$RB1_act == 'mut', 'RBmut',metadata$Group)

table(metadata$Group) # RBmut: 384 / RBupstream: 762 / WT: 986

# Metadata with mutations, genomic profile and cnv
#write.table(metadata, '/Users/asanchezb/Desktop/metadata_280226.tsv', sep = '\t', row.names = F)

# ---------------- Clean metadata ---------------- 

metadata <- metadata %>%
  mutate(Subtype = case_when(
    ModelSubtypeFeatures %in% c('basal TNBC', 'basal_A TNBC', 'basal_B TNBC','luminal TNBC', 'TNBC') ~ 'TNBC',
    ModelSubtypeFeatures %in% c('ER+', 'luminal ER, PR+', 'PR+') ~ 'ER+',
    ModelSubtypeFeatures %in% c('HER2+','basal_A HER2+') ~ 'HER2+',
    ModelSubtypeFeatures %in% c('luminal ER+, HER2+', 'luminal ER+, PR+, HER2+','ER+, HER2+') ~ 'ER+ HER2+',
    ModelSubtypeFeatures == '' ~ NA_character_,
    TRUE ~ ModelSubtypeFeatures  # cualquier otro valor se mantiene
  ))


# ---------------- CDKs & Cyclins ---------------- 

CDKS <- c('ModelID','CDK1','CDK2','CDK3','CDK4','CDK5','CDK6','CDK7','CDK8','CDK9','CDK10','CDK12','CDK13','CDK14','CDK15','CDK16','CDK17','CDK18','CDK19','CDK20') # ['CDK11', "CDK11A", "CDK11B"] not in index
CYCLINS <- c("ModelID","CCNA1", "CCNA2", "CCNB1", "CCNB2", "CCNB3", "CCNC", "CCND1", "CCND2", "CCND3", "CCNE1", "CCNE2", "CCNF","CCNG1", "CCNG2", "CCNH", "CCNI", "CCNI2", "CCNJ", "CCNJL", "CCNK", "CCNL1", "CCNL2", "CCNO", "CCNP","CCNT1", "CCNT2", "CCNY", "CCNYL1","CNPPD1", "CCNQ", "CNTD1", "CABLES1", "CABLES2","CDK5R1", "CDK5R2","SPDYA", "SPDYE1","SPDYE3","SPDYE4") ## [''CCNYL2'] is a pseudogene; "CCNYL1B" not in index  # as in Quandt et al., 2019
ALL <- c('ModelID','CDK1','CDK2','CDK3','CDK4','CDK5','CDK6','CDK7','CDK8','CDK9','CDK10','CDK12','CDK13','CDK14','CDK15','CDK16','CDK17','CDK18','CDK19','CDK20',"CCNA1", "CCNA2", "CCNB1", "CCNB2", "CCNB3", "CCNC", "CCND1", "CCND2", "CCND3", "CCNE1", "CCNE2", "CCNF","CCNG1", "CCNG2", "CCNH", "CCNI", "CCNI2", "CCNJ", "CCNJL", "CCNK", "CCNL1", "CCNL2", "CCNO", "CCNP","CCNT1", "CCNT2", "CCNY", "CCNYL1","CNPPD1", "CNTD1", "CABLES1", "CABLES2","CDK5R1", "CDK5R2","SPDYA", "SPDYE1","SPDYE4")#  "CCNQ", "SPDYE3"




#############################################################################################################
# 2. DATA EXPLORATION
#############################################################################################################


dep_CDKs <- dependency_data[,ALL]

#write.table(dep_CDKs, '25Q3_RESULTS/CDKs_Cyclins_dependencies.csv', sep = ',', row.names = F)

# Merge
select <- c('ModelID','StrippedCellLineName','OncotreeLineage','OncotreePrimaryDisease','ModelSubtypeFeatures')
merged <- merge(metadata[,select],dep_CDKs, by.x = 'ModelID', by.y = 'ModelID')


# ---------------- FIG 1A ---------------- 


genes <- setdiff(colnames(merged),
                 c('ModelID','StrippedCellLineName',
                   'OncotreeLineage','OncotreePrimaryDisease',
                   'ModelSubtypeFeatures'))

mean_df <- merged %>%
  mutate(OncotreeLineage = ifelse(OncotreePrimaryDisease == "Non-Cancerous",
                                  "Non-Cancerous",
                                  OncotreeLineage)) %>%
  pivot_longer(all_of(genes),
               names_to = "Gene",
               values_to = "Dependency") %>%
  group_by(OncotreeLineage, Gene) %>%
  summarise(mean_dep = mean(Dependency, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(Gene = factor(Gene, levels = genes))

# Media Non-Cancerous dentro del mismo objeto
non_cancer_df <- mean_df %>%
  filter(OncotreeLineage == "Non-Cancerous") %>%
  select(Gene, non_cancer_mean = mean_dep)

png('/Users/asanchezb/Desktop/Fig1A.png', width = 2800, height = 900, res = 250)

ggplot(mean_df, aes(Gene, mean_dep)) +
  geom_point(aes(fill = mean_dep),
             shape = 21,
             size = 2.5,
             stroke = 0.5,
             color = "black",
             alpha = 0.95,
             position = position_jitter(width = 0.2)) +
  
  geom_hline(yintercept = c(-1, -0.5, 0),
             linetype = "dashed",
             color = c("#A90C38", "#D59B9B", "#C8BCDF"),
             linewidth = 0.6) +
  
  geom_segment(data = non_cancer_df,
               aes(x = as.numeric(Gene) - 0.3,
                   xend = as.numeric(Gene) + 0.3,
                   y = non_cancer_mean,
                   yend = non_cancer_mean),
               inherit.aes = FALSE,
               color = "#53B6F9",
               linewidth = 0.5) +
  
  scale_fill_gradientn(
    colours = c("#A90C38", "white", "#0B0B8C"),
    values = scales::rescale(c(-3, 0, 1)),
    limits = c(-3, 1),
    name = "Gene median"
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_text(size = 12)
  ) +
  
  labs(
    x = "Gene",
    y = "MeanDep Lineage"
  )

dev.off()


# ---------------- FIG S1A ---------------- 

class <- read_excel('/Users/asanchezb/Downloads/Primarydisease_classified.xlsx', sheet = 1) |>
  setNames(c('OncotreePrimaryDisease','Class'))

dep_df <- merged %>%
  left_join(class, by = "OncotreePrimaryDisease")

genes <- setdiff(colnames(dep_df),
                 c('ModelID','StrippedCellLineName',
                   'OncotreeLineage','OncotreePrimaryDisease',
                   'ModelSubtypeFeatures','Class'))

# Matriz directamente en un pipeline
heat_mat <- dep_df %>%
  pivot_longer(all_of(genes),
               names_to = "Gene",
               values_to = "Dependency") %>%
  group_by(Class, Gene) %>%
  summarise(mean_dep = mean(Dependency, na.rm = TRUE),
            #sd   = sd(Dependency, na.rm = TRUE),
            # n        = sum(!is.na(Dependency)),
            .groups = "drop") %>%
  pivot_wider(names_from = Gene,
              values_from = mean_dep) %>%
  column_to_rownames("Class") %>%
  as.matrix()

# Orden opcional por Class (si quieres mantenerlo)
heat_mat <- heat_mat[order(rownames(heat_mat)), genes]

# Colores
class_colors <- c("CNS/Neuro"  = '#cab2d6',
                  "Epithelial" = '#fb9a99',
                  "Germ/Others"  = '#a6cee3',
                  "Hematopoietic" = '#fdd0a2',
                  "Mesenchymal"  = '#8B864E',
                  "NeuroEndocrine"= '#b2df8a')

col_fun <- colorRamp2(c(-2.5, 0, 1),
                      c("#A90C38", "white", "#0B0B8C"))

png('/Users/asanchezb/Desktop/FigS1A.png', width = 3000, height = 600, res = 250)

Heatmap(heat_mat,
        name = "Dependency",
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        column_names_rot = 90,
        border = TRUE)

dev.off()

colores_paired_25 <- c(
  "#1f78b4", "#a6cee3",  # Azul / Azul pastel
  "#33a02c", "#b2df8a",  # Verde / Verde pastel
  "#e31a1c", "#fb9a99",  # Rojo / Rojo pastel
  "#ff7f00", "#fdbf6f",  # Naranja / Naranja pastel
  "#6a3d9a", "#cab2d6",  # Púrpura / Lila pastel
  "#b15928", "#fdd0a2",  # Marrón / Marrón pastel
  "#458B74", "#66CDAA",  # Azul pastel / Azul oscuro
  "#8B8878", "#CDC8B1",  # Verde pastel / Verde oscuro
  "#8B7500", "#CDAD00",  # Rojo pastel / Rojo oscuro
  "#8B3A62", "#CD6090",  # Naranja pastel / Naranja oscuro
  "#8B864E", "#CDC673",  # Lila pastel / Púrpura oscuro
  "#8B4C39", "#EE8262" ,"#B8860B"     # Marrón pastel / Marrón oscuro
)



# ---------------- FIG S1B ---------------- 

# ----------------------
# Resumen por PrimaryDisease y Gene
heat_df <- long_df %>%
  left_join(class, by = "OncotreePrimaryDisease") %>%
  mutate(Dependency = as.numeric(as.character(Dependency))) %>%  # asegura numeric
  group_by(OncotreePrimaryDisease, Class, Gene) %>%
  summarise(mean_dep = mean(Dependency, na.rm = TRUE),
            .groups = "drop")

# Crear matriz para heatmap
# Crear matriz solo con los genes
heat_mat <- heat_df %>%
  select(-Class) %>%                      # quitar la columna Class
  pivot_wider(names_from = Gene, values_from = mean_dep) %>%
  column_to_rownames("OncotreePrimaryDisease") %>%
  as.matrix()   

# Filtrar solo columnas que existen en heat_mat (evita vacío)
heat_mat <- heat_mat[, intersect(genes, colnames(heat_mat))]

# Ordenar filas por Class y alfabéticamente
row_order <- heat_df %>%
  select(OncotreePrimaryDisease, Class) %>%
  distinct() %>%
  arrange(Class, OncotreePrimaryDisease)

heat_mat <- heat_mat[intersect(row_order$OncotreePrimaryDisease, rownames(heat_mat)), ]

# Paleta de colores
col_fun <- colorRamp2(c(-2.5, 0, 1), c("#A90C38", "white", "#0B0B8C"))

# Anotación de clases usando class_colors previamente definido
row_ha <- rowAnnotation(
  Class = row_order$Class,
  col = list(Class = class_colors),
  show_annotation_name = TRUE
)

# Exportar heatmap
png('/Users/asanchezb/Desktop/FigS1B.png', width = 3000, height = 3000, res = 250)

Heatmap(heat_mat,
        name = "Dependency",
        col = col_fun,
        cluster_rows = FALSE,       # respetar orden por Class
        cluster_columns = FALSE,
        row_names_side = "left",
        column_names_rot = 90,
        border = TRUE,
        left_annotation = row_ha)

dev.off()


# ---------------- FIG 1C ---------------- 



# --- Datos ejemplo ---
radar_df <- dep_df %>%
  group_by(Class) %>%
  summarise(
    CDK4 = mean(CDK4, na.rm = TRUE),
    CDK6 = mean(CDK6, na.rm = TRUE)#,
    #CCND3 = mean(CCND3, na.rm = TRUE)
  ) %>%
  ungroup()

# Transponer: filas = CDK, columnas = clase
radar_plot_df <- t(radar_df[,2:3])
colnames(radar_plot_df) <- radar_df$Class
rownames(radar_plot_df) <- c("CDK4","CDK6")#,"CCND3")
radar_plot_df <- radar_plot_df*-1

# --- Establecer Max y Min dinámicos basados en los datos reales ---
max_vals <- rep(3, ncol(radar_plot_df))
min_vals <- rep(0, ncol(radar_plot_df))

radar_plot_df <- as.data.frame(rbind(max_vals, min_vals, radar_plot_df))

# Colores
cdk_colors <- c("#739FD3","#E08A4E","#85B15F")

# --- Dibujar radar plot ---
png('/Users/asanchezb/Desktop/test.png', width = 2000, height = 2000, res = 280)
radarchartcirc(
  radar_plot_df,
  pcol = cdk_colors,
  plwd = 2,
  plty= 1,
  seg=6,
  pfcol = alpha(cdk_colors, 0.3),
  cglcol = "grey",
  cglty = 1,
  axistype = 1,
  axislabcol = "black",
  caxislabels = seq(0, 3, length.out = 7), # 7 líneas exactas
  vlcex = 1.2
)

legend(
  "topright",
  legend = rownames(radar_plot_df)[-c(1,2)],
  col = cdk_colors,
  lty = 1,
  lwd = 2,
  bty = "n"
)
dev.off()

# ---------------- FIG 1F ---------------- 

new_df <- merge(metadata,dep_CDKs, by.x = 'ModelID', by.y = 'ModelID')
new_df <- new_df[new_df$OncotreeLineage == 'Breast',]

new_df <- new_df %>%
  mutate(Subtype = case_when(
    ModelSubtypeFeatures %in% c('basal TNBC', 'basal_A TNBC', 'basal_B TNBC','luminal TNBC', 'TNBC') ~ 'TNBC',
    ModelSubtypeFeatures %in% c('ER+', 'luminal ER, PR+', 'PR+') ~ 'ER+',
    ModelSubtypeFeatures %in% c('HER2+','basal_A HER2+') ~ 'HER2+',
    ModelSubtypeFeatures %in% c('luminal ER+, HER2+', 'luminal ER+, PR+, HER2+','ER+, HER2+') ~ 'ER+ HER2+',
    ModelSubtypeFeatures == '' ~ NA_character_,
    TRUE ~ ModelSubtypeFeatures  # cualquier otro valor se mantiene
  ))

new_df

# --- Datos ejemplo ---
radar_df <- new_df %>%
  group_by(Subtype) %>%
  summarise(
    CDK4 = mean(CDK4, na.rm = TRUE),
    CDK6 = mean(CDK6, na.rm = TRUE)
    #CCND1 = mean(CCND1, na.rm = TRUE),
    #CCND2 = mean(CCND2, na.rm = TRUE),
    #CCND3 = mean(CCND3, na.rm = TRUE)
  ) %>%
  ungroup()

# Transponer: filas = CDK, columnas = clase
radar_plot_df <- t(radar_df[,2:3])
colnames(radar_plot_df) <- radar_df$Subtype
rownames(radar_plot_df) <- c("CDK4","CDK6")
rownames(radar_plot_df) <- c("CCND1","CCND2","CCND3")
radar_plot_df <- radar_plot_df*-1

# --- Establecer Max y Min dinámicos basados en los datos reales ---
max_vals <- rep(3, ncol(radar_plot_df))
min_vals <- rep(0, ncol(radar_plot_df))

radar_plot_df <- as.data.frame(rbind(max_vals, min_vals, radar_plot_df))
radar_plot_df <- radar_plot_df[,1:4]
# Colores
#cdk_colors <- c("#CDC8B1","#EE6A50")#,"#2E8B57")

# --- Dibujar radar plot ---
png('/Users/asanchezb/Desktop/test.png', width = 2000, height = 2000, res = 280)
radarchartcirc(
  radar_plot_df,
  pcol = cdk_colors,
  plwd = 2,
  plty= 1,
  seg=6,
  pfcol = alpha(cdk_colors, 0.3),
  cglcol = "grey",
  cglty = 1,
  axistype = 1,
  axislabcol = "black",
  caxislabels = seq(0, 3, length.out = 7), # 7 líneas exactas
  vlcex = 1
)

legend(
  "topright",
  legend = rownames(radar_plot_df)[-c(1,2)],
  col = cdk_colors,
  lty = 1,
  lwd = 2,
  bty = "n"
)
dev.off()



# ---------------- FIG 1D ---------------- 

genes_plot <- c("CDK4","CDK6","CCND1","CCND2","CCND3")


mean_df <- dep_df %>%
  select(OncotreeLineage, CDK4, CDK6, CCND1, CCND2, CCND3) %>%
  pivot_longer(cols = c(CDK4, CDK6, CCND1, CCND2, CCND3),
               names_to = "Gene",
               values_to = "Dependency") %>%
  group_by(OncotreeLineage, Gene) %>%
  summarise(mean_dep = mean(Dependency, na.rm = TRUE),
            .groups = "drop")


mean_df$Gene <- factor(mean_df$Gene, levels = genes_plot)


# 3️⃣ Plot usando SOLO medias
png('/Users/asanchezb/Desktop/boxplot_lineage_means.png',
    width = 2800, height = 2000, res = 280)

ggplot(mean_df, aes(x = Gene, y = mean_dep)) +
  
  geom_boxplot(fill = "grey90",
               color = "black",
               width = 0.6,
               outlier.shape = NA) +
  
  # Puntos con borde negro y relleno por dependencia
  geom_point(aes(fill = OncotreeLineage),
             size = 3,
             shape = 21,
             alpha = 0.9,
             stroke = 0.5,
             color = "black",
             position = position_jitter(width = 0.15, height = 0)) +
  
  geom_hline(yintercept = -1,
             linetype = "dashed",
             color = "#A90C38",
             linewidth = 0.6) +
  
  geom_hline(yintercept = -0.5,
             linetype = "dashed",
             color = "orange",
             linewidth = 0.6) +
  
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "#C8BCDF",
             linewidth = 0.6) +
  
  #scale_fill_brewer(palette = "Set2", name = "Tumor Class") +
  # Paleta discreta para Class
  scale_fill_manual(values = colores_paired_25,
                    name = "Lineage") +

  
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.major.x = element_blank()
  ) +
  
  labs(
    x = "Gene",
    y = "Mean Dep per Tumor"
  )

dev.off()

mean_df <- dep_df %>%
  #select(OncotreePrimaryDisease, Class, CDK4, CDK6, CCND1, CCND2, CCND3) %>%
  pivot_longer(cols = all_of(genes),
 # pivot_longer(cols = c(CDK4, CDK6, CCND1, CCND2, CCND3),
               names_to = "Gene",
               values_to = "Dependency")# %>%
  #group_by( Gene) %>%
  #summarise(mean_dep = mean(Dependency, na.rm = TRUE),
   #         .groups = "drop")


mean_df$Gene <- factor(mean_df$Gene, levels = ALL)

png('/Users/asanchezb/Desktop/boxplot_lineage_means.png',
    width = 3500, height = 1200, res = 250)

  ggplot(mean_df, aes(x = Gene, y = Dependency)) +
  
  geom_boxplot(fill = "grey90",
               color = "black",
               width = 0.6,
               outlier.shape = NA) +
  
  # Puntos coloreados por Class
  geom_point(aes(fill = Dependency),
             size = 0.5,
             shape = 21,
             alpha = 0.5,
             stroke = 0.5,
             color = "black",
             position = position_jitter(width = 0.15, height = 0)) +
  
  geom_hline(yintercept = -1,
             linetype = "dashed",
             color = "#A90C38",
             linewidth = 0.6) +
  
  geom_hline(yintercept = -0.5,
             linetype = "dashed",
             color = "orange",
             linewidth = 0.6) +
  
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "#C8BCDF",
             linewidth = 0.6) +
  
  # Paleta discreta para Class
  #scale_fill_manual(values = class_colors,
                 #   name = "Class") +
    # Tu escala personalizada
   scale_fill_gradientn(
      colours = c("#A90C38", "white", "#0B0B8C"),
      values = scales::rescale(c(-3, 0, 1)),
      limits = c(-4.5, 1),
      name = "Gene median"
    ) +
  
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5),
    panel.grid.major.x = element_blank()
  ) +
  
  labs(
    x = "Gene",
    y = "Dependency all cell lines"
  )

dev.off()

write.table(heat_df,'/Users/asanchezb/Desktop/mean_primary_class_25Q3.tsv', sep = "\t")



#############################################################################################################
# 3. CORRELATION
#############################################################################################################

# ---------------- Pancancer ---------------- 

res <- correlation_asb(merged,
                       mode = "Pancancer",
                       type = NULL,
                       exc = select,
                       path = '25Q3_RESULTS',
                       plot = TRUE,
                       save = TRUE)

# ---------------- Lineage ---------------- 

dir.create('25Q3_RESULTS/Lineage')

merged <- merged %>%
  mutate(OncotreeLineage = gsub("/", "-", OncotreeLineage))

results_list <- list()

for(i in unique(merged$OncotreeLineage)){
  # Skip groups less than 7
  n_rows <- nrow(merged[merged$OncotreeLineage == i, ])
  if(n_rows < 7) {
    message("Skipping lineage ", i, " (less than 7 cell lines)")
    next  # salta a la siguiente iteración del for
  } else {print(i)}
  
  res <- correlation_asb(merged,
                         mode = "Cluster",
                         type = i,
                         exc = select,
                         path = '25Q3_RESULTS/Lineage',
                         plot = TRUE,
                         save = TRUE)
  results_list[[i]] <- res
  
}

# Combine results

lineage_corr <- do.call(rbind, results_list)
write.table(lineage_corr, '25Q3_RESULTS/Lineage/Lineage_table_corr.tsv', sep = "\t", row.names = FALSE, quote = FALSE)

# ---------------- Primary ---------------- 

dir.create('25Q3_RESULTS/Primary')

table(merged$OncotreePrimaryDisease)

merged <- merged %>%
  mutate(OncotreePrimaryDisease = gsub("/", "-", OncotreePrimaryDisease))

results_list <- list()

for(i in unique(merged$OncotreePrimaryDisease)){
  # Skip groups less than 7
  n_rows <- nrow(merged[merged$OncotreePrimaryDisease == i, ])
  if(n_rows < 7) {
    message("Skipping tumor type ", i, " (less than 7 cell lines)")
    next  # salta a la siguiente iteración del for
  } else {print(i)}
  
  res <- correlation_asb(merged,
                         mode = "Primary",
                         type = i,
                         exc = select,
                         path = '25Q3_RESULTS/Primary',
                         plot = TRUE,
                         save = TRUE)
  results_list[[i]] <- res
  
}

# Combine results

primary_corr <- do.call(rbind, results_list)
write.table(primary_corr, '25Q3_RESULTS/Primary/Primary_table_corr.tsv', sep = "\t", row.names = FALSE, quote = FALSE)

# ---------------- Clusters ---------------- 

dir.create('25Q3_RESULTS/Cluster')
results_list <- list()

merged_clu <- merge(merged, met_filt[,c('ModelID','cluster')], by='ModelID')
select2 <- c("ModelID","StrippedCellLineName","OncotreeLineage","OncotreePrimaryDisease","ModelSubtypeFeatures","cluster")

for(i in unique(merged_clu$cluster)){
  # Skip groups less than 7
  n_rows <- nrow(merged_clu[merged_clu$cluster == i, ])
  if(n_rows < 7) {
    message("Skipping cluster ", i, " (less than 7 cell lines)")
    next  # salta a la siguiente iteración del for
  } else {print(i)}
  
  res <- correlation_asb(merged_clu,
                         mode = "Cluster",
                         type = i,
                         exc = select2,
                         path = '25Q3_RESULTS/Cluster',
                         plot = TRUE,
                         save = TRUE)
  results_list[[i]] <- res
  
}

# Combine results

primary_corr <- do.call(rbind, results_list)
write.table(primary_corr, '25Q3_RESULTS/Cluster/Cluster_table_corr.tsv', sep = "\t", row.names = FALSE, quote = FALSE)


# ---------------- Plotting specific correlations ---------------- 

# To select lineage
merged_sel <- merged[merged$OncotreeLineage =='Ovary/Fallopian Tube',]
merged_sel <- merged[merged$OncotreeLineage =='Breast',]
merged_sel <- merged[merged$OncotreeLineage =='Myeloid',]

png('25Q3_RESULTS/corr_CDK2_CCNE1_ovarian.png', width = 1000, height = 1000, res = 250)
ggscatter(merged_sel, x = "CDK2", y = "CCNE1", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = FALSE, cor.method = "pearson",
          xlab = "CDK2", ylab = "CCNE1") +
  xlim(c(-2.7,0.5))+
  ylim(c(-2.3,0.5))
dev.off()

png('25Q3_RESULTS/corr_CDK4_CCND1_breast.png', width = 1000, height = 1000, res = 250)
ggscatter(merged_sel, x = "CDK4", y = "CCND1", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "CDK4", ylab = "CCND1") +
  xlim(c(-4.5,0.5))+
  ylim(c(-4.5,0.5))
dev.off()

png('25Q3_RESULTS/corr_CDK9_CCNK_breast.png', width = 1000, height = 1000, res = 250)
ggscatter(merged_sel, x = "CDK9", y = "CCNK", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = FALSE, cor.method = "pearson",
          xlab = "CDK9", ylab = "CCNK") +
  xlim(c(-2.5,0.5))+
  ylim(c(-2.5,0.5))
dev.off()

png('25Q3_RESULTS/corr_CDK6_CCND3_myeloid.png', width = 1000, height = 1000, res = 250)
ggscatter(merged_sel, x = "CDK6", y = "CCND3", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = FALSE, cor.method = "pearson",
          xlab = "CDK6", ylab = "CCND3") +
  xlim(c(-3,0.5))+
  ylim(c(-2.5,0.5))
dev.off()


# ---------------- Plotting specific correlations by lienage---------------- 

res <- read_excel('/Users/asanchezb/Desktop/PROJECTS/CODEPENDENCIES_CDKS_CYCLINS/PAPER/Supplementary_tables/TableS2_CDKs_Cyclins_correlation.xlsx', sheet = 4)

# Volcan plot
res$delabel <- ifelse((res$Corr > 0.2 & res$Padj < 0.05), 'Positive', ifelse((res$Corr < -0.2 & res$Padj < 0.05), 'Negative','NA'))
res$gene_lab <- ifelse((abs(res$Corr) > 0.2 & res$Padj < 0.05), res$Pair, '')
res$Target <-  ifelse(( abs(res$Corr) > 0.2 & res$Padj < 0.05), 'Yes', 'No')

max <- max(res$Corr)
min <- min(res$Corr)
val <- ifelse(abs(min) > abs(max), abs(min),abs(max))


png('/Users/asanchezb/Desktop/test.png',width=2000, height=1000, res = 300)
p <- ggplot(data=res, aes(x=Tumor_type, y=Corr, label=gene_lab)) + 
  geom_point(aes(shape=delabel,fill=delabel, size=delabel,color=Target)) +
  theme_minimal() +
  scale_shape_manual(values=c(20,21,21))+
  #geom_text_repel(box.padding = 0.2,  
  #                segment.size = 0.2,  
   #             segment.color = "black",
    #              max.overlaps = 100,
     #             size = 2) + 
  scale_fill_manual(values = c("Positive" = "#FF6666",  
                               "Negative" = "#1874CD"), 
                    na.value = "#E0E0E0" ) +
  scale_color_manual(values=c("#E0E0E0","black"))+
  scale_size_manual(values=c(0.8,1.8,1.8))+
  geom_vline(xintercept=c(0),col="#C0C0C0", size=0.2) +
  ylim(c(-val,val)) +
  xlab("Lineage") + 
  ylab("-log10(Padj)") +
  theme(panel.grid = element_blank(),
        plot.margin = margin(3, 3, 3, 3, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color = "black"))+
  guides(color = FALSE, shape = FALSE, fill=FALSE, size=FALSE) 
print(p)
dev.off()

png('/Users/asanchezb/Desktop/test.png',width=1400, height=2000, res = 300)
p <- ggplot(data=res, aes(x=Corr, y=Tumor_type, label=gene_lab)) + 
  geom_point(aes(shape=delabel,fill=delabel, size=delabel,color=Target)) +
  theme_minimal() +
  scale_shape_manual(values=c(20,21,21))+
  #geom_text_repel(box.padding = 0.2,  
  #                segment.size = 0.2,  
  #             segment.color = "black",
  #              max.overlaps = 100,
  #             size = 2) + 
  scale_fill_manual(values = c("Positive" = "#FF6666",  
                               "Negative" = "#1874CD"), 
                    na.value = "#E0E0E0" ) +
  scale_color_manual(values=c("#E0E0E0","black"))+
  scale_size_manual(values=c(0.8,1.8,1.8))+
  geom_vline(xintercept=c(0),col="#C0C0C0", size=0.2) +
  xlim(c(-val,val)) +
  ylab("Lineage") + 
  xlab("-log10(Padj)") +
  theme(panel.grid = element_blank(),
        plot.margin = margin(3, 3, 3, 3, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color = "black"))+
  guides(color = FALSE, shape = FALSE, fill=FALSE, size=FALSE) 
print(p)
dev.off()


# ---------------- Determine optimal n of clusters ---------------- 

# Using CDKs and Cyclins (Only CDKs and only Cyclins work worst)

# Determine the optimal number of clusters
png('25Q3_RESULTS/silhouette_CDKs_Cyclins.png',width=1200, height=1000, res = 250)
fviz_nbclust(mat, kmeans, method = "silhouette") + 
  labs(subtitle = "Elbow Method")
dev.off()  

# Determine the optimal number of clusters
png('25Q3_RESULTS/gap_stat_CDKs_Cyclins.png',width=1200, height=1000, res = 250)
fviz_nbclust(mat, kmeans, method = "gap_stat") + 
  labs(subtitle = "Elbow Method")
dev.off()  

# Determine the optimal number of clusters
png('25Q3_RESULTS/wss_CDKs_Cyclins.png',width=1200, height=1000, res = 250)
fviz_nbclust(mat, kmeans, method = "wss") + 
  labs(subtitle = "Elbow Method")
dev.off()  

# 6 cluster selected by gap_stat and wss methods, 2 by silhouette


# ---------------- Clustering  CDKs and Cyclins ---------------- 

# seed the clustering to be reproducible
set.seed(123)
# number of clusters selected
k <- 6
# kmeans clustering
km_res <- kmeans(mat, centers = k, nstart = 25)
# Append cluster assignment to the data
cluster_assignment <- km_res$cluster
cluster_colors <- c("1" = "#CD3333", "2" = "#1874CD", "3" = "#CD9B1D", "4" = "#CD6090", "5" = "#EE9572", "6" = "#7A378B")

# Colors for mutations 
unique_lineages <- unique(met_filt$RB1_status)
color_palette <- c( "#EEE9E9","#424242")  
mut_colors <- setNames(color_palette[seq_along(unique_lineages)], unique_lineages)

# Colors for lineage
unique_lineages <- unique(met_filt$OncotreePrimaryDisease)
color_palette <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique_lineages))
lineage_colors <- setNames(color_palette[seq_along(unique_lineages)], unique_lineages)

# Colors for dependency
col_fun <- colorRamp2(c(-5, 0, 5), c("#0B0B8C", "white", "#A90C38")) #"#0B0B8C", "white", "#A90C38"

# Plot heatmap
png('25Q3_RESULTS/Pancancer_CDKs_Cyclins_2.png',width=900, height=1500, res = 100) #width=2000, height=500, res = 150
draw(
  Heatmap(
    mat, 
    name = "Gene effect",
    show_row_names = TRUE, 
    show_column_names = TRUE, 
    column_names_rot = 55,
    column_names_gp = gpar(fontsize = 8.5),
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    #column_split = cluster_assignment,
    row_title = "Genes",
    column_title = "Cell Lines",
    col = col_fun
  ),
  show_heatmap_legend = FALSE,
  show_annotation_legend = FALSE
)
dev.off()

#############################################################################################################
# 3. CLUSTERING
#############################################################################################################

# Filter metadata 
met_filt<- metadata[metadata$ModelID %in% merged$ModelID,]

# Mantain only dependencies values
merged_filt <- merged %>%
  dplyr::select(-all_of(select))

# Remove SPDYE1 due to NA values
merged_filt <- merged_filt %>%
  select(-SPDYE1)

# ---------------- Determine optimal n of clusters ---------------- 

# Using CDKs and Cyclins (Only CDKs and only Cyclins work worst)

# Determine the optimal number of clusters
png('25Q3_RESULTS/silhouette_CDKs_Cyclins.png',width=1200, height=1000, res = 250)
fviz_nbclust(merged_filt, kmeans, method = "silhouette") + 
  labs(subtitle = "Elbow Method")
dev.off()  

# Determine the optimal number of clusters
png('25Q3_RESULTS/gap_stat_CDKs_Cyclins.png',width=1200, height=1000, res = 250)
fviz_nbclust(merged_filt, kmeans, method = "gap_stat") + 
  labs(subtitle = "Elbow Method")
dev.off()  

# Determine the optimal number of clusters
png('25Q3_RESULTS/wss_CDKs_Cyclins.png',width=1200, height=1000, res = 250)
fviz_nbclust(merged_filt, kmeans, method = "wss") + 
  labs(subtitle = "Elbow Method")
dev.off()  

# 6 cluster selected by gap_stat and wss methods, 2 by silhouette


# ---------------- Clustering  CDKs and Cyclins ---------------- 

# seed the clustering to be reproducible
set.seed(123)
# number of clusters selected
k <- 6
# kmeans clustering
km_res <- kmeans(merged_filt, centers = k, nstart = 25)
# Append cluster assignment to the data
cluster_assignment <- km_res$cluster
cluster_colors <- c("1" = "#CD3333", "2" = "#1874CD", "3" = "#CD9B1D", "4" = "#CD6090", "5" = "#EE9572", "6" = "#7A378B")

# Colors for mutations 
unique_lineages <- unique(met_filt$RB1_status)
color_palette <- c( "#EEE9E9","#424242")  
mut_colors <- setNames(color_palette[seq_along(unique_lineages)], unique_lineages)

# Colors for lineage
unique_lineages <- unique(met_filt$OncotreePrimaryDisease)
color_palette <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique_lineages))
lineage_colors <- setNames(color_palette[seq_along(unique_lineages)], unique_lineages)

# Colors for dependency
col_fun <- colorRamp2(c(-3, 0, 1), c("#A90C38", "white", "#0B0B8C")) #"#0B0B8C", "white", "#A90C38"

group_colors <- c(
  "RBmut" = "#724F13",
  "RBupstream" = "#D59D3B",
  "WT" = "#E8C584"
)

# Plot heatmap
png('25Q3_RESULTS/FigS3A.png',width=2000, height=1400, res = 150) #width=2000, height=500, res = 150
draw(
  Heatmap(
    t(merged_filt), 
    name = "Gene effect",
    show_row_names = TRUE, 
    show_column_names = FALSE, 
    column_names_rot = 55,
    column_names_gp = gpar(fontsize = 8.5),
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    column_split = cluster_assignment,
    row_title = "Genes",
    column_title = "Cell Lines",
    top_annotation = HeatmapAnnotation(
      Cluster = as.factor(cluster_assignment),
      Subtype = met_filt$OncotreePrimaryDisease,
      Group = met_filt$Group,
      TP53 = met_filt$TP53_status,
      RB1 = met_filt$RB1_act,
      WNT_pathway = met_filt$WNT_pathway,
      PI3K_pathway = met_filt$PI3,
      RB1_upstream = met_filt$Rb_pathway,
      RBL1_2_status = met_filt$RBL1_2_status,
      col = list(
        Cluster = cluster_colors,
        Subtype = lineage_colors,
        Group = group_colors,
        WNT_pathway = mut_colors,
        TP53 = mut_colors,
        RB1 = mut_colors,
        PI3K_pathway = mut_colors,
        RB1_upstream = mut_colors,
        RBL1_2_status = mut_colors
      )
    ),
    col = col_fun
  ),
  show_heatmap_legend = FALSE,
  show_annotation_legend = FALSE
)
dev.off()

# ---------------- Most variable genes selection ---------------- 

# Variance by gene
gene_var <- apply(merged_filt, 2, var, na.rm = TRUE)
# Order 
gene_var_sorted <- sort(gene_var, decreasing = TRUE)
# top 15
head(gene_var_sorted, n = 15)
df <- data.frame(
  Gene = names(gene_var_sorted),
  Var  = gene_var_sorted
)
write.table(df, '25Q3_RESULTS/CDKs_Cyclins_variance.csv', sep = ',', row.names = F)

# Selected
# Removig essential genes CDK1, CDK7, CCNA2,CDK9,CCNH,CCNK
# CCNE2 was added due to similarity to CCNE1
main <- c('CDK2','CDK4','CDK6','CCND1','CCND2','CCND3', 'CCNE1', 'CCNE2','CCNL1','CCNC','CDK8','CCNB1','CCNF','CDK12','CDK13')


#selected genes
merged_sel <- merged_filt[,main]


#---------------- Determine optimal n of clusters selected ---------------- 
  
# Using CDKs and Cyclins (Only CDKs and only Cyclins work worst)
  
# Determine the optimal number of clusters
png('25Q3_RESULTS/silhouette_selected_2.png',width=1200, height=1000, res = 250)
fviz_nbclust(merged_sel, kmeans, method = "silhouette") + 
  labs(subtitle = "silhouette Method")
dev.off()  

# Determine the optimal number of clusters
png('25Q3_RESULTS/gap_stat_selected_2.png',width=1200, height=1000, res = 250)
fviz_nbclust(merged_sel, kmeans, method = "gap_stat") + 
  labs(subtitle = "gap_stat Method")
dev.off()  

# Determine the optimal number of clusters
png('25Q3_RESULTS/wss_selected_2.png',width=1200, height=1000, res = 250)
fviz_nbclust(merged_sel, kmeans, method = "wss") + 
  labs(subtitle = "wss Method")
dev.off()  

# 6 cluster selected due to biological sense

# ---------------- Clustering selected ---------------- 

# seed the clustering to be reproducible
set.seed(124)
# number of clusters selected
k <- 6
# kmeans clustering
km_res <- kmeans(merged_sel, centers = k, nstart = 25)
# Append cluster assignment to the data
cluster_assignment <- km_res$cluster
cluster_colors <- c("1" = "#CD3333", "2" = "#1874CD", "3" = "#CD9B1D", "4" = "#CD6090", "5" = "#EE9572", "6" = "#7A378B")
# Add to metadata
met_filt$cluster <- cluster_assignment
table(met_filt$cluster)
#write.table(met_filt[,c('ModelID','cluster')], '25Q3_RESULTS/cluster_assigment.csv', sep = ',', row.names = F)
# Colors for mutations 
unique_lineages <- unique(met_filt$RB1_status)
color_palette <- c( "#EEE9E9","#424242")  
mut_colors <- setNames(color_palette[seq_along(unique_lineages)], unique_lineages)

# Colors for lineage
unique_lineages <- unique(met_filt$OncotreePrimaryDisease)
color_palette <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique_lineages))
lineage_colors <- setNames(color_palette[seq_along(unique_lineages)], unique_lineages)

# Colors for dependency
col_fun <- colorRamp2(c(-3, 0, 1), c("#A90C38", "white", "#0B0B8C")) #"#0B0B8C", "white", "#A90C38"

# Plot heatmap
png('25Q3_RESULTS/Fig3S.png',width=2000, height=800, res = 150) #width=2000, height=500, res = 150
draw(
  Heatmap(
    t(merged_sel), 
    name = "Gene effect",
    show_row_names = TRUE, 
    show_column_names = FALSE, 
    column_names_rot = 55,
    column_names_gp = gpar(fontsize = 8.5),
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    column_split = cluster_assignment,
    row_title = "Genes",
    column_title = "Cell Lines",
    top_annotation = HeatmapAnnotation(
      Cluster = as.factor(cluster_assignment),
      Subtype = met_filt$OncotreePrimaryDisease,
      Group = met_filt$Group,
      TP53 = met_filt$TP53_status,
      RB = met_filt$RB1_act,
      WNT_pathway = met_filt$WNT_pathway,
      PI3K_pathway = met_filt$PI3,
      RB_pathway = met_filt$Rb_pathway,
      RBL1_2_status = met_filt$RBL1_2_status,
      col = list(
        Cluster = cluster_colors,
        Subtype = lineage_colors,
        Group = group_colors,
        WNT_pathway = mut_colors,
        TP53 = mut_colors,
        RB = mut_colors,
        PI3K_pathway = mut_colors,
        RB_pathway = mut_colors,
        RBL1_2_status = mut_colors
      )
    ),
    col = col_fun
  ),
  #show_heatmap_legend = FALSE,
  show_annotation_legend = FALSE
)
dev.off()


# ---------------- Re-Clustering selected ---------------- 

rownames(merged_sel) <- met_filt$ModelID
merged_new <- merged_sel[rownames(merged_sel) %in% met_filt[met_filt$cluster == 6,]$ModelID,]
#merged_new <- merged_sel[rownames(merged_sel) %in% met_filt[met_filt$cluster %in% c(4,5),]$ModelID,]


# Determine the optimal number of clusters
png('25Q3_RESULTS/silhouette_selected_3.png',width=1200, height=1000, res = 250)
fviz_nbclust(merged_new, kmeans, method = "silhouette") + 
  labs(subtitle = "silhouette Method")
dev.off()  

# Determine the optimal number of clusters
png('25Q3_RESULTS/gap_stat_selected_3.png',width=1200, height=1000, res = 250)
fviz_nbclust(merged_new, kmeans, method = "gap_stat") + 
  labs(subtitle = "gap_stat Method")
dev.off()  

# Determine the optimal number of clusters
png('25Q3_RESULTS/wss_selected_3.png',width=1200, height=1000, res = 250)
fviz_nbclust(merged_new, kmeans, method = "wss") + 
  labs(subtitle = "wss Method")
dev.off()  

# seed the clustering to be reproducible
set.seed(124)
# number of clusters selected
k <- 6
# kmeans clustering
km_res <- kmeans(merged_new, centers = k, nstart = 25)
# Append cluster assignment to the data
cluster_assignment <- km_res$cluster
cluster_colors <- c("1" = "#CD3333", "2" = "#1874CD", "3" = "#CD9B1D", "4" = "#CD6090", "5" = "#EE9572", "6" = "#7A378B")
# Add to metadata
met_new <- met_filt[met_filt$cluster == 6,]
#met_new <- met_filt[met_filt$cluster %in% c(4,5), ]
met_new$cluster <- cluster_assignment
table(met_new$cluster)
#write.table(met_filt[,c('ModelID','cluster')], '25Q3_RESULTS/cluster_assigment.csv', sep = ',', row.names = F)
# Colors for mutations 
unique_lineages <- unique(met_new$RB1_status)
color_palette <- c( "#EEE9E9","#424242")  
mut_colors <- setNames(color_palette[seq_along(unique_lineages)], unique_lineages)

# Colors for lineage
unique_lineages <- unique(met_new$OncotreePrimaryDisease)
color_palette <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique_lineages))
lineage_colors <- setNames(color_palette[seq_along(unique_lineages)], unique_lineages)

# Colors for dependency
col_fun <- colorRamp2(c(-3, 0, 1), c("#A90C38", "white", "#0B0B8C")) #"#0B0B8C", "white", "#A90C38"

# Plot heatmap
png('25Q3_RESULTS/Pancancer_CDKs_Cyclins_sel_other_4.png',width=1200, height=800, res = 150) #width=2000, height=500, res = 150
draw(
  Heatmap(
    t(merged_new), 
    name = "Gene effect",
    show_row_names = TRUE, 
    show_column_names = FALSE, 
    column_names_rot = 55,
    column_names_gp = gpar(fontsize = 8.5),
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    column_split = cluster_assignment,
    row_title = "Genes",
    column_title = "Cell Lines",
    top_annotation = HeatmapAnnotation(
      Cluster = as.factor(cluster_assignment),
      Subtype = met_new$OncotreePrimaryDisease,
      TP53 = met_new$TP53_status,
      RB = met_new$RB1_status,
      WNT_pathway = met_new$WNT_pathway,
      PI3K_pathway = met_new$PI3K_pathway,
      CDK4_6_pathway = met_new$CDK4_6_pathway,
      RBL1_2_status = met_new$RBL1_2_status,
      col = list(
        Cluster = cluster_colors,
        Subtype = lineage_colors,
        WNT_pathway = mut_colors,
        TP53 = mut_colors,
        RB = mut_colors,
        PI3K_pathway = mut_colors,
        CDK4_6_pathway = mut_colors,
        RBL1_2_status = mut_colors
      )
    ),
    col = col_fun
  ),
  #show_heatmap_legend = FALSE,
  show_annotation_legend = FALSE
)
dev.off()

# rename

met_new$cluster[met_new$cluster == 2] <- 7
met_new$cluster[met_new$cluster == 3] <- 8
met_new$cluster[met_new$cluster == 4] <- 9
met_new$cluster[met_new$cluster == 5] <- 10
met_new$cluster[met_new$cluster == 6] <- 11
met_new$cluster[met_new$cluster == 1] <- 6

idx <- met_filt$cluster == 6

met_filt$cluster[idx] <- 
  met_new$cluster[ match(met_filt$ModelID[idx],
                         met_new$ModelID) ]

table(met_filt$cluster)
cluster_colors <- c("1" = "#CD3333", "2" = "#1874CD", "3" = "#CD9B1D", "4" = "#7A378B", "5" = "#EE9572", "6" = "#17301D","7" = "#264F30", "8" = "#3E804F","9" = "#72BB84", "10" = "#B0D9BA", "11" = "#CFE8D5")
# Plot heatmap
png('25Q3_RESULTS/Fig3A.png',width=2000, height=900, res = 150) #width=2000, height=500, res = 150
draw(
  Heatmap(
    t(merged_sel), 
    name = "Gene effect",
    show_row_names = TRUE, 
    show_column_names = FALSE, 
    column_names_rot = 55,
    column_names_gp = gpar(fontsize = 8.5),
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    column_split = met_filt$cluster,
    row_title = "Genes",
    column_title = "Cell Lines",
    top_annotation = HeatmapAnnotation(
      Cluster = met_filt$cluster,
      #Subtype = met_filt$OncotreePrimaryDisease,
      Group = met_filt$Group,
      TP53 = met_filt$TP53_status,
      RB = met_filt$RB1_act,
      WNT_pathway = met_filt$WNT_pathway,
      PI3K_pathway = met_filt$PI3,
      RB_pathway = met_filt$Rb_pathway,
      RBL1_2_status = met_filt$RBL1_2_status,
      col = list(
        Cluster = cluster_colors,
        Subtype = lineage_colors,
        Group = group_colors,
        WNT_pathway = mut_colors,
        TP53 = mut_colors,
        RB = mut_colors,
        PI3K_pathway = mut_colors,
        RB_pathway = mut_colors,
        RBL1_2_status = mut_colors
      )
    ),
    col = col_fun
  ),
  #show_heatmap_legend = FALSE,
  show_annotation_legend = FALSE
)
dev.off()

# rename 

met_filt$cluster[met_filt$cluster == 1] <- 'A'
met_filt$cluster[met_filt$cluster == 2] <- 'B'
met_filt$cluster[met_filt$cluster == 3] <- 'C'
met_filt$cluster[met_filt$cluster == 4] <- 'D'
met_filt$cluster[met_filt$cluster == 5] <- 'E'
met_filt$cluster[met_filt$cluster == 6] <- 'F1'
met_filt$cluster[met_filt$cluster == 7] <- 'F2'
met_filt$cluster[met_filt$cluster == 8] <- 'F3'
met_filt$cluster[met_filt$cluster == 9] <- 'F4'
met_filt$cluster[met_filt$cluster == 10] <- 'F5'
met_filt$cluster[met_filt$cluster == 11] <- 'F6'

table(met_filt$cluster)

# save metadata filtering
write.table(met_filt, '/Users/asanchezb/Desktop/met_filt_030326.tsv', sep = '\t')

met_filt <- read.table('/Users/asanchezb/Desktop/met_filt_030326.tsv',sep = '\t',header = TRUE)

raw_drug_screen$ModelID <- rownames(raw_drug_screen)
met_filt_drug_new <- merge(met_filt_drug, 
                       raw_drug_screen,
                       by = "ModelID",
                       all.x = TRUE)
# Add drugs
col_drug <- colorRamp2(c(0.75, 1), c("#8B7500","white")) #"#0B0B8C", "white", "#A90C38"

# Plot heatmap
png('25Q3_RESULTS/Fig_drugs.png',width=2000, height=900, res = 150) #width=2000, height=500, res = 150
draw(
  Heatmap(
    t(merged_sel), 
    name = "Gene effect",
    show_row_names = TRUE, 
    show_column_names = FALSE, 
    column_names_rot = 55,
    column_names_gp = gpar(fontsize = 8.5),
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    column_split = met_filt$cluster,
    row_title = "Genes",
    column_title = "Cell Lines",
    top_annotation = HeatmapAnnotation(
      Cluster = met_filt_drug$cluster,
      #Subtype = met_filt$OncotreePrimaryDisease,
      Group = met_filt_drug$Group,
      TP53 = met_filt_drug$TP53_status,
      RB = met_filt_drug$RB1_act,
      WNT_pathway = met_filt_drug$WNT_pathway,
      PI3K_pathway = met_filt_drug$PI3,
      RB_pathway = met_filt_drug$Rb_pathway,
      RBL1_2_status = met_filt_drug$RBL1_2_status,
      Palbociclib = met_filt_drug$palbociclib,
      Abemaciclib = met_filt_drug$abemaciclib,
      Ribociclib = met_filt_drug$ribociclib,
      col = list(
        Cluster = cluster_colors,
        Subtype = lineage_colors,
        Group = group_colors,
        WNT_pathway = mut_colors,
        TP53 = mut_colors,
        RB = mut_colors,
        PI3K_pathway = mut_colors,
        RB_pathway = mut_colors,
        RBL1_2_status = mut_colors,
        Palbociclib = col_drug,
        Abemaciclib = col_drug,
        Ribociclib = col_drug
      )
    ),
    col = col_fun
  ),
  #show_heatmap_legend = FALSE,
  show_annotation_legend = FALSE
)
dev.off()


# filtering NA
met_filt_drug_only <- met_filt_drug_new[
  !is.na(met_filt_drug_new$C0.15_uM) ,
  ]

met_filt_drug_only <- met_filt_drug[
  !is.na(met_filt_drug$palbociclib) &
    !is.na(met_filt_drug$abemaciclib) &
    !is.na(met_filt_drug$ribociclib),
  ]

merged_sel_drug <- merged_sel[rownames(merged_sel)%in% met_filt_drug_only$ModelID,]

# Crear data frame auxiliar para ordenar
ordering_df <- data.frame(
  cluster = met_filt_drug_only$cluster,
  palbo = met_filt_drug_only$C0.6_uM
)

# Orden por cluster y luego por palbociclib
ord <- order(ordering_df$cluster, ordering_df$C0.15_uM)

# Reordenar todo de forma consistente
merged_sel_drug <- merged_sel_drug[ord,]
met_filt_drug_only <- met_filt_drug_only[ord, ]

max(met_filt_drug_only$C0.15_uM)
min(met_filt_drug_only$C0.15_uM)

max(met_filt_drug_only$C0.6_uM)
min(met_filt_drug_only$C0.6_uM)

col_drug <- colorRamp2(c(-1, 1), c("#8B7500","white")) #"#0B0B8C", "white", "#A90C38"
# Plot heatmap
png('25Q3_RESULTS/Fig_drugs_2.png',width=2000, height=900, res = 150) #width=2000, height=500, res = 150
draw(
  Heatmap(
    t(merged_sel_drug), 
    name = "Gene effect",
    show_row_names = TRUE, 
    show_column_names = FALSE, 
    column_names_rot = 55,
    column_names_gp = gpar(fontsize = 8.5),
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    column_split = met_filt_drug_only$cluster,
    row_title = "Genes",
    column_title = "Cell Lines",
    top_annotation = HeatmapAnnotation(
      Cluster = met_filt_drug_only$cluster,
      #Subtype = met_filt$OncotreePrimaryDisease,
      #Group = met_filt_drug_only$Group,
      TP53 = met_filt_drug_only$TP53_status,
      RB1 = met_filt_drug_only$RB1_act,
      WNT_pathway = met_filt_drug_only$WNT_pathway,
      PI3K_pathway = met_filt_drug_only$PI3,
      RB_upstream = met_filt_drug_only$Rb_pathway,
      RBL1_2_status = met_filt_drug_only$RBL1_2_status,
      #Palbociclib = met_filt_drug_only$palbociclib,
      #Abemaciclib = met_filt_drug_only$abemaciclib,
      Palbo_0.15 = met_filt_drug_only$C0.15_uM,
      Palbo_0.60 = met_filt_drug_only$C0.6_uM,
      #Ribociclib = met_filt_drug_only$ribociclib,
      col = list(
        Cluster = cluster_colors,
        #Subtype = lineage_colors,
        #Group = group_colors,
        WNT_pathway = mut_colors,
        TP53 = mut_colors,
        RB1 = mut_colors,
        PI3K_pathway = mut_colors,
        RB_upstream = mut_colors,
        RBL1_2_status = mut_colors,
        Palbociclib = col_drug,
        Abemaciclib = col_drug,
        Ribociclib = col_drug,
        Palbo_0.15 = col_drug,
        Palbo_0.60 = col_drug
      )
    ),
    col = col_fun
  ),
  #show_heatmap_legend = FALSE,
  show_annotation_legend = FALSE
)
dev.off()

ord <- order(met_filt_drug_only$C0.6_uM)

merged_sel_drug <- merged_sel_drug[ord,]
met_filt_drug_only <- met_filt_drug_only[ord, ]
col_drug <- colorRamp2(c(0.75, 1), c("#8B7500","white")) #"#0B0B8C", "white", "#A90C38"
# Plot heatmap
png('25Q3_RESULTS/Fig_drugs_all_TEST.png',width=2000, height=900, res = 150) #width=2000, height=500, res = 150
draw(
  Heatmap(
    t(merged_sel_drug), 
    name = "Gene effect",
    show_row_names = TRUE, 
    show_column_names = FALSE, 
    column_names_rot = 55,
    column_names_gp = gpar(fontsize = 8.5),
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    #column_split = met_filt_drug_only$cluster,
    row_title = "Genes",
    column_title = "Cell Lines",
    top_annotation = HeatmapAnnotation(
      Cluster = met_filt_drug_only$cluster,
      #Subtype = met_filt$OncotreePrimaryDisease,
      #Group = met_filt_drug_only$Group,
      TP53 = met_filt_drug_only$TP53_status,
      RB1 = met_filt_drug_only$RB1_act,
      WNT_pathway = met_filt_drug_only$WNT_pathway,
      PI3K_pathway = met_filt_drug_only$PI3,
      RB_upstream = met_filt_drug_only$Rb_pathway,
      RBL1_2_status = met_filt_drug_only$RBL1_2_status,
      Palbo_0.15 = met_filt_drug_only$C0.15_uM,
      Palbo_0.60 = met_filt_drug_only$C0.6_uM,
      #Palbociclib = met_filt_drug_only$palbociclib,
      #Abemaciclib = met_filt_drug_only$abemaciclib,
      #Ribociclib = met_filt_drug_only$ribociclib,
      col = list(
        Cluster = cluster_colors,
        #Subtype = lineage_colors,
        #Group = group_colors,
        WNT_pathway = mut_colors,
        TP53 = mut_colors,
        RB1 = mut_colors,
        PI3K_pathway = mut_colors,
        RB_upstream = mut_colors,
        RBL1_2_status = mut_colors,
        Palbociclib = col_drug,
        Abemaciclib = col_drug,
        Ribociclib = col_drug,
        Palbo_0.15 = col_drug,
        Palbo_0.60 = col_drug
      )
    ),
    col = col_fun
  ),
  #show_heatmap_legend = FALSE,
  show_annotation_legend = TRUE
)
dev.off()

#############################################################################################################
# 4. MUTATION ANALYSIS
#############################################################################################################

met_filt$cluster_2 <- ifelse(met_filt$cluster == 'D', 'D', 'rest')

# ---------------- TMB ---------------- 

# Relevant mutation by cell line
mutation_by_model <- mutaciones_relevants %>%
  group_by(ModelID) %>%
  summarise(num = n())

# Longitud del exoma en megabases (Mb)
lg_genome_Mb <- 50  # Ajusta según corresponda

# Calcular el TMB por cada línea celular
mutation_by_model <- mutation_by_model %>%
  mutate(TMB = num / lg_genome_Mb)

met_filt_2 <- merge(met_filt, mutation_by_model, by.x = 'ModelID',by.y = 'ModelID')

# ---------------- Plotting Pancancer ---------------- 

# Variables a graficar
vars <- c("TMB","MSIScore","Ploidy","LoHFraction","CIN","Aneuploidy")

# Función para calcular densidad en la mediana
calculate_density_at_median <- function(data, var, median_value) {
  dens <- density(data[[var]], from = min(data[[var]], na.rm=TRUE),
                  to   = max(data[[var]], na.rm=TRUE),
                  na.rm = TRUE)
  approx(dens$x, dens$y, xout = median_value, rule = 2)$y
}

# Loop por cada variable
for (var in vars) {
  
  message("Procesando: ", var)
  
  # Filtrar datos sin NA para la variable
  df_clean <- met_filt_2 %>% filter(!is.na(.data[[var]]))
  
  # Calcular medianas y densidades
  medianas <- df_clean %>%
    group_by(cluster_2) %>%
    summarise(median_value = median(.data[[var]], na.rm = TRUE)) %>%
    rowwise() %>%
    mutate(density_at_median = calculate_density_at_median(
      df_clean %>% filter(cluster_2 == cluster_2),
      var,
      median_value
    ))
  
  # Nombre del archivo
  file_name <- paste0("25Q3_RESULTS/", var, "_pancancer.png")
  
  png(file_name, width=1000, height=1000, res=220)
  
  print(
    ggplot() +
      # Densidad cluster 6
      geom_density(
        data = subset(met_filt_2, cluster_2 == 'D'),
        aes(x = .data[[var]], y = ..density..),
        fill = "#7A378B", alpha = 0.6
      ) +
      geom_segment(
        data = medianas %>% filter(cluster_2 == 'D'),
        aes(x = median_value, xend = median_value,
            y = 0, yend = density_at_median - 0.01),
        linetype = "dashed", color = "black", size = 0.5
      ) +
      # Etiqueta de la mediana cluster 6
      geom_text(
        data = medianas %>% filter(cluster_2 == 'D'),
        aes(x = median_value, 
            y = density_at_median + 0.02,
            label = paste0(round(median_value, 2))),
        size = 4
      ) +
      
      # Densidad grupo rest invertida
      geom_density(
        data = subset(met_filt_2, cluster_2 == 'rest'),
        aes(x = .data[[var]], y = -..density..),
        fill = "#7F7F7F", alpha = 0.6
      ) +
      geom_segment(
        data = medianas %>% filter(cluster_2 == 'rest'),
        aes(x = median_value, xend = median_value,
            y = 0, yend = -density_at_median),
        linetype = "dashed", color = "black", size = 0.5
      ) +
      # Etiqueta de la mediana cluster rest
      geom_text(
        data = medianas %>% filter(cluster_2 == 'rest'),
        aes(x = median_value,
            y = -density_at_median - 0.02,
            label = paste0(round(median_value, 2))),
        size = 4
      ) +
      
      xlab(var) +
      ylab("Density") +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        plot.margin = margin(10,10,10,10,"pt")
      )
  )
  
  dev.off()
}



# ----------------  Plotting by Lineage ----------------

calculate_density_at_median <- function(data, var, median_value) {
  dens <- density(data[[var]],
                  from = min(data[[var]], na.rm = TRUE),
                  to   = max(data[[var]], na.rm = TRUE),
                  na.rm = TRUE)
  approx(dens$x, dens$y, xout = median_value, rule = 2)$y
}
plot_density_by_category <- function(df, category_name, category_value, out_dir) {
  
  vars <- c("TMB","MSIScore","Ploidy","LoHFraction","CIN","Aneuploidy")
  
  for (var in vars) {
    
    message("Procesando: ", category_name, " = ", category_value, " | ", var)
    
    df_clean <- df %>% filter(!is.na(.data[[var]]))
    
    # Necesitamos ambos clusters
    if (!all(c("6", "rest") %in% df_clean$cluster_2)) next
    
    # Calcular medianas y densidades
    medianas <- df_clean %>%
      group_by(cluster_2) %>%
      summarise(median_value = median(.data[[var]], na.rm = TRUE)) %>%
      rowwise() %>%
      mutate(density_at_median = calculate_density_at_median(
        df_clean %>% filter(cluster_2 == cluster_2),
        var,
        median_value
      ))
    
    file_name <- file.path(
      out_dir,
      paste0(var, "_", gsub("[^A-Za-z0-9]", "_", category_value), ".png")
    )
    
    png(file_name, width = 1000, height = 1000, res = 220)
    
    print(
      ggplot() +
        geom_density(
          data = subset(df_clean, cluster_2 == "6"),
          aes(x = .data[[var]], y = ..density..),
          fill = "#7A378B", alpha = 0.6
        ) +
        geom_segment(
          data = medianas %>% filter(cluster_2 == "6"),
          aes(x = median_value, xend = median_value,
              y = 0, yend = density_at_median),
          linetype = "dashed", color = "black"
        ) +
        geom_text(
          data = medianas %>% filter(cluster_2 == "6"),
          aes(x = median_value,
              y = density_at_median + 0.02,
              label = round(median_value, 2)),
          size = 4
        ) +
        geom_density(
          data = subset(df_clean, cluster_2 == "rest"),
          aes(x = .data[[var]], y = -..density..),
          fill = "#7F7F7F", alpha = 0.6
        ) +
        geom_segment(
          data = medianas %>% filter(cluster_2 == "rest"),
          aes(x = median_value, xend = median_value,
              y = 0, yend = -density_at_median),
          linetype = "dashed", color = "black"
        ) +
        geom_text(
          data = medianas %>% filter(cluster_2 == "rest"),
          aes(x = median_value,
              y = -density_at_median - 0.02,
              label = round(median_value, 2)),
          size = 4
        ) +
        xlab(var) +
        ylab("Density") +
        theme_minimal() +
        theme(
          panel.grid = element_blank(),
          axis.line = element_line(color = "black")
        )
    )
    
    dev.off()
  }
}
safe_name <- function(x) {
  x %>%
    gsub("/", "_", .) %>%
    gsub(" ", "_", .) %>%
    gsub("[^A-Za-z0-9_]", "", .)
}

dire <- "CODEPENDENCIES_CDKS_CYCLINS/25Q3_RESULTS/Genomics/Lineage"
dir.create(dire, recursive = TRUE, showWarnings = FALSE)

for (lin in unique(met_filt$OncotreeLineage)) {
  
  df_lin <- met_filt %>% filter(OncotreeLineage == lin)
  if (nrow(df_lin) < 7) next
  
  lin_safe <- safe_name(lin)
  out_dir <- file.path(dire, lin_safe)
  dir.create(out_dir, showWarnings = FALSE)
  
  plot_density_by_category(
    df = df_lin,
    category_name = "Lineage",
    category_value = lin_safe,
    out_dir = out_dir
  )
}



dire <- "CODEPENDENCIES_CDKS_CYCLINS/25Q3_RESULTS/Genomics/Primary"
dir.create(dire, recursive = TRUE, showWarnings = FALSE)

for (lin in unique(met_filt$OncotreePrimaryDisease)) {
  
  df_lin <- met_filt %>% filter(OncotreePrimaryDisease == lin)
  if (nrow(df_lin) < 7) next
  
  lin_safe <- safe_name(lin)
  out_dir <- file.path(dire, lin_safe)
  dir.create(out_dir, showWarnings = FALSE)
  
  plot_density_by_category(
    df = df_lin,
    category_name = "Primary",
    category_value = lin_safe,
    out_dir = out_dir
  )
}


#############################################################################################################
# 4. CLUSTERING CHARACTERIZATION
#############################################################################################################


# ---------------- Cluster distribution ---------------- 

table(met_filt$OncotreeLineage,met_filt$cluster)
table(met_filt$OncotreePrimaryDisease,met_filt$cluster)

# create data frame
metadata_filt_pct <- met_filt %>%
  group_by(OncotreeLineage, cluster) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(OncotreeLineage) %>%
  mutate(percentage = count / sum(count)) %>%
  ungroup()
# remove < 7 cell lines
metadata_filt_pct <- metadata_filt_pct %>%
  group_by(OncotreeLineage) %>%
  filter(sum(count) >= 7) %>%
  ungroup()

cluster_colors <- c("A" = "#CD3333", "B" = "#1874CD", "C" = "#CD9B1D", "D" = "#7A378B", "E" = "#EE9572", 
                    "F1" = "#17301D","F2" = "#264F30", "F3" = "#3E804F","F4" = "#72BB84", "F5" = "#B0D9BA", "F6" = "#CFE8D5")

# create barplot
png('25Q3_RESULTS/primary_cluster_distribution.png',width=2000, height=1200, res = 220)
ggplot(metadata_filt_pct, aes(x = OncotreeLineage, y = percentage, fill = as.factor(cluster))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = cluster_colors) + 
  theme_minimal() +
  labs(x = "Lineage", y = "Cell lines", fill = "Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid = element_blank(),
        plot.margin = margin(10, 10, 10, 10, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))
dev.off()

# create barplot
png('25Q3_RESULTS/primary_cluster_distribution.png',width=2000, height=1200, res = 220)
ggplot(metadata_filt_pct, aes(x = OncotreeLineage, y = percentage, fill = as.factor(cluster))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = cluster_colors) + 
  theme_minimal() +
  labs(x = "Lineage", y = "Cell lines", fill = "Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid = element_blank(),
        plot.margin = margin(10, 10, 10, 10, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))
dev.off()


# Contamos cuántas líneas hay por cluster y linaje
metadata_filt_pct_2 <- met_filt %>%
  group_by(cluster, OncotreeLineage) %>%
  summarise(count = n(), .groups = 'drop') %>%
  ungroup()

# Definir linajes mayoritarios (>=7 líneas en total)
major_lineages <- metadata_filt_pct_2 %>%
  group_by(OncotreeLineage) %>%
  summarise(total_count = sum(count), .groups = 'drop') %>%
  filter(total_count >= 7) %>%
  pull(OncotreeLineage)

# Reemplazar los linajes pequeños por "Others"
metadata_filt_pct_2 <- metadata_filt_pct_2 %>%
  mutate(OncotreeLineage = ifelse(OncotreeLineage %in% major_lineages, 
                                  OncotreeLineage, "Others"))

# Recalcular porcentaje dentro de cada cluster
metadata_filt_pct_2 <- metadata_filt_pct_2 %>%
  group_by(cluster, OncotreeLineage) %>%
  summarise(count = sum(count), .groups = 'drop') %>%  # sumar los que son "Others"
  group_by(cluster) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

write.table(metadata_filt_pct_2, '25Q3_RESULTS/primary_tumor_cluster_percentaje.tsv', sep = '\t')

colores_paired_25 <- c(
  "#1f78b4", "#a6cee3",  # Azul / Azul pastel
  "#33a02c", "#b2df8a",  # Verde / Verde pastel
  "#e31a1c", "#fb9a99",  # Rojo / Rojo pastel
  "#ff7f00", "#fdbf6f",  # Naranja / Naranja pastel
  "#6a3d9a", "#cab2d6",  # Púrpura / Lila pastel
  "#b15928", "#fdd0a2",  # Marrón / Marrón pastel
  "#458B74", "#66CDAA",  # Azul pastel / Azul oscuro
  "#8B8878", "#CDC8B1",  # Verde pastel / Verde oscuro
  "#8B7500", "#CDAD00",  # Rojo pastel / Rojo oscuro
  "#8B3A62", "#CD6090",  # Naranja pastel / Naranja oscuro
  "#8B864E", "#CDC673",  # Lila pastel / Púrpura oscuro
  "#8B4C39", "#EE8262" ,"#B8860B"     # Marrón pastel / Marrón oscuro
)

png('25Q3_RESULTS/cluster_lineage_distribution.png',
    width=2000, height=1200, res=220)

ggplot(metadata_filt_pct_2, aes(x = as.factor(cluster), y = percentage, fill = OncotreeLineage)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = colores_paired_25) +  # define tus colores por linaje
  theme_minimal() +
  labs(x = "Cluster", y = "Percentage of cell lines", fill = "Oncotree Lineage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        plot.margin = margin(10, 10, 10, 10, "pt"),
        axis.line = element_line(color = "black"))

dev.off()




# ---------------- Differential expression ---------------- 

breast_cell_lines <- metadata[metadata$OncotreeLineage == 'Breast',]$ModelID
lung_cell_lines <- metadata[metadata$OncotreeLineage == 'Lung',]$ModelID

# cluster D
D_cell_lines <- met_filt[met_filt$cluster == 'D',]$ModelID

# match data
cts <- t(expression_data)

# if tumor type
cts <- cts[,colnames(cts) %in% breast_cell_lines]
cts <- cts[,colnames(cts) %in% lung_cell_lines]

cts <- cts[,colnames(cts) %in% D_cell_lines]

coldata <- met_filt[met_filt$ModelID %in% colnames(cts),]
table(coldata$cluster)
table(met_filt$cluster)

names <- coldata$ModelID
rownames(coldata) <- names
cts <- cts[,colnames(cts) %in% rownames(coldata)]
all(rownames(coldata) %in% colnames(cts))

# order data
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

coldata$cluster3 <- ifelse(coldata$cluster=='F1'|coldata$cluster=='F2'|coldata$cluster=='F3'|coldata$cluster=='F4'|coldata$cluster=='F5'|coldata$cluster=='F6','F','rest')

# Differential expression
coldata$cluster3 <- ifelse(coldata$cluster=='C','C','rest')
coldata$cluster3 <- ifelse(coldata$RB1_status=='mut','mut','WT')


class(coldata$cluster3)
coldata$cluster3 <- as.factor(coldata$cluster3)
coldata$cluster3 <- relevel(coldata$cluster3, ref = "rest")
table(coldata$cluster3)

# expression matrix log(TPM+1)
mat <- as.matrix(cts)

# metadata, convert labels
coldata$cluster3 <- factor(coldata$cluster3,
                           levels = c("F", "rest"),
                           labels = c("F", "rest"))

# design model
design <- model.matrix(~ 0 + coldata$cluster3)
colnames(design) <- levels(coldata$cluster3)

fit <- lmFit(mat, design)

# contrast
contrast.matrix <- makeContrasts(
  F_vs_rest = F - rest,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# results
results <- topTable(fit2, coef = "F_vs_rest", number = Inf)
results$score <- results$logFC * -log10(results$adj.P.Val)
results$gene <- rownames(results)
write.table(results, '25Q3_RESULTS/F_vs_rest_breast.tsv', sep = '\t')

# Plotting
res_ord <- results[order(results$adj.P.Val),]
head(res_ord)
summary(res_ord)

sum(res_ord$adj.P.Val < 0.05, na.rm = TRUE) #### 6167 genes

res_acq <- results

list_genes_int <- c('CDKN2A','CDKN2B','CCNE1','CCNE2','RB1','CCND1','CDK6','E2F3','INSM1','BEX2','CCND2','HES6','SKP2','CCNA1','MYCL','STAT6','MYEOV','RAC2','STEAP1')

list_genes_int <- c(
  "CCND1", "CAV1", "CAV2", "FOSL1", "YAP1",
  "MYOF", "STEAP1", "MYEOV", "RB1","E2F3",
  "KIF5C", "INSM1", "SEZ6", "CHGA", "HES6",
  "BEX1", "CDK5R1", "DLL3", "SEZ6L",
  "DPYSL5", "CCNE2"
)

list_genes_int <- c(
  "EEF1A2", "RB1", "AGR2",
  "CDKN2A", "PLIN2", "KANK4", "CMPK2", "DZIP1",
  "CHI3L1", "GLIPR2", "GFPT2", "NRN1",
  "CHI3L2", "CCNE1"
)

FC <- 0.7 # PANCANCER 0.7, LUNG 1
res_acq$delabel <- ifelse((res_acq$logFC > FC & res_acq$adj.P.Val < 0.05), 'Up', ifelse((res_acq$logFC < -FC & res_acq$adj.P.Val < 0.05), 'Down','NA'))
res_acq$gene_lab <- ifelse(( abs(res_acq$logFC) > FC & res_acq$adj.P.Val < 0.05), rownames(res_acq), '')
res_acq$gene_lab <- ifelse(rownames(res_acq) %in% list_genes_int, rownames(res_acq), '')

res_acq$Target <- ifelse(( abs(res_acq$logFC) > FC & res_acq$adj.P.Val < 0.05),'Yes', 'No')

res_df <- as.data.frame(res_acq)
res_df$genes <- rownames(res_df)

png('25Q3_RESULTS/volcan_plot_cluster6_vs_rest_lung_3.png', width = 1200, height = 1200, res = 250)
ggplot(data=res_df, aes(x=logFC, y=-log10(adj.P.Val), label=gene_lab)) + 
  geom_point(aes(shape=delabel,fill=delabel, size=delabel,color=Target)) +
  theme_minimal() +
  scale_shape_manual(values=c(21,20,21))+
  geom_text_repel(box.padding = 0.7,  # Adjust this value to control the nudge towards the left
                  segment.size = 0.2,  # Adjust the segment size to control the line segment's length
                  segment.color = "black",
                  max.overlaps = Inf,
                  size = 2# Adjust segment color to match the theme
  ) + 

  #scale_fill_manual(values=c("#0080FF", '#E0E0E0',"#FF6666")) +
  scale_fill_manual(values=c("#7F7F7F", '#E0E0E0',"#7A378B")) +
  scale_color_manual(values=c("#E0E0E0","black"))+
  scale_size_manual(values=c(1.7,0.9,1.7))+
  geom_vline(xintercept=c(0),col="#C0C0C0", size=0.2) +
  #geom_hline(yintercept=-log10(0.05), col="black") +
  #ggtitle('Responders VS Non-responders') +
  xlim(c(-3.6,3.6)) +
  xlab("FC") + 
  ylab("-log10(padj)") +
  theme(panel.grid = element_blank(),
        plot.margin = margin(30, 30, 30, 30, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))+
  guides(color = FALSE, shape = FALSE, fill=FALSE, size=FALSE) 
dev.off()

write.table(met_filt[,c('ModelID','cluster')],'/Users/asanchezb/Desktop/clusters.tsv',sep='\t', row.names = F)
library(dplyr)

rb1_lineage_pct <- coldata %>%
  group_by(RB1_status, OncotreePrimaryDisease) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(RB1_status) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

colores_paired_25 <- c(
  "#1f78b4", "#a6cee3",  # Azul / Azul pastel
  "#33a02c", "#b2df8a",  # Verde / Verde pastel
  "#e31a1c", "#fb9a99",  # Rojo / Rojo pastel
  "#ff7f00", "#fdbf6f",  # Naranja / Naranja pastel
  "#6a3d9a", "#cab2d6",  # Púrpura / Lila pastel
  "#b15928", "#fdd0a2",  # Marrón / Marrón pastel
  "#458B74", "#66CDAA",  # Azul pastel / Azul oscuro
  "#8B8878", "#CDC8B1",  # Verde pastel / Verde oscuro
  "#8B7500", "#CDAD00",  # Rojo pastel / Rojo oscuro
  "#8B3A62", "#CD6090",  # Naranja pastel / Naranja oscuro
  "#8B864E", "#CDC673",  # Lila pastel / Púrpura oscuro
  "#8B4C39", "#EE8262" ,"#B8860B",
  "#1f78b4", "#a6cee3",  # Azul / Azul pastel
  "#33a02c", "#b2df8a",  # Verde / Verde pastel
  "#e31a1c", "#fb9a99",  # Rojo / Rojo pastel
  "#ff7f00", "#fdbf6f",  # Naranja / Naranja pastel
  "#6a3d9a", "#cab2d6",  # Púrpura / Lila pastel
  "#b15928", "#fdd0a2",  # Marrón / Marrón pastel
  "#458B74"  # Azul pastel / Azul oscuro# Marrón pastel / Marrón oscuro
)


png('25Q3_RESULTS/twt.png', width = 1300, height = 1200, res = 170)

ggplot(rb1_lineage_pct, aes(x = as.factor(RB1_status), y = percentage, fill = OncotreePrimaryDisease)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = colores_paired_25) +  # define tus colores por linaje
  theme_minimal() +
  labs(x = "Cluster", y = "Percentage of cell lines", fill = "Oncotree Lineage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        plot.margin = margin(10, 10, 10, 10, "pt"),
        axis.line = element_line(color = "black"))

dev.off()


# ---------------- CNV ---------------- 

cnv_relevants <- merge(cnv_data,met_filt,by.x='ModelID',by.y='ModelID')

gene_cols <- colnames(cnv_data)[-1]   # quitando ModelID si era la primera

cnv_long <- cnv_relevants %>%
  pivot_longer(
    cols = all_of(gene_cols),
    names_to = "Gene",
    values_to = "CNV"
  )

cnv_long <- cnv_long[,c('ModelID','cluster','Gene','CNV')]

# load genomic data
cnv_long_2 <- merge(cnv_long, genomic_profile, by = "ModelID")

cnv_long_2$CNV_rel <- log2(cnv_long_2$CNV / cnv_long_2$Ploidy)
cnv_long_2$CNV_rel_2 <- log2(cnv_long_2$CNV / cnv_long_2$Ploidy + 1)

cnv_long_2 <- cnv_long_2 %>%
  mutate(CNV_class = case_when(
    CNV_rel >= 1       ~ "Amplification",
    CNV_rel >= 0.5      ~ "Gain",
    CNV_rel <= -1     ~ "Loss",
    CNV_rel <= -0.5     ~ "Loss",
    TRUE                ~ "Normal"
  ))

table(cnv_long_2$cluster,cnv_long_2$CNV_class)

cnv_gene_summary <- cnv_long_2 %>%
  group_by(cluster, Gene) %>%
  summarise(
    mean_CNV_abs = mean(CNV, na.rm = TRUE),
    mean_CNV_rel = mean(CNV_rel, na.rm = TRUE),
    mean_CNV_rel_2 = mean(CNV_rel_2, na.rm = TRUE)
  ) %>%
  ungroup()

# Filtramos Amplification y Gain, luego calculamos porcentaje
cnv_amp_summary <- cnv_long_2 %>%
  group_by(cluster, Gene) %>%
  summarise(
    total = n(), 
    n_amp = sum(CNV_class %in% c( "Amplification"))
  ) %>%
  mutate(perc_amp = n_amp / total * 100) %>%
  ungroup()

# Filtramos Amplification y Gain, luego calculamos porcentaje
cnv_loss_summary <- cnv_long_2 %>%
  group_by(cluster, Gene) %>%
  summarise(
    total = n(), 
    n_amp = sum(CNV_class %in% c( "Loss"))
  ) %>%
  mutate(perc_amp = n_amp / total * 100) %>%
  ungroup()



png('25Q3_RESULTS/test_2.png', width = 2000, height = 1200, res = 250)
ggplot(cnv_gene_summary,
       aes(x = cluster,
           y = mean_CNV_rel_2,
           fill = cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
  scale_fill_manual(values = cluster_colors) +
  scale_color_manual(values = cluster_colors) +
  #ylim(c(0,100))+
  #geom_text_repel(
   # data = cnv_gene_summary,
   #aes(label = Gene),
  #size = 2.5,
  #box.padding = 0.9,
  #point.padding = 0.2,
  #show.legend = FALSE
  #) +
  theme_minimal()
dev.off()

write.table(cnv_gene_summary, '25Q3_RESULTS/mean_cnv_by_cluster_rel.tsv', sep = '\t', row.names = F)
write.table(cnv_long_2, '25Q3_RESULTS/cnv_by_celline_rel.tsv', sep = '\t',row.names = F)

# ---------------- Mutations ---------------- 

# Contar mutaciones relevantes por cada línea celular (ModelID)
mutaciones_por_modelo <- mutaciones_relevants %>%
  group_by(ModelID) %>%
  summarise(num_mutaciones_relevantes = n())

# Longitud del exoma en megabases (Mb)
longitud_genoma_Mb <- 50  # Ajusta según corresponda

# Calcular el TMB por cada línea celular
mutaciones_por_modelo <- mutaciones_por_modelo %>%
  mutate(TMB = num_mutaciones_relevantes / longitud_genoma_Mb)


mutaciones_relevants <- merge(mutaciones_relevants,met_filt,by.x='ModelID',by.y='ModelID')

# Número de líneas celulares únicas por cluster
lineas_por_cluster <- mutaciones_relevants %>%
  distinct(ModelID, cluster) %>%
  group_by(cluster) %>%
  summarise(total_lineas = n())

# Contar cuántas líneas celulares en cada cluster tienen una mutación en cada gen
mutaciones_por_cluster <- mutaciones_relevants %>%
  distinct(ModelID, HugoSymbol, cluster) %>%
  group_by(cluster, HugoSymbol) %>%
  summarise(lineas_con_mutacion = n())

porcentaje_mutaciones <- mutaciones_por_cluster %>%
  left_join(lineas_por_cluster, by = "cluster") %>%
  mutate(porcentaje = (lineas_con_mutacion / total_lineas) * 100) %>%
  select(cluster, HugoSymbol, porcentaje)
write.table(porcentaje_mutaciones, '25Q3_RESULTS/percentage_mut_cluster.tsv', sep = '\t', row.names = F)
write.table(metadata_filt_pct, '25Q3_RESULTS/percentage_cluster_tumor.tsv', sep = '\t',row.names = F)

genes_interes <- c("TP53", "CDKN2A","RB1","KRAS")

png('25Q3_RESULTS/test.png', width = 2000, height = 1200, res = 250)
ggplot(porcentaje_mutaciones,
       aes(x = cluster,
           y = porcentaje,
           fill = cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
  scale_fill_manual(values = cluster_colors) +
  scale_color_manual(values = cluster_colors) +
  ylim(c(0,100))+
  #geom_text_repel(
  #  data = subset(porcentaje_mutaciones, HugoSymbol %in% genes_interes),
   # aes(label = HugoSymbol),
    #size = 2.5,
    #box.padding = 0.9,
    #point.padding = 0.2,
    #show.legend = FALSE
  #) +
  geom_hline(yintercept = 25,
             linetype = "dashed",
             color = "grey60",
             linewidth = 0.4) +
  theme_minimal()
dev.off()

genes_interes <- c("TP53", "CDKN2A","RB1")

df_plot <- porcentaje_mutaciones %>%
  filter(HugoSymbol %in% genes_interes)

gene_colors <- c(
  RB1     = "#739683",
  CDKN2A = "#B6CEC1",
  TP53   = "#31493D"
)

png('25Q3_RESULTS/percentage_mut.png',width=700, height=600, res = 160)
ggplot(df_plot, aes(x = cluster, y = porcentaje, fill = HugoSymbol)) +
  geom_col(position = position_dodge(width = 0.85)) +
  scale_fill_manual(values = gene_colors) +
  labs(
    x = "Cluster",
    y = "% cell lines mutated",
    fill = "Gene"
  ) +
  ylim(c(0,100))+
  theme_classic()
dev.off()


ggplot(df_plot, aes(x = cluster, y = porcentaje, fill = HugoSymbol)) +
  geom_col(position = position_dodge(width = 1)) +
  labs(
    x = "Cluster",
    y = "% cell lines mutated",
    fill = "Gene"
  ) +
  ylim(c(0,100))+
  theme_classic()


# Crear una matriz binaria donde las filas son genes y las columnas son líneas celulares

df_wide <- mutaciones_relevants %>%
  distinct(ModelID, HugoSymbol) %>%
  mutate(mutated = 1) %>%
  pivot_wider(names_from = ModelID, values_from = mutated, values_fill = 0)

# Asegurarte de que las columnas son de tipo carácter
mutaciones_relevants$ModelID <- as.character(mutaciones_relevants$ModelID)
mutaciones_relevants$HugoSymbol <- as.character(mutaciones_relevants$HugoSymbol)
mutaciones_relevants$VariantInfo <- as.character(mutaciones_relevants$VariantInfo)



# ---------------- Oncoprint ---------------- 


df_wide <- mutaciones_relevants %>%
  distinct(ModelID, HugoSymbol, VariantType) %>%
  pivot_wider(names_from = ModelID, values_from = VariantType, values_fill = list(""))

# Convertir el dataframe a matriz
mat <- as.matrix(df_wide[, -1])
rownames(mat) <- df_wide$HugoSymbol

# Crear el gráfico oncoprint


col = c("deletion" = "#1874CD", "insertion" = "#CD5555", "SNV" = "#CDAD00",'substitution' = '#2E8B57')

# just for demonstration
alter_fun = list(
  background = alter_graphic("rect", fill = "#CCCCCC"),   
  deletion = alter_graphic("rect", fill = col["deletion"]),
  insertion = alter_graphic("rect", fill = col["insertion"]),
  SNV = alter_graphic("rect",  fill = col["SNV"]), #height = 0.33,
  substitution = alter_graphic("rect",  fill = col["substitution"])#height = 0.33,
)


column_title = "Breast Cancer Cell Lines"
heatmap_legend_param = list(title = "Alternations", at = c("deletion", "insertion", "SNV","substitution"), 
                            labels = c("deletion", "insertion", "SNV","substitution"))

breast_cell_lines <- metadata[metadata$OncotreeLineage == 'Breast',]$ModelID
lung_cell_lines <- metadata[metadata$OncotreeLineage == 'Lung',]$ModelID

driver_genes <- read.csv('/Users/asanchezb/Desktop/PROJECTS/CDK_PREDICT_CDK4_6_BRCA/RESOURCES/breast-cancer-driver-genes.txt', sep = '\t',header = F)
driver_genes <- read.csv('/Users/asanchezb/Desktop/PROJECTS/CDK_PREDICT_CDK4_6_BRCA/RESOURCES/lung-cancer-driver-genes.txt', sep = '\t',header = F)

mat <- mat[,colnames(mat) %in% breast_cell_lines]
mat <- mat[,colnames(mat) %in% lung_cell_lines]
mat <- mat[,colnames(mat) %in% D_cell_lines]
mat <- mat[row.names(mat) %in% driver_genes$V1,]

#int_genes <- c('RB1','PTEN','CDKN2A','GATA3','MAP3K1','KRAS','BRAF','CCND1','KMT2C','ATM','BRCA1','BRCA2','TP53','PIK3CA','NF1','ARID1A','AR','HRAS','KMT2D','MEN1','PPP2R1A')

# Convierte la columna 'cluster' a un factor (opcional, pero a menudo útil para anotaciones)
clus <- factor(clusters[clusters$ModelID %in% breast_cell_lines,]$cluster)

# Crea la anotación de calor
cluster_colors <- c("1" = "#CD3333", "2" = "#1874CD", "5" = "#EE9572", "6" = "#7A378B")
cluster_colors <- c("1" = "#CD3333", "2" = "#1874CD", "3" = "#CD9B1D", "4" = "#CD6090", "5" = "#EE9572", "6" = "#7A378B")


# Filtra la información del cluster para las columnas de la matriz
annotation_data <- met_filt[met_filt$ModelID %in% colnames(mat), "cluster"]
annotation_data <- met_filt[met_filt$ModelID %in% colnames(mat), "RB1_status"]
tmb_breast <- mutaciones_por_modelo[mutaciones_por_modelo$ModelID %in% colnames(mat),]$TMB
MSI_breast <- genomic_profile[genomic_profile$ModelID %in% colnames(mat),'MSIScore']

# Contar cuántas muestras tienen alguna alteración por gen
alter_count <- rowSums(mat != "" & !is.na(mat))

# Seleccionar los top 100 genes
top100_genes <- names(sort(alter_count, decreasing = TRUE))[1:100]

# Subset de la matriz
mat_top100 <- mat[top100_genes, ]

tmb_col_fun <- colorRamp2(
  c(min(tmb_breast, na.rm = TRUE),
    max(tmb_breast, na.rm = TRUE)/2,
    max(tmb_breast, na.rm = TRUE)),
  c("white", "#D8BFD8", "#5D478B")  # lila → violeta → índigo
)

msi_col_fun <- colorRamp2(
  c(min(MSI_breast, na.rm = TRUE),
    max(MSI_breast, na.rm = TRUE)/2,
    max(MSI_breast, na.rm = TRUE)),
  c("white", "#B4EEB4", "#698B69")  # verde claro → medio → oscuro
)

top_annotation <- HeatmapAnnotation(
  cbar = anno_oncoprint_barplot(),
  
  Cluster = anno_simple(annotation_data,
                        col = cluster_colors,
                        na_col = "grey"),
  
  TMB = anno_simple(tmb_breast,
                    col = tmb_col_fun,
                    na_col = "grey90"),
  
  MSI = anno_simple(MSI_breast,
                    col = msi_col_fun,
                    na_col = "grey90")
)

# Crear la anotación de calor
top_annotation <- HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                                    Cluster = anno_simple(annotation_data, col = cluster_colors, na_col = "grey"),
                                    TMB = tmb_breast,
                                    MSI = MSI_breast # Opcional: color para valores NA
)


png('25Q3_RESULTS/breast_drivers.png',width=2000, height=750, res = 120)
#png('25Q3_RESULTS/breast_drivers.png',width=1100, height=1250, res = 120)
oncoPrint(
  mat,
  alter_fun = alter_fun, 
  col = col, 
  column_title = column_title, 
  heatmap_legend_param = heatmap_legend_param,
  remove_empty_columns = TRUE, 
  remove_empty_rows = TRUE,
  pct_side = "right", 
  row_names_side = "left",
  top_annotation = top_annotation,
  column_split = factor(met_filt$cluster[match(colnames(mat), met_filt$ModelID)])
)
dev.off()


# ---------------- Coocurrence ---------------- 


#Generate Extended Data Figure 4a
m <- as.data.frame(t(ifelse(mat == "", 0,1)))

m$Cluster <- annotation_data
#retain driver genes mutated in at least X samples
m <- m[,colSums(m)>=20]
m$Cluster<-factor(m$Cluster)

png('25Q3_RESULTS/breast_coocurrence.png',width=1000, height=700, res = 160)
upset(
  m, 
  nsets = (ncol(m)-1),  # Asumiendo que la última columna es 'Cluster'
  point.size = 2,
  mb.ratio = c(0.5, 0.5),
  nintersects = NA,
  sets.bar.color = "#56B4E9",
  matrix.color = "#474954",
  shade.color = "#7284A8",
  text.scale = 1.2,
  sets.x.label = "Number of Cases",
  queries = list(
    list(
      query = elements, 
      query.name = "1", 
      params = list(c("Cluster","1","2","5","6")), 
      color = "#CD3333", 
      active = TRUE
    ),
    list(
      query = elements, 
      query.name = "2", 
      params = list(c("Cluster","2","5","6")), 
      color = "#CD9B1D", 
      active = TRUE
    ),
    list(
      query = elements, 
      query.name = "5", 
      params = list(c("Cluster","5","6")), 
      color = "#EE9572", 
      active = TRUE
    ),
    list(
      query = elements, 
      query.name = "6", 
      params = list(c("Cluster","6")), 
      color = "#7A378B", 
      active = TRUE
    )
  ),
  query.legend = "top",
  order.by = "freq",
  empty.intersections = NULL
)


dev.off()



# ---------------- Association with clusters ---------------- 

# Mutations associated with Cluster6
rcb<-factor(m$Cluster,levels=c("A","B","C","D","E","F1","F2","F3","F4","F5","F6"),labels=c(0,0,0,1,0,0,0,0,0,0,0))
df.pcr<-character()  
for (g in colnames(m)[1:(ncol(m)-1)]){
  gene       <- factor(m[,g])
  model      <- glm(formula = rcb ~ gene, family = "binomial")
  lrm        <- summary(model)
  pval       <- lrm$coefficients[c(2:nrow(lrm$coefficients)),"Pr(>|z|)"]
  odds.ratio <- exp(lrm$coefficients[c(2:nrow(lrm$coefficients)),"Estimate"])
  conf       <- suppressMessages(confint(model))
  ci.low     <- exp(conf[c(2:nrow(conf)),1])
  ci.hi      <- exp(conf[c(2:nrow(conf)),2])
  md         <- data.frame(pval,odds.ratio,ci.low,ci.hi)
  md$variable<- g
  df.pcr <- rbind(df.pcr,md)
}
df.pcr$class<-"Cluster D"
write.table(df.pcr, '25Q3_RESULTS/Cluster_D_mut_association.tsv', sep = '\t')

mutAssoc <- df.pcr
mutAssoc$col <- "INSIG"
mutAssoc[mutAssoc$odds.ratio > 1 & mutAssoc$pval <= 0.05,"col"] <- "UP"
mutAssoc[mutAssoc$odds.ratio < 1 & mutAssoc$pval <= 0.05,"col"] <- "DOWN"
mutAssoc<-mutAssoc[order(mutAssoc$pval),] 
mutAssoc$variable<-factor(mutAssoc$variable,levels=unique(mutAssoc$variable))
mutAssoc$class<-factor(mutAssoc$class,levels=c("Cluster 6"))
mutAssoc$p_adj <- p.adjust(mutAssoc$pval, method = "BH")
mutAssoc<-mutAssoc[mutAssoc$pval <= 0.05,] 
mutAssoc$variable <- factor(mutAssoc$variable,
                        levels = mutAssoc$variable[order(mutAssoc$odds.ratio,
                                                     decreasing = TRUE)])
#mutAssoc<-mutAssoc[abs(mutAssoc$odds.ratio) > 1 & mutAssoc$pval <= 0.05,] 

top10 <- mutAssoc %>%
  filter(pval <= 0.05) %>%
  group_by(col) %>%
  slice_min(pval, n = 5, with_ties = FALSE) %>%
  ungroup()

top10$variable <- factor(top10$variable,
                            levels = top10$variable[order(top10$odds.ratio,
                                                             decreasing = TRUE)])

png('25Q3_RESULTS/log_odds_ratio_pancancer_clusterD.png',width=1000, height=600, res = 200)
eFig4a_inset <-
  ggplot(data=top10, aes(y=variable, x=log(odds.ratio))) + 
  geom_errorbarh(aes(xmax = log(ci.hi), xmin = log(ci.low),color=col), size = .5, height = .2) +
  geom_point(aes(color=col),shape=16,size=3)+
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dotted") +
  labs (x="Log odds ratio", y="") + 
  scale_color_manual(values=c("#7F7F7F","gray","#7A378B"),breaks = c("DOWN","INSIG","UP"))+
  scale_x_continuous(limits=c(-10,10))+
  guides(color="none",size="none")+
  #facet_grid(~class,scales="free_y",space = "free")+
  theme_manuscript(base_size = figure_font_size)+
  theme(legend.position = "bottom",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line = element_blank(),
        strip.background =element_blank(),
        strip.text = element_text(size = 12),
        axis.text.y = element_text(face="italic"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.3,0.3,0.3,0.3), "lines"))
eFig4a_inset
dev.off()




##################################
# 7. RB pathway
##################################


df_plot <- met_filt %>%
  group_by(Group) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(percentage = n / sum(n) * 100)

df_plot <- met_filt %>%
  group_by(cluster, Group) %>%   # Agrupamos por linaje y grupo
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster) %>%          # Normalizamos dentro de cada linaje
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()

group_colors <- c(
  "RBmut" = "#724F13",
  "RBupstream" = "#D59D3B",
  "WT" = "#E8C584"
)

# Filtrar WT
df_plot_noWT <- df_plot %>%
  filter(Group != "WT")
df_plot_noWT_complete <- df_plot_noWT %>%
  tidyr::complete(cluster, Group, fill = list(percentage = 0))

df_plot_filtered <- df_plot_noWT_complete %>%
  group_by(cluster) %>%
  filter(sum(n) >= 5) %>%  # suma de cell lines por cluster >= 5
  ungroup()

png('25Q3_RESULTS/test_8.png',
    width = 700, height = 1000, res = 250)

ggplot(df_plot_filtered, aes(x = Group, y = percentage, fill = Group)) + 
  geom_bar(stat = "identity",position = position_dodge(width = 0.99)) + 
  theme_minimal() +
  ylim(c(0,100.5)) +
  labs(x = "Group", y = "Percentage", fill = "Group") +
  scale_fill_manual(values = group_colors)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        plot.margin = margin(10, 10, 10, 10, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))

ggplot(df_plot_filtered, aes(x = Group, y = percentage, fill = Group)) + 
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.9),  # aumenta la separación
           width = 0.9) +                           # hace las barras más delgadas
  theme_minimal() +
  ylim(c(0,100.5)) +
  labs(x = "Cluster", y = "Percentaje", fill = "Group") +
  scale_fill_manual(values = group_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        panel.grid = element_blank(),
        plot.margin = margin(10, 10, 10, 10, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))

dev.off()

# by subtype
met_filt_breast <- met_filt %>% filter(OncotreeLineage == "Breast")
met_filt_breast <- met_filt %>%filter(OncotreeLineage == "Lung")

table(met_filt_breast$OncotreeCode)
table(met_filt_breast$DepmapModelType)
table(met_filt_breast$OncotreePrimaryDisease)


met_filt_breast <- met_filt_breast %>%
  mutate(Subtype = case_when(
    OncotreePrimaryDisease %in% c('Lung Neuroendocrine Tumor') ~ 'SCLC',
    OncotreePrimaryDisease %in% c('Non-Small Cell Lung Cancer') ~ 'NSCLC',
    TRUE ~ OncotreePrimaryDisease  # cualquier otro valor se mantiene
  ))

table(met_filt_breast$Subtype)
table(met_filt_breast$ModelSubtypeFeatures.x)

df_plot <- met_filt_breast %>%
  group_by(Subtype, Group) %>%   # Agrupamos por linaje y grupo
  summarise(n = n(), .groups = "drop") %>%
  group_by(Subtype) %>%          # Normalizamos dentro de cada linaje
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()

# Filtrar WT
df_plot_noWT <- df_plot %>%
  filter(Group != "WT")%>%
  filter(Subtype != "NA")
df_plot_noWT_complete <- df_plot_noWT %>%
  tidyr::complete(Subtype, Group, fill = list(percentage = 0))

df_plot_filtered <- df_plot_noWT_complete %>%
  group_by(Subtype) %>%
  filter(sum(n) >= 5) %>%  # suma de cell lines por cluster >= 5
  ungroup()


png('25Q3_RESULTS/test_8.png',
    width = 800, height = 1200, res = 250)

ggplot(df_plot_filtered, aes(x = Subtype, y = percentage, fill = Group)) + 
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.9),  # aumenta la separación
           width = 0.9) +                           # hace las barras más delgadas
  theme_minimal() +
  ylim(c(0,100.5)) +
  labs(x = "Cluster", y = "Percentaje", fill = "Group") +
  scale_fill_manual(values = group_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        panel.grid = element_blank(),
        plot.margin = margin(10, 10, 10, 10, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))

dev.off()


##################################
# 8. RECLUSTERING BREAST
##################################

# ---------------- Re-Clustering selected ---------------- 

rownames(merged_sel) <- met_filt$ModelID
merged_new <- merged_sel[rownames(merged_sel) %in% met_filt[met_filt$OncotreeLineage == 'Breast',]$ModelID,]
#merged_new <- merged_sel[rownames(merged_sel) %in% met_filt[met_filt$cluster %in% c(4,5),]$ModelID,]


# Determine the optimal number of clusters
png('25Q3_RESULTS/silhouette_selected_3.png',width=1200, height=1000, res = 250)
fviz_nbclust(merged_new, kmeans, method = "silhouette") + 
  labs(subtitle = "silhouette Method")
dev.off()  

# Determine the optimal number of clusters
png('25Q3_RESULTS/gap_stat_selected_3.png',width=1200, height=1000, res = 250)
fviz_nbclust(merged_new, kmeans, method = "gap_stat") + 
  labs(subtitle = "gap_stat Method")
dev.off()  

# Determine the optimal number of clusters
png('25Q3_RESULTS/wss_selected_3.png',width=1200, height=1000, res = 250)
fviz_nbclust(merged_new, kmeans, method = "wss") + 
  labs(subtitle = "wss Method")
dev.off()  

# seed the clustering to be reproducible
set.seed(124)
# number of clusters selected
k <- 5
# kmeans clustering
km_res <- kmeans(merged_new, centers = k, nstart = 25)
# Append cluster assignment to the data
cluster_assignment <- km_res$cluster
cluster_colors <- c("A" = "#CD3333", "B" = "#1874CD", "C" = "#CD9B1D", "D" = "#7A378B", "E" = "#EE9572", 
                    "F1" = "#17301D","F2" = "#264F30", "F3" = "#3E804F","F4" = "#72BB84", "F5" = "#B0D9BA", "F6" = "#CFE8D5")
cluster_colors2 <- c("1" = "#CD3333", "2" = "#1874CD", "3" = "#CD9B1D", "4" = "#CD6090", "5" = "#EE9572", "6" = "#7A378B")
# Add to metadata
met_new <- met_filt[met_filt$OncotreeLineage == 'Breast',]
#met_new <- met_filt[met_filt$cluster %in% c(4,5), ]
met_new$cluster2 <- cluster_assignment
table(met_new$cluster)
#write.table(met_filt[,c('ModelID','cluster')], '25Q3_RESULTS/cluster_assigment.csv', sep = ',', row.names = F)
# Colors for mutations 
unique_lineages <- unique(met_new$RB1_status)
color_palette <- c( "#EEE9E9","#424242")  
mut_colors <- setNames(color_palette[seq_along(unique_lineages)], unique_lineages)


met_new <- met_new %>%
  mutate(Subtype = case_when(
    OncotreePrimaryDisease %in% c('Lung Neuroendocrine Tumor') ~ 'SCLC',
    OncotreePrimaryDisease %in% c('Non-Small Cell Lung Cancer') ~ 'NSCLC',
    OncotreePrimaryDisease %in% c('SMARCA4-deficient undifferentiated tumor') ~ 'Other',
    TRUE ~ OncotreePrimaryDisease  # cualquier otro valor se mantiene
  ))

# Colors for lineage
unique_lineages <- unique(met_new$Subtype)
color_palette <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique_lineages))
lineage_colors <- setNames(color_palette[seq_along(unique_lineages)], unique_lineages)

# Asegurarse de que los niveles coincidan y reemplazar NA por "Unknown" o similar
met_new$Subtype <- as.character(met_new$Subtype)
met_new$Subtype[is.na(met_new$Subtype)] <- "Unknown"

# Definir colores con nombres correctos
lineage_colors <- c(
  "HER2+" = "#B68FAB",
  "ER+" = "#CFAEB5",
  "TNBC" = "#513671",
  "ER+ HER2+" = "#814B95",
  "Unknown" = "#F9F2E8"
)

# Definir colores con nombres correctos
lineage_colors <- c(
"Other" = "#814B95",
  "NSCLC" = "#CFAEB5",
  "SCLC" = "#513671"
)

# Colors for dependency
col_fun <- colorRamp2(c(-3, 0, 1), c("#A90C38", "white", "#0B0B8C")) #"#0B0B8C", "white", "#A90C38"

# Plot heatmap
png('25Q3_RESULTS/Pancancer_CDKs_Cyclins_sel_other_6.png',width=1700, height=800, res = 150) #width=2000, height=500, res = 150
draw(
  Heatmap(
    t(merged_new), 
    name = "Gene effect",
    show_row_names = TRUE, 
    show_column_names = FALSE, 
    column_names_rot = 55,
    column_names_gp = gpar(fontsize = 8.5),
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    #column_split = met_new$cluster,
    column_split = cluster_assignment,
    row_title = "Genes",
    column_title = "Cell Lines",
    top_annotation = HeatmapAnnotation(
      Cluster2 = as.factor(cluster_assignment),
      Cluster = met_new$cluster,
      Subtype = met_new$Subtype,
      Group = met_new$Group,
      TP53 = met_new$TP53_status,
      RB1 = met_new$RB1_status,
      WNT_pathway = met_new$WNT_pathway,
      PI3K_pathway = met_new$PI3K_pathway,
      RB1_upstream = met_new$Rb_pathway,
      RBL1_2_status = met_new$RBL1_2_status,
      col = list(
        Cluster = cluster_colors,
        Cluster2 = cluster_colors2,
        Subtype = lineage_colors,
        Group = group_colors,
        WNT_pathway = mut_colors,
        TP53 = mut_colors,
        RB1 = mut_colors,
        PI3K_pathway = mut_colors,
        RB1_upstream = mut_colors,
        RBL1_2_status = mut_colors
      )
    ),
    col = col_fun
  ),
  #show_heatmap_legend = FALSE,
  show_annotation_legend = TRUE
)
dev.off()

# see expression

ht <- Heatmap(
  t(merged_new),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_split = met_new$cluster,
  col = col_fun
)

ht_drawn <- draw(ht)
col_order_heatmap <- unlist(column_order(ht_drawn))
ordered_samples <- rownames(merged_new)[col_order_heatmap]

sel <- c('CDKN2A', 'CDKN2B', 'RB1', 'RBL1', 'RBL2', 'CCNE1', 'CCNE2', 'E2F1', 'E2F2', 'E2F3')
common_samples <- intersect(ordered_samples, rownames(expression_data))
new_df_ordered <- expression_data[common_samples,sel]

# Colors for dependency
col_fun <- colorRamp2(c(-2,0,2), c("#EEDC82","white","#5D478B")) #"#0B0B8C", "white", "#A90C38"

# Colors for mutations 
unique_lineages <- unique(met_new$RB1_loss)
color_palette <- c( "#EEE9E9","#424242")  
cnv_colors <- setNames(color_palette[seq_along(unique_lineages)], unique_lineages)

# Add to metadata
met_exp <- met_filt[match(rownames(new_df_ordered), met_filt$ModelID), ]
png('25Q3_RESULTS/Pancancer_CDKs_Cyclins_sel_other_6.png',width=1700, height=600, res = 150) #width=2000, height=500, res = 150
draw(
  Heatmap(
    t(scale(t(t(new_df_ordered)))), 
    name = "Gene effect",
    show_row_names = TRUE, 
    show_column_names = FALSE, 
    column_names_rot = 55,
    column_names_gp = gpar(fontsize = 8.5),
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    column_split = met_exp$cluster,
    #column_split = cluster_assignment,
    row_title = "Genes",
    column_title = "Cell Lines",
    top_annotation = HeatmapAnnotation(
      #Cluster2 = as.factor(cluster_assignment),
      Cluster = met_exp$cluster,
      Subtype = met_exp$Subtype,
      Group = met_exp$Group,
      TP53 = met_exp$TP53_status,
      RB1mut = met_exp$RB1_status,
      RB1loss = met_exp$RB1_loss,
      WNT_pathway = met_exp$WNT_pathway,
      PI3K_pathway = met_exp$PI3K_pathway,
      RB1_upstream = met_exp$Rb_pathway,
      RBL1_2_status = met_exp$RBL1_2_status,
      col = list(
        Cluster = cluster_colors,
        #Cluster2 = cluster_colors2,
        Subtype = lineage_colors,
        Group = group_colors,
        WNT_pathway = mut_colors,
        TP53 = mut_colors,
        RB1mut = mut_colors,
        RB1loss = cnv_colors,
        PI3K_pathway = mut_colors,
        RB1_upstream = mut_colors,
        RBL1_2_status = mut_colors
      )
    ),
    col = col_fun
  ),
  #show_heatmap_legend = FALSE,
  show_annotation_legend = TRUE
)
dev.off()

##################################
# 9. MARKERS BY CLUSTERS
##################################

cts <- t(expression_data) 
coldata <- met_filt[met_filt$ModelID %in% colnames(cts),]

names <- coldata$ModelID
rownames(coldata) <- names
cts <- cts[,colnames(cts) %in% rownames(coldata)]
all(rownames(coldata) %in% colnames(cts))

# order data
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

coldata$cluster <- factor(
  coldata$cluster,
  levels = c("A","B","C","D","E","F1","F2","F3","F4","F5","F6")
)

design <- model.matrix(~0 + cluster, data = coldata)
colnames(design) <- levels(coldata$cluster)

fit <- lmFit(cts, design)

cats <- levels(coldata$cluster)

contrast_list <- sapply(cats, function(ct){
  others <- setdiff(cats, ct)
  paste0(ct, " - (", paste(others, collapse="+"), ")/", length(others))
})

names(contrast_list) <- paste0(cats, "_vs_rest")

contrast.matrix <- makeContrasts(
  contrasts = contrast_list,
  levels = design
)
colnames(contrast.matrix) <- paste0(levels(coldata$cluster), "_vs_rest")
colnames(contrast.matrix)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# logFC matrix
logFC_mat <- fit2$coefficients

# p-values
pval_mat <- fit2$p.value

# adjusted p-values (FDR)
padj_mat <- p.adjust(pval_mat, method="BH")
padj_mat <- matrix(padj_mat,
                   nrow = nrow(pval_mat),
                   dimnames = dimnames(pval_mat))
score_mat <- logFC_mat * -log10(padj_mat)

genes <- rownames(logFC_mat)

final_table <- data.frame(gene = genes)

for (ct in cats) {
  
  coef_name <- paste0(ct, "_vs_rest")
  
  final_table[[paste0(ct, "_FC")]]   <- logFC_mat[, coef_name]
  final_table[[paste0(ct, "_pval")]] <- pval_mat[, coef_name]
  final_table[[paste0(ct, "_padj")]] <- padj_mat[, coef_name]
  final_table[[paste0(ct, "_score")]] <- score_mat[, coef_name]
}

write.table(final_table,
            "25Q3_RESULTS/cluster_vs_rest_wide_table.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# select top markers
# ejemplo: top 50 genes positivos por cluster
top_markers <- lapply(cats, function(ct){

  tmp <- final_table[order(-final_table[[paste0(ct, "_score")]]), ]
  
  tmp <- tmp[tmp[[paste0(ct, "_FC")]] > 0, ]
  
  head(tmp$gene, 100)
})

names(top_markers) <- cats

# cats ya contiene los nombres de los clusters
# top_markers ya está definido como lista de vectores de genes

# convertir lista a data.frame ancho

max_genes <- max(sapply(top_markers, length))  # normalmente 50

# Crear tabla: filas = Rank, columnas = Clusters
marker_table <- data.frame(matrix(NA, nrow = max_genes, ncol = length(top_markers)))
rownames(marker_table) <- paste0("Rank_", 1:max_genes)
colnames(marker_table) <- names(top_markers)

for(ct in names(top_markers)){
  genes <- top_markers[[ct]]
  marker_table[1:length(genes), ct] <- genes
}

# Guardar en archivo
write.table(marker_table,
            file = "25Q3_RESULTS/top100_markers_per_cluster.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# Guardar en archivo
write.table(met_filt,
            file = "25Q3_RESULTS/metadata_200226.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)




