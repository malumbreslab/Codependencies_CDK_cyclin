# Manuscript:   A reference map of CDK-cyclin dependencies across human cancers
# Project:      Codependencies
# Part:         B. PRISM screening
# Author:       Agustín Sánchez-Belmonte
# Date:         06/02/24
# Last update:  04/02/26

# ---------------- INDEX ---------------- 

# 0. SET UP
# 1. LOAD DATA
#   a. PRISM data
#   b. Metadata
#   c. Dependency
#   d. Mutations
#   e. Merge and filter
# 2. CORRELATION ANALYSIS
# 3. PRIMS PROFILE TUMOR TYPE
# 4. CORRELATION KO AND PRISM
# 5. CORRELATION EXPRESSION AND PRISM
# 6. AUC VALUES IN CLUSTERS


#############################################################################################################
# 0. SET UP
#############################################################################################################

# ---------------- PATHS ---------------- 

setwd("CODEPENDENCIES_CDKS_CYCLINS/")
DATA="CODEPENDENCIES_CDKS_CYCLINS/19Q4_data"

# ---------------- PACKAGES ---------------- 

library(readxl)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(tidyr)
library(RColorBrewer)
library(dplyr)
library(tibble)



#############################################################################################################
# 1. LOAD DATA
#############################################################################################################

# PRISM Data were downloaded from DepMap 19Q4 (DepMap, 2019).

# ---------------- PRISM SCREEN DATA ---------------- 


#load drug screen data & sample info --> download from https://depmap.org/repurposing/
drug_screen<-read.csv(paste0(DATA,"/secondary-screen-dose-response-curve-parameters.csv"),stringsAsFactors = F)

# raw data of logfold change to different doses
raw_drug_screen<-read.csv('PRISM_SCREENING/data/secondary-screen-replicate-collapsed-logfold-change.csv',stringsAsFactors = F)
names_cell <- raw_drug_screen$X

# Data were downloaded from DepMap 25Q3 (DepMap, 2025).

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


# ---------------- METADATA ---------------- 
cell_line_info <- read.csv('25Q3_data/Model.csv', check.names = F)
cell_line_info <- read.csv('/Users/asanchezb/Desktop/metadata_dep_200226.tsv', check.names = F, sep = '\t')

# ---------------- DEPENDENCY DATA ---------------- 
dependency_data <- read.csv('25Q3_data/CRISPRGeneEffect.csv', check.names = F)

names <- c()
for(i in colnames(dependency_data)){
  new <- strsplit(i, " ")[[1]][1]
  names <- c(names,new)
}

colnames(dependency_data) <- names
names(dependency_data)[1] <- 'ModelID'

ALL <- c('ModelID','CDK1','CDK2','CDK3','CDK4','CDK5','CDK6','CDK7','CDK8','CDK9','CDK10','CDK12','CDK13','CDK14','CDK15','CDK16','CDK17','CDK18','CDK19','CDK20',"CCNA1", "CCNA2", "CCNB1", "CCNB2", "CCNB3", "CCNC", "CCND1", "CCND2", "CCND3", "CCNE1", "CCNE2", "CCNF","CCNG1", "CCNG2", "CCNH", "CCNI", "CCNI2", "CCNJ", "CCNJL", "CCNK", "CCNL1", "CCNL2", "CCNO", "CCNP","CCNT1", "CCNT2", "CCNY", "CCNYL1","CNPPD1", "CNTD1", "CABLES1", "CABLES2","CDK5R1", "CDK5R2","SPDYA", "SPDYE1","SPDYE4")#  "CCNQ", "SPDYE3"
dep_CDKs <- dependency_data[,ALL]

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

# ---------------- MERGE AND FILTER DATA ---------------- 

drug_screen <- merge(drug_screen, cell_line_info,by.x='depmap_id',by.y='ModelID', all.x = TRUE)

met_filt_drug <- merge(met_filt, 
                     drug_screen_CDKS_wide,
                     by.x = "ModelID",
                     by.y = "depmap_id",
                     all.x = TRUE)


# select drug
table(drug_screen$name)
table(drug_screen$screen_id)

# CDK inhibitors
CDK_inh=unique(drug_screen$name[grep("CDK",drug_screen$target)])
drug_screen_CDKS <-drug_screen[drug_screen$name %in% CDK_inh,]
table(drug_screen_CDKS$screen_id) # HTS002 MTS010
table(drug_screen_CDKS$target)

drug_screen_CDKS <- drug_screen_CDKS[drug_screen_CDKS$screen_id == 'MTS010',]
drug_screen_CDKS <- drug_screen_CDKS[,c('depmap_id','name','auc')]

drug_screen_CDKS <- drug_screen_CDKS %>%
  filter(!is.na(depmap_id))

drug_screen_CDKS_wide <- drug_screen_CDKS %>%
  pivot_wider(names_from = name, values_from = auc)

min(drug_screen_CDKS$auc)
max(drug_screen_CDKS$auc)

# merge 
NAMES <-dep_CDKs$ModelID
dependencia_genes_adjusted <- dep_CDKs[,-c(1)]
row.names(dependencia_genes_adjusted) <- NAMES

merge_dep_drug <- merge(dependencia_genes_adjusted,drug_screen_CDKS_wide, by.x = 0, by.y = 'depmap_id')
NAMES2 <- merge_dep_drug$Row.names
merge_dep_drug <- merge_dep_drug[,-c(1)]
row.names(merge_dep_drug) <- NAMES2


#############################################################################################################
# 2. CORRELATION
#############################################################################################################

# ---------------- CORRELATION ---------------- 

# Calculate corr and pvalues
res <- cor(merge_dep_drug, method = "spearman",use = "pairwise.complete.obs")
#res <- cor(merge_dep_drug, use = "pairwise.complete.obs")
testRes = cor.mtest(merge_dep_drug, conf.level = 0.95, method = "spearman")

nCDKs <- length(ALL)-1
nCDKs_2 <- nCDKs+1

total_genes = dim(res)[2]
no_cdks_genes = total_genes - nCDKs
# Pivot tables
res = as.data.frame(res[1:nCDKs,nCDKs_2:total_genes]) #53
res$CDK <- rownames(res)
new_res <- pivot_longer(res,cols = 1:no_cdks_genes, names_to = 'DRUG', values_to = 'Corr')

pval = as.data.frame( testRes$p[1:nCDKs,nCDKs_2:total_genes])
pval$CDK <- rownames(pval)
new_pval <- pivot_longer(pval,cols = 1:no_cdks_genes, names_to = 'DRUG', values_to = 'Pvalue')

# merge
new_res$Pvalue <- new_pval$Pvalue

new_res$CDK <- factor(new_res$CDK, levels=rev(rownames(res)))
new_res$DRUG <- factor(new_res$DRUG, levels=colnames(res))

dependency_means <- as.data.frame(colMeans(merge_dep_drug[sapply(merge_dep_drug, is.numeric)], na.rm=TRUE))
colnames(dependency_means) <- 'Dep'
dependency_means$CDK <- rownames(dependency_means)
dependency_means <- dependency_means[1:nCDKs,]
new_res_2 <- merge(new_res,dependency_means, by.x = 'CDK', by.y = 'CDK')


new_res_2$Pval <- ifelse(new_res_2$Pvalue < 0.05 , '< 0.05','No sign')

dep = 'Gene'

# ---------------- SIMILARITY ---------------- 


# 1. matriz DRUG x CDK
mat <- new_res_2 %>%
  select(DRUG, CDK, Corr) %>%
  pivot_wider(names_from = CDK, values_from = Corr) %>%
  column_to_rownames("DRUG")%>%
  as.matrix()

# Clustering jerárquico
hc_corr <- hclust(as.dist(1 - cor(t(mat), use = "pairwise.complete.obs")))

# Plot
png(filename = '/Users/asanchezb/Desktop/test.png', width = 1500, height = 1100, res = 200)
plot(hc_corr)
plot(hc_corr, labels = FALSE, main = "Dendrograma de DRUGS", xlab = "", ylab = "1 - Correlación")
dev.off()

# 3. orden
drug_order <- rownames(mat)[hc_corr$order]

# 4. aplicar orden
new_res_2$DRUG <- factor(new_res_2$DRUG, levels = drug_order)


# ---------------- PLOTTING ---------------- 

png(filename = '/Users/asanchezb/Desktop/test.png', width = 2500, height = 3200, res = 305)#1300#2500
p <- ggplot(new_res_2, aes(DRUG, CDK))+
  geom_point(aes(colour = Pval, size=Dep*-1, shape=Pval, fill=Corr)) + # -log10(Pvalue) # size=Dep*-1
  scale_size(paste0(dep, ' dep'), limits = c(-0.3,3))+ #limits = c(-0.3,3)
  theme_bw() +
  scale_color_manual(values=c("black","#E0E0E0"))+    
  scale_fill_gradientn(colours = rev(brewer.pal(name = 'RdBu', n=11)), limits = c(-0.5,0.5), oob = scales::squish, name = 'Correlation') +
  scale_x_discrete(position = "top") +
  cowplot::theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  #theme(axis.ticks = element_blank()) +
  #geom_text(label=new_res_2$sign, nudge_y = -0.05, size = 2.5) +
  #scale_shape_manual(name = "Pvalue", values = unique(new_res_2$sign)) 
  scale_shape_manual(values=c(21,21))
print(plot(p))
dev.off()


png('/Users/asanchezb/Desktop/corr_test_2.png', width = 1000, height = 1000, res = 250)
ggscatter(merge_dep_drug, x = "CCNE1", y = "palbociclib", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "CDK2", ylab = "palbociclib") 
dev.off()



drug_screen_CDKS_wide[drug_screen_CDKS_wide == 'NULL'] <- 'NA'
is.na(drug_screen_CDKS_wide)
drug_screen_abema <-drug_screen[drug_screen$name =="abemaciclib",]
drug_screen_abema <-drug_screen_abema[drug_screen_abema$screen_id =="MTS010",]
table(drug_screen_abema$screen_id)


#############################################################################################################
# 3. PRIMS SCREEN PROFILE
#############################################################################################################

# raw data of logfold change to different doses
raw_drug_screen<-read.csv('PRISM_SCREENING/data/secondary-screen-replicate-collapsed-logfold-change.csv',stringsAsFactors = F)
names_cell <- raw_drug_screen$X
# differetn scrrening 'screen_id'
# BRD-K51313569-001-07-8 - MTS010 palbo
# BRD-K51313569-003-03-3 - HTS002 palbo

# BRD-K33622447-066-01-9 - MTS010/HTS002 abema

raw_drug_screen <- raw_drug_screen[, grepl("BRD.K51313569.003.03.3", colnames(raw_drug_screen))] # palbociclib HTS002
raw_drug_screen <- raw_drug_screen[, grepl("BRD.K51313569.001.07.8", colnames(raw_drug_screen))] # palbociclib MTS010

raw_drug_screen <- raw_drug_screen[, grepl("BRD.K33622447.066.01.9", colnames(raw_drug_screen))] # abemaciclib 
raw_drug_screen <- raw_drug_screen[, grepl("MTS010", colnames(raw_drug_screen))] # abemaciclib HTS002

raw_drug_screen <- raw_drug_screen[, grepl("BRD.K36788280.001.01.2", colnames(raw_drug_screen))] # ribociclib 
raw_drug_screen <- raw_drug_screen[, grepl("MTS010", colnames(raw_drug_screen))] # ribociclib HTS002

rownames(raw_drug_screen) <- names_cell

# na remove data
raw_drug_screen <- na.omit(raw_drug_screen)

colnames(raw_drug_screen) <- c('0.0006 uM', '0.002 uM','0.01 uM','0.04 uM','0.15 uM', '0.6 uM','10 uM','2.5 uM')
colnames(raw_drug_screen) <- c('C0.0006_uM', 'C0.002_uM','C0.01_uM','C0.04_uM','C0.15_uM', 'C0.6_uM','C10_uM','C2.5_uM')
raw_drug_screen <- raw_drug_screen[,c('C0.0006_uM','C0.002_uM','C0.01_uM','C0.04_uM','C0.15_uM', 'C0.6_uM','C2.5_uM','1C0_uM')]

#
#drug_screen_abema <- drug_screen_abema[drug_screen_abema$depmap_id %in% rownames(raw_drug_screen),]
#raw_drug_screen <- raw_drug_screen[rownames(raw_drug_screen) %in% drug_screen_abema$depmap_id ,]

# Obtener el orden basado en la columna deseada de df2
orden <- match(drug_screen_abema$depmap_id, row.names(raw_drug_screen))

# Reordenar df1
drug_screen_abema <- drug_screen_abema[order(orden), ]
drug_screen_abema$depmap_id

# scale data between 0 - 1
#matriz_datos_escala <- (as.matrix(raw_drug_screen) - min(as.matrix(raw_drug_screen))) / (max(as.matrix(raw_drug_screen)) - min(as.matrix(raw_drug_screen)))
matriz_datos_escala <- apply(raw_drug_screen, 2, function(x) (x - min(raw_drug_screen, na.rm = TRUE)) / (max(raw_drug_screen, na.rm = TRUE) - min(raw_drug_screen, na.rm = TRUE)))

#rownames(drug_screen_abema) <- drug_screen_abema$depmap_id

# heatmap: distance and clustering method
png('/Users/asanchezb/Desktop/test.png', width = 500, height = 1100, res = 170)
pheatmap(matriz_datos_escala, 
         cluster_cols = FALSE, 
         show_dendrogram = FALSE,
         #color = colorRampPalette(c("#2E5A87", "#FFF6F4", "#CD2626"))(100), 
         color = colorRampPalette(c("#CD2626", "#FFF6F4", "#2E5A87"))(100), 
         clustering_distance_rows = 'manhattan', clustering_method = "ward.D",
         show_rownames = FALSE)

dev.off()

#############################################################################################################
# 3. PRIMS PROFILE TUMOR TYPE
#############################################################################################################


# dependency data of cell lines which are in PRIMS screening
cell_exp_filt_drug<-dependency_data[rownames(dependency_data_new)%in%drug_screen$depmap_id,]

# merge info with drug screen
merge_drug_info <- merge(drug_screen, cell_line_info, by.x='depmap_id',by.y='ModelID')
merge_drug_info <- drug_screen
#merge_drug_info <- merge_drug_info[merge_drug_info$screen_id == 'HTS002',]
merge_drug_info <- merge_drug_info[merge_drug_info$screen_id == 'MTS010',]

# Select drug
merge_drug_info_palbo <- merge_drug_info[merge_drug_info$name == 'palbociclib',]
merge_drug_info_palbo <- merge_drug_info[merge_drug_info$name == 'ribociclib',]
merge_drug_info_palbo <- merge_drug_info[merge_drug_info$name == 'abemaciclib',]



# Calcular el cuartil inferior
q25 <- quantile(merge_drug_info_palbo$auc, 0.5, na.rm = TRUE)
q75 <- quantile(merge_drug_info_palbo$auc, 0.6, na.rm = TRUE)
# Calcular medianas 
medianas <- merge_drug_info_palbo %>%
  group_by(OncotreeLineage) %>%
  summarize(mediana_auc = median(auc))

medianas <- medianas[order(medianas$mediana_auc,decreasing = T),]
merge_drug_info_palbo <- merge(merge_drug_info_palbo, medianas, by.x = 'OncotreeLineage',by.y='OncotreeLineage', all.x=TRUE)
#merge_drug_info_palbo$label <- ifelse(merge_drug_info_palbo$mediana_auc > lower_quartile, 'Resistant','Sensitive')

merge_drug_info_palbo$label <- dplyr::case_when(
  merge_drug_info_palbo$mediana_auc < q25 ~ "Sensitive",
  merge_drug_info_palbo$mediana_auc > q75 ~ "Resistant",
  TRUE                            ~ "Non effect"
)
merge_drug_info_palbo$label <- factor(
  merge_drug_info_palbo$label,
  levels = c("Sensitive", "Non effect","Resistant")
)


merge_drug_info_palbo$OncotreeLineage <- factor(merge_drug_info_palbo$OncotreeLineage, levels = medianas$OncotreeLineage)

merge_drug_info_palbo_filt <- merge_drug_info_palbo %>%
  group_by(OncotreeLineage) %>%
  filter(n() >= 5) %>%
  ungroup()

merge_drug_info_palbo_filt$OncotreeLineage <-
  droplevels(merge_drug_info_palbo_filt$OncotreeLineage)

merge_drug_info_palbo_filt$OncotreeLineage <- factor(
  merge_drug_info_palbo_filt$OncotreeLineage,
  levels = medianas$OncotreeLineage
)

table(merge_drug_info_palbo_filt$OncotreeLineage)
png('/Users/asanchezb/Desktop/test.png', width = 1400, height = 1000, res = 150)

ggplot(merge_drug_info_palbo_filt,
       aes(x = OncotreeLineage,
           y = auc,
           fill = label,
           group = OncotreeLineage)) +
  
  geom_boxplot(position = position_dodge(width = 0.8),
               width = 0.6) +
  
  geom_jitter(position = position_jitterdodge(
    jitter.width = 0.1,
    dodge.width = 0.8),
    size = 0.5)+
  xlab('') + 
  ylab('Palbociclib (AUC)') +
  scale_fill_manual(
    values = c(
      "Sensitive"  = "#E99A98",
      "Non effect" = "#FFFFF0",
      "Resistant"  = "#6A81A0"
    )
  ) +
  geom_hline(yintercept = q25, col = "black", linetype = "dashed") +
  geom_hline(yintercept = q75, col = "black", linetype = "dashed") +
  theme_bw() +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour = 'black'),
    panel.grid = element_blank(),
    plot.margin = margin(30, 30, 30, 30, "pt"),
    axis.line = element_line(color = "black")
  ) +
  #geom_boxplot(width = 0.1) +
  guides(fill = FALSE)


dev.off()


#############################################################################################################
# 4. CORRELATION KO AND PRISM
#############################################################################################################

# use MTS010

table(drug_screen_CDKS$screen_id) # HTS002 MTS010
drug_screen_CDKS <-drug_screen[drug_screen$name %in% CDK_inh,]

drug_screen_CDKS <- drug_screen_CDKS[drug_screen_CDKS$screen_id == 'MTS010',]
drug_screen_CDKS <- drug_screen_CDKS[,c('depmap_id','name','auc')]

drug_screen_CDKS <- drug_screen_CDKS %>%
  filter(!is.na(depmap_id))

drug_screen_CDKS_wide <- drug_screen_CDKS %>%
  pivot_wider(names_from = name, values_from = auc)

# merge 
NAMES <-dependency_data$ModelID
dependencia_genes_adjusted <- dependency_data[,-c(1)]
row.names(dependencia_genes_adjusted) <- NAMES

merge_dep_drug <- merge(dependencia_genes_adjusted,drug_screen_CDKS_wide[,c('depmap_id','ribociclib')], by.x = 0, by.y = 'depmap_id')
NAMES2 <- merge_dep_drug$Row.names
merge_dep_drug <- merge_dep_drug[,-c(1)]
row.names(merge_dep_drug) <- NAMES2

# correlation
palbo <- merge_dep_drug$palbociclib
palbo <- merge_dep_drug$abemaciclib
palbo <- merge_dep_drug$ribociclib


genes_df <- merge_dep_drug[, setdiff(colnames(merge_dep_drug), "palbociclib")]
genes_df <- merge_dep_drug[, setdiff(colnames(merge_dep_drug), "abemaciclib")]
genes_df <- merge_dep_drug[, setdiff(colnames(merge_dep_drug), "ribociclib")]

rho <- sapply(
  genes_df,
  function(x) cor(x, palbo, method = "spearman", use = "pairwise.complete.obs")
)

pvals <- sapply(
  genes_df,
  function(x) cor.test(x, palbo, method = "spearman", use = "pairwise.complete.obs")$p.value
)

res_df <- data.frame(
  gene = names(rho),
  spearman_rho = rho,
  p_value = pvals
)

res_df$FDR <- p.adjust(res_df$p_value, method = "fdr")
res_df <- res_df[order(res_df$spearman_rho), ]
res_acq <- res_df

FC <- 0.1

res_acq$delabel <- ifelse((res_acq$spearman_rho > FC & res_acq$FDR < 0.05), 'Up', ifelse((res_acq$spearman_rho < -FC & res_acq$FDR < 0.05), 'Down','NA'))
res_acq$gene_lab <- ifelse(( abs(res_acq$spearman_rho) > FC & res_acq$FDR < 0.05), res_acq$gene, '')


res_acq$Target <- ifelse(( abs(res_acq$spearman_rho) > FC & res_acq$FDR < 0.05),'Yes', 'No')

res_df <- as.data.frame(res_acq)
res_df$genes <- rownames(res_df)

# plotting
png('/Users/asanchezb/Desktop/test.png', width = 1200, height = 1200, res = 250)
ggplot(data=res_df, aes(x=spearman_rho, y=-log10(FDR), label=gene_lab)) + 
  geom_point(aes(shape=delabel,fill=delabel, size=delabel,color=Target)) +
  theme_minimal() +
  scale_shape_manual(values=c(21,20,21))+
  geom_text_repel(box.padding = 0.7,  # Adjust this value to control the nudge towards the left
                  segment.size = 0.2,  # Adjust the segment size to control the line segment's length
                  segment.color = "black",
                  max.overlaps = Inf,
                  size = 2.5# Adjust segment color to match the theme
  ) + # box.padding = 0.5,max.overlaps = Inf
  #scale_fill_manual(values=c("#0080FF", '#E0E0E0',"#FF6666")) +
  scale_fill_manual(values=c("#6A81A0", '#E0E0E0',"#E99A98")) +
  scale_color_manual(values=c("#E0E0E0","black"))+
  scale_size_manual(values=c(1.7,0.9,1.7))+
  geom_vline(xintercept=c(0),col="#C0C0C0", size=0.2) +
  #geom_hline(yintercept=-log10(0.05), col="black") +
  #ggtitle('Responders VS Non-responders') +
  xlim(c(-0.85,0.85)) +
  xlab("FC") + 
  ylab("-log10(padj)") +
  theme(panel.grid = element_blank(),
        plot.margin = margin(30, 30, 30, 30, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))+
  guides(color = FALSE, shape = FALSE, fill=FALSE, size=FALSE) 
dev.off()


#############################################################################################################
# 5. CORRELATION EXPRESSION AND PRISM
#############################################################################################################

# use MTS010

table(drug_screen_CDKS$screen_id) # HTS002 MTS010
drug_screen_CDKS <-drug_screen[drug_screen$name %in% CDK_inh,]

drug_screen_CDKS <- drug_screen_CDKS[drug_screen_CDKS$screen_id == 'MTS010',]
drug_screen_CDKS <- drug_screen_CDKS[,c('depmap_id','name','auc')]

drug_screen_CDKS <- drug_screen_CDKS %>%
  filter(!is.na(depmap_id))

drug_screen_CDKS_wide <- drug_screen_CDKS %>%
  pivot_wider(names_from = name, values_from = auc)

# merge 

merge_dep_drug <- merge(expression_data,drug_screen_CDKS_wide[,c('depmap_id','ribociclib')], by.x = 0, by.y = 'depmap_id')
NAMES2 <- merge_dep_drug$Row.names
merge_dep_drug <- merge_dep_drug[,-c(1)]
row.names(merge_dep_drug) <- NAMES2

# correlation
palbo <- merge_dep_drug$palbociclib
palbo <- merge_dep_drug$abemaciclib
palbo <- merge_dep_drug$ribociclib


genes_df <- merge_dep_drug[, setdiff(colnames(merge_dep_drug), "palbociclib")]
genes_df <- merge_dep_drug[, setdiff(colnames(merge_dep_drug), "abemaciclib")]
genes_df <- merge_dep_drug[, setdiff(colnames(merge_dep_drug), "ribociclib")]

rho <- sapply(
  genes_df,
  function(x) cor(x, palbo, method = "spearman", use = "pairwise.complete.obs")
)

pvals <- sapply(
  genes_df,
  function(x) cor.test(x, palbo, method = "spearman", use = "pairwise.complete.obs")$p.value
)

res_df <- data.frame(
  gene = names(rho),
  spearman_rho = rho,
  p_value = pvals
)

res_df$FDR <- p.adjust(res_df$p_value, method = "fdr")
res_df <- res_df[order(res_df$spearman_rho), ]
res_acq <- res_df

list_genes_int <- c(
  "RB1", "CDKN2D",
  "CDKN2A", "CDKN2B", "CCNE1", "CCNE2", "CCNG1",
  "CCNP", "ESD", "TACSTD2"
)

FC <- 0.1
res_acq$spearman_rho <- res_acq$spearman_rho*-1

res_acq$delabel <- ifelse((res_acq$spearman_rho > FC & res_acq$FDR < 0.05), 'Up', ifelse((res_acq$spearman_rho < -FC & res_acq$FDR < 0.05), 'Down','NA'))
res_acq$gene_lab <- ifelse(( abs(res_acq$spearman_rho) > FC & res_acq$FDR < 0.05), res_acq$gene, '')
res_acq$gene_lab <- ifelse(res_acq$gen %in% list_genes_int, res_acq$gen, '')


res_acq$Target <- ifelse(( abs(res_acq$spearman_rho) > FC & res_acq$FDR < 0.05),'Yes', 'No')

res_df <- as.data.frame(res_acq)
res_df$genes <- rownames(res_df)

# plotting
png('/Users/asanchezb/Desktop/test.png', width = 1200, height = 1200, res = 250)
ggplot(data=res_df, aes(x=spearman_rho, y=-log10(FDR), label=gene_lab)) + 
  geom_point(aes(shape=delabel,fill=delabel, size=delabel,color=Target)) +
  theme_minimal() +
  scale_shape_manual(values=c(21,20,21))+
  geom_text_repel(box.padding = 0.7,  # Adjust this value to control the nudge towards the left
                  segment.size = 0.2,  # Adjust the segment size to control the line segment's length
                  segment.color = "black",
                  max.overlaps = Inf,
                  size = 2.5# Adjust segment color to match the theme
  ) + # box.padding = 0.5,max.overlaps = Inf
  #scale_fill_manual(values=c("#0080FF", '#E0E0E0',"#FF6666")) +
  scale_fill_manual(values=c("#6A81A0", '#E0E0E0',"#E99A98")) +
  scale_color_manual(values=c("#E0E0E0","black"))+
  scale_size_manual(values=c(1.7,0.9,1.7))+
  geom_vline(xintercept=c(0),col="#C0C0C0", size=0.2) +
  #geom_hline(yintercept=-log10(0.05), col="black") +
  #ggtitle('Responders VS Non-responders') +
  xlim(c(-0.35,0.35)) +
  xlab("FC") + 
  ylab("-log10(padj)") +
  theme(panel.grid = element_blank(),
        plot.margin = margin(30, 30, 30, 30, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))+
  guides(color = FALSE, shape = FALSE, fill=FALSE, size=FALSE) 
dev.off()


#############################################################################################################
# 6. AUC VALUES IN CLUSTERS
#############################################################################################################


#use MTS010
# read cluster assignment
cluster <- read_excel('/Users/asanchezb/Desktop/PROJECTS/CODEPENDENCIES_CDKS_CYCLINS/PAPER/Supplementary_tables/TableS3_Cluster_assigment_cell_lines.xlsx',sheet = 1)
table(cluster$cluster)

merge_dep_drug <- merge(cluster,drug_screen_CDKS_wide[,c('depmap_id','ribociclib')], by.x = 'ModelID', by.y = 'depmap_id')
NAMES2 <- merge_dep_drug$Row.names
merge_dep_drug <- merge_dep_drug[,-c(1)]
row.names(merge_dep_drug) <- NAMES2



merge_dep_drug$cluster <- factor(
  merge_dep_drug$cluster,
  levels = c(1, 2, 3, 4, 5, 6)
)
cluster_colors <- c("1" = "#CD3333", "2" = "#1874CD", "3" = "#CD9B1D", "4" = "#CD6090", "5" = "#EE9572", "6" = "#7A378B")


table(merge_dep_drug$cluster)
png('/Users/asanchezb/Desktop/test.png', width = 800, height = 500, res = 200)
ggplot(met_filt_drug, aes(x = cluster, y = ribociclib, fill = cluster)) +
  #geom_violin(trim = FALSE, width = 1.5, alpha = 0.8) +
  geom_boxplot(width = 0.7, outlier.shape = NA, color = "black") +
  geom_jitter(
    width = 0.1,
    size = 0.6,
    alpha = 0.6,
    color = "black"
  ) +
  xlab("Cluster") +
  ylab("Ribociclib (AUC)") +
  ylim(c(0.5,1))+
  scale_fill_manual(values = cluster_colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 11, colour = "black"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  guides(fill = "none")




dev.off()



#############################################################################################################
# 7. CORRELATION MUTATIONS AND PRISM
#############################################################################################################

df_wide <- mutaciones_relevants %>%
  distinct(ModelID, HugoSymbol) %>%
  mutate(mutated = 1) %>%
  pivot_wider(
    names_from = HugoSymbol,
    values_from = mutated,
    values_fill = 0
  )

# use MTS010
table(drug_screen_CDKS$screen_id) # HTS002 MTS010
drug_screen_CDKS <-drug_screen[drug_screen$name %in% CDK_inh,]

drug_screen_CDKS <- drug_screen_CDKS[drug_screen_CDKS$screen_id == 'MTS010',]
drug_screen_CDKS <- drug_screen_CDKS[,c('depmap_id','name','auc')]

drug_screen_CDKS <- drug_screen_CDKS %>%
  filter(!is.na(depmap_id))

drug_screen_CDKS_wide <- drug_screen_CDKS %>%
  pivot_wider(names_from = name, values_from = auc)

# merge 

merge_dep_drug <- merge(df_wide,drug_screen_CDKS_wide[,c('depmap_id','ribociclib')], by.x = 'ModelID', by.y = 'depmap_id')
NAMES2 <- merge_dep_drug$Row.names
merge_dep_drug <- merge_dep_drug[,-c(1)]
row.names(merge_dep_drug) <- NAMES2

# correlation
palbo <- merge_dep_drug$palbociclib
palbo <- merge_dep_drug$abemaciclib
palbo <- merge_dep_drug$ribociclib


genes_df <- merge_dep_drug[, setdiff(colnames(merge_dep_drug), "palbociclib")]
genes_df <- merge_dep_drug[, setdiff(colnames(merge_dep_drug), "abemaciclib")]
genes_df <- merge_dep_drug[, setdiff(colnames(merge_dep_drug), "ribociclib")]

rho <- sapply(
  genes_df,
  function(x) cor(x, palbo, method = "spearman", use = "pairwise.complete.obs")
)

pvals <- sapply(
  genes_df,
  function(x) cor.test(x, palbo, method = "spearman", use = "pairwise.complete.obs")$p.value
)

res_df <- data.frame(
  gene = names(rho),
  spearman_rho = rho,
  p_value = pvals
)

res_df$FDR <- p.adjust(res_df$p_value, method = "fdr")
res_df <- res_df[order(res_df$spearman_rho), ]
res_acq <- res_df

list_genes_int <- c(
  "RB1", "CDKN2D",
  "CDKN2A", "CDKN2B", "CCNE1", "CCNE2", "CCNG1",
  "CCNP", "ESD", "TACSTD2"
)

FC <- 0.1
res_acq$spearman_rho <- res_acq$spearman_rho*-1

res_acq$delabel <- ifelse((res_acq$spearman_rho > FC & res_acq$FDR < 0.05), 'Up', ifelse((res_acq$spearman_rho < -FC & res_acq$FDR < 0.05), 'Down','NA'))
res_acq$gene_lab <- ifelse(( abs(res_acq$spearman_rho) > FC & res_acq$FDR < 0.05), res_acq$gene, '')
res_acq$gene_lab <- ifelse(res_acq$gen %in% list_genes_int, res_acq$gen, '')


res_acq$Target <- ifelse(( abs(res_acq$spearman_rho) > FC & res_acq$FDR < 0.05),'Yes', 'No')

res_df <- as.data.frame(res_acq)
res_df$genes <- rownames(res_df)

# plotting
png('/Users/asanchezb/Desktop/test.png', width = 1200, height = 1200, res = 250)
ggplot(data=res_df, aes(x=spearman_rho, y=-log10(FDR), label=gene_lab)) + 
  geom_point(aes(shape=delabel,fill=delabel, size=delabel,color=Target)) +
  theme_minimal() +
  scale_shape_manual(values=c(21,20,21))+
  geom_text_repel(box.padding = 0.7,  # Adjust this value to control the nudge towards the left
                  segment.size = 0.2,  # Adjust the segment size to control the line segment's length
                  segment.color = "black",
                  max.overlaps = Inf,
                  size = 2.5# Adjust segment color to match the theme
  ) + # box.padding = 0.5,max.overlaps = Inf
  #scale_fill_manual(values=c("#0080FF", '#E0E0E0',"#FF6666")) +
  scale_fill_manual(values=c("#6A81A0", '#E0E0E0',"#E99A98")) +
  scale_color_manual(values=c("#E0E0E0","black"))+
  scale_size_manual(values=c(1.7,0.9,1.7))+
  geom_vline(xintercept=c(0),col="#C0C0C0", size=0.2) +
  #geom_hline(yintercept=-log10(0.05), col="black") +
  #ggtitle('Responders VS Non-responders') +
  xlim(c(-0.35,0.35)) +
  xlab("FC") + 
  ylab("-log10(padj)") +
  theme(panel.grid = element_blank(),
        plot.margin = margin(30, 30, 30, 30, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))+
  guides(color = FALSE, shape = FALSE, fill=FALSE, size=FALSE) 
dev.off()


















