# Manuscript:   A reference map of CDK-cyclin dependencies across human cancers
# Project:      Codependencies
# Part:         C. Efficacy and selectivity
# Author:       Agustín Sánchez-Belmonte
# Date:         04/09/23
# Last update:  08/02/26

# Methods: Kenichi Shimada John A Bachman Jeremy L Muhlich Timothy J Mitchison (2021) shinyDepMap, a tool to identify targetable cancer genes and their functional connections from Cancer Dependency Map data eLife 10:e57116.
# github: https://github.com/kenichi-shimada/depmap-analysis/tree/various_mix

# Using the same methods but with updated data (25Q3) and only CRISPR data



# ---------------- INDEX ---------------- 


# 0. SET UP
# 1. LOAD DATA
# 2. FILTER BY TUMOR TYPE
# 3. THRESHOLD FOR EFFICACY
# 4. EFFICACY 
# 5. SELECTIVITY 
# 6. CORRELATIONS Q5TH AND Q95TH
# 7. CORRELATIONS SELECTIVITY AND EFFICACY
# 8. PLOT INDIVIDUAL DEPENDENCIES
# 9. COMPUTE SPEARMAN CORRELATION COEFFICIENTS
# 10. REDUCTIONAL DIMENSION
# 11. PLOTTING



#############################################################################################################
# 0. SET UP
#############################################################################################################


# ---------------- PATHS ---------------- 

setwd("CODEPENDENCIES_CDKS_CYCLINS/EFFICACY_AND_SELECTIVITY/")
PLOTS <- 'CODEPENDENCIES_CDKS_CYCLINS/MEDIAN_SELECTIVITY/'
#dir.create(PLOTS, showWarnings = TRUE)

# ---------------- PACKAGES ---------------- 


library(L1pack)
library(org.Hs.eg.db)
library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(readxl)
library(svglite)



# ---------------- FUNCTIONS ---------------- 

#Functions
find.threshold <- function(dep.scores,#breaks=100,
                           side.fit=c("upper","lower"),side.thres=c("upper","lower"),
                           ps=1e-3,plot=TRUE,...){
  
  # xh <- hist(dep.scores,col="lightblue",border=NA,breaks=seq.br,plot=TRUE,...)
  d <- density(dep.scores)
  mod <- d$x[which.max(d$y)]
  
  if(side.fit[1]=="upper"){
    val <- (dep.scores-mod)
    val <- val[val > 0] 
    sd <- sqrt(sum(val^2)/(length(val)-1))
  }else{
    val <- (dep.scores-mod)
    val <- val[val < 0] 
    sd <- sqrt(sum(val^2)/(length(val)-1))
  }
  
  if(!is.na(ps)){
    if(side.thres[1]=="upper"){
      ths <- qnorm(1-ps,mean=mod,sd=sd)
    }else{
      ths <- mod * 2 - qnorm(1-ps,mean=mod,sd=sd)
    }
    names(ths) <- ps
  }
  if(plot){
    plot(d$x,d$y,type="n",...)
    polygon(d$x,d$y,col="lightblue",border=NA)
    ybr <- dnorm(d$x,mean=mod,sd=sd)
    f <- max(d$y)/max(ybr)
    lines(d$x,f*ybr,col=4,lwd=2)
    if(!is.na(ps))abline(v=ths,col=2)
  }
  
  return(ths)
}







#############################################################################################################
# 1. LOAD DATA
#############################################################################################################


# ---------------- DEPENDENCY ---------------- 

dependency_data <- read.csv('25Q3_data/CRISPRGeneEffect.csv', check.names = F)

names <- c()
for(i in colnames(dependency_data)){
  new <- strsplit(i, " ")[[1]][1]
  names <- c(names,new)
}

colnames(dependency_data) <- names
names(dependency_data)[1] <- 'ModelID'

cell_names <- dependency_data$ModelID
dependency_data_new <- dependency_data[,-c(1)]
row.names(dependency_data_new) <- cell_names

# ---------------- METADATA ---------------- 

info <- read.csv('25Q3_data/Model.csv', check.names = F)
info <- read_excel('Supplementary_tables/TableS3_Cluster_assigment_cell_lines.xlsx', sheet = 1)

# ---------------- MERGE ---------------- 

info <- info[info$ModelID %in% rownames(dependency_data_new),]
dependency_data_new<- dependency_data_new[rownames(dependency_data_new) %in% info$ModelID,]
dependency_data_new<- dependency_data_new[info$ModelID,]



#############################################################################################################
# 2. FILTER BY TUMOR TYPE
#############################################################################################################

# Choose tumor type: Breast or Lung
tumor_type = 'Breast'
cluster = 'A'

filt_cell_lines <- info[info$OncotreeLineage == tumor_type,]$ModelID # TUMOR TYPE
filt_cell_lines <- info[info$cluster == cluster,]$ModelID # CLUSTER
dependency_data_sub<- dependency_data_new[rownames(dependency_data_new) %in% filt_cell_lines,]

IC50_data <- IC50_data[IC50_data$Lineage ==tumor_type,]

# fill NA values in Lilly data with Finn et al. paper (only BRCA in palbo)
IC50_data$Sensitivity_Lilly_Palbo <- ifelse(is.na(IC50_data$Sensitivity_Lilly_Palbo), IC50_data$Sensitivity_finn_palbo, IC50_data$Sensitivity_Lilly_Palbo)
sum(info$StrippedCellLineName %in% IC50_data$Cell_line)

# merge fata with info
merge_data <- merge(IC50_data, info,by.x='Cell_line', by.y='StrippedCellLineName')
sum(merge_data$ModelID %in% rownames(dependency_data_sub))

# Separate in resistant and sensitive
# Sensitivity_Lilly_Palbo, Sensitivity_Lilly_Ribo, Sensitivity_Lilly_Abema
Resistant_cell_lines <- merge_data[merge_data$Sensitivity_Lilly_Palbo == 'Resistant',]$ModelID
Sensitive_cell_lines <- merge_data[merge_data$Sensitivity_Lilly_Palbo == 'Sensitive',]$ModelID

resistant_dep <- dependency_data_sub[rownames(dependency_data_sub) %in% Resistant_cell_lines,]
sensitive_dep <- dependency_data_sub[rownames(dependency_data_sub) %in% Sensitive_cell_lines,] 



#############################################################################################################
# 3. THRESHOLD FOR EFFICACY
#############################################################################################################

## cutoff values (quantiles) for efficacy and selectivity
qtl <- c(0.5,1,2.5,5,10,25)
qtl2 <- 100-qtl
#cr.ratio <- paste0("r",seq(0,1,.2))

## find essential and tumor suppressor genes
sc.ths <- is.ts <- is.ess <- list()
n.ts.gns <- n.ess.gns <- array(NA,c(1,length(qtl)))

is.ts[[1]] <- is.ess[[1]] <- list()

ef <- t(resistant_dep)
ef <- t(dependency_data)
ef <- na.omit(ef) # added for new version
  
## compute thresholds for essential and tumor suppressor genes, based on data distribution (To)
sc <- as.vector(ef)
sc <- sc[!is.na(sc)]
  
p.th <- 1e-3

pdf(paste0(PLOTS,"eff_th.svg"), width = 5, height = 5)
sc.ths[[1]] <- find.threshold(sc,
                              xlab="dependency score",plot=TRUE,
                              main="", ps=c(p.th,1-p.th),
                              side.fit="upper",side.thres="lower")
dev.off() 


k = 1
# Lower ratio 
pdf(paste0(PLOTS,"lower_ratio.pdf"),width=5,height=5)
par(mar=c(5,5,3,3))
for(i in 1:6){
    ## efficacy -qtl[i]-th quantile score
    qth <- apply(ef,1,function(x)quantile(x,qtl[i]/100,na.rm=T))
    
    ##
    ths <- find.threshold(qth,
                          xlab="dependency score",plot=TRUE,
                          main=paste0(qtl[i]," percentile"), ps=0,
                          side.fit="upper",side.thres="lower",
                          xlim=c(-3,2))
    abline(v=sc.ths[[1]],col=2)
    
    is.ess[[1]][[i]] <- idx <- qth < sc.ths[[1]][1]
    n.ess.gns[1,i] <- n.genes <- sum(idx)
    
}
dev.off()

# Upper ratio 
pdf(paste0(PLOTS,"upper_ratio.pdf"),width=5,height=5)
par(mar=c(5,5,3,3))

for(i in 1:6){
  qth2 <- apply(ef,1,function(x)quantile(x,qtl2[i]/100,na.rm=T))
  ths <- find.threshold(qth2,
                        xlab="dependency score",plot=TRUE,
                        main=paste(qtl2[i]," percentile"),ps=0,
                        side.fit="upper",side.thres="upper",
                        xlim=c(-3,2))
  abline(v=sc.ths[[k]],col=2)
  
  is.ts[[k]][[i]] <- idx <- qth2 > sc.ths[[k]][2]
  n.ts.gns[k,i] <- n.genes <- sum(idx)
  
}
dev.off()

  
# Select essential genes
names(sc.ths) <- names(is.ess) <- names(is.ts) <- 'CRISPR'
is.ess <- lapply(is.ess,function(x){
  names(x) <- paste0("q",qtl)
  return(x)
})
is.ts <- lapply(is.ts,function(x){
  names(x) <- paste0("q",qtl)
  return(x)
})
rownames(n.ess.gns) <- rownames(n.ts.gns) <- 'CRISPR'
colnames(n.ess.gns) <- colnames(n.ts.gns) <- paste0("q",qtl)

is.non <- lapply(colnames(n.ess.gns),function(r)!is.ess$CRISPR[[r]] & !is.ts$CRISPR[[r]])
names(is.non) <- colnames(n.ess.gns)

# Data fraem with genes of each group in q1 and CRISPR (low, high and non), here you choose the qth (before:q1, now:q5)
lo.tab <- sapply(is.ess,function(x)x[[4]])
hi.tab <- sapply(is.ts,function(x)x[[4]])
non.tab <- as.data.frame(is.non$q1)
colnames(non.tab) <- 'CRISPR'

# Amount og gene en each group depend of quartiles
sapply(is.ess$CRISPR,sum)
sapply(is.ts$CRISPR,sum)
sapply(is.non,sum)

# Saving objects
save(is.non,is.ess,is.ts,sc.ths,n.ess.gns,n.ts.gns,lo.tab,hi.tab,non.tab,
     file="../25Q3_RESULTS/find_essential_genes.rda")



#############################################################################################################
# 4. EFFICACY
#############################################################################################################

# Depend on the quartile selected the amount of essential genes change
quartile_selected = 0.25 #0.25 #selected 0.25 due to the low number of cell lines by group/pancancer 0.05/ subgroup 0.25
quartile_selected2 = 1-quartile_selected

# Calculate quartiles
qth <- apply(ef,1,function(x)quantile(x,quartile_selected,na.rm=T))
qth2 <- apply(ef,1,function(x)quantile(x,quartile_selected2,na.rm=T))

# regression model
lf1 <- l1fit(qth,qth2)

# Plot correlation between 5th quartile and 95th quartile
pdf(paste0(PLOTS,"qth_qth2_correlation.pdf"),width=8,height=8)
smoothScatter(qth,qth2,pch=20,
              main='Quartile correlation')
abline(lf1,col=2)
abline(h=0,v=0)
dev.off()

# Calculate dispersion
rth <- qth2-qth
# regression model
fit  <- l1fit(qth,rth)
#plot
pdf(paste0(PLOTS,"fitted_rth_correlation.pdf"),width=8,height=8)
smoothScatter(fit$fitted.values,rth,pch=20,main='Fitted values and dispersion correlation')
abline(a=0,b=1)
abline(h=0,v=0)
dev.off()





#############################################################################################################
# 5. SELECTIVITY
#############################################################################################################

# Calculate selesctivity dispersion/fitted values-1
sel <- rth/fit$fitted.values-1
ef.sel <- data.frame(eff=qth,sel=sel)
ef.sel$Gene <- rownames(ef.sel)
ef.sel$Score <- ef.sel$sel*(-ef.sel$eff)
Lung_sen.ef.sel <- ef.sel
write.table(ef.sel,paste0(PLOTS,"Pancancer_sen.ef.sel.tsv"), sep = '\t', row.names = F)

min.qth <- min(qth)
#x <- hist(qth,breaks=100)
mid.qth <- x$mids[which.max(x$density)]


syms <- c('CDK1','CDK2','CDK3','CDK4','CDK5','CDK6','CDK7','CDK8','CDK9','CDK10','CDK12','CDK13','CDK14','CDK15','CDK16','CDK17','CDK18','CDK19','CDK20','CCND1','CCND2','CCND3','CCNC','CCNE1','CCNE2','CCNF','CCNA1','CCNA2','CCNB1','CCNB2','CCNB3','CCNG1','CCNG2','CCNH','CCNI','CCNI2','CCNJ','CCNJL','CCNK','CCNL1','CCNL2','CCNO','CCNQ','CCNP','CCNT1','CCNT2','CCNY','CCNYL1','CNPPD1','CNTD1','SPDYA','SPDYE1','SPDYE3','SPDYE4','CABLES1','CABLES2','CDK5R1','CDK5R2') #,'CDKL1','CDKL2','CDKL3','CDKL4','CDKL5','KIF18A','KIF18B','ERBB2'
# Ordenar por score descendente
# Seleccionar top 10
top10_hits <- head(ef.sel[order(-ef.sel$Score),'Gene' ], 10)
qth.1 <- qth[syms]
qth2.1 <- qth2[syms]


#############################################################################################################
# 6. CORRELATIONS Q5TH AND Q95TH
#############################################################################################################

# load data
#load("../24Q2_RESULTS//eff_sel.rda")
#save(sel.th,ef.sel,ef.sel.1,
     #file="eff_and_sel/eff_sel.rda")

# A. CDKs

non_sel <- c('CDK1','CDK7','CDK9','CDK12')
selective <- c('CDK2','CDK4','CDK6','CDK8','CDK13')
others <- c('CDK3','CDK5','CDK10','CDK14','CDK15','CDK16','CDK17','CDK18','CDK19','CDK20')
selected <- c(2,4,6,8)
selected2 <- c(1,9,11)

pdf(paste0(PLOTS,"qth_qth2_CDKs.pdf"),width=8,height=8)
# png("qth_qth2_v2.png",width=800,height=800)
par(mar=c(5,5,3,3))
cols <- densCols(data.frame(qth,qth2),colramp = colorRampPalette(blues9[-(1:3)]))
plot(qth,qth2,pch=20,col=cols,type="n",
     xlab="5th percentile score",ylab="95th percentile score")
#cex.lab=2,cex.axis=2)
abline(h=0,v=0,lty=1,col="grey80")
abline(a=0,b=1,lty=2)
points(qth,qth2,pch=20,col=cols,cex=.5)
abline(lf1$coefficients,col=2)
cols <- brewer.pal(11,"Spectral")[c(1,4,7)]
points(qth.1[non_sel],qth2.1[non_sel],bg=cols[1],col=1,cex=1,pch=21)
points(qth.1[selective],qth2.1[selective],bg=cols[2],col=1,cex=1,pch=21)
points(qth.1[others],qth2.1[others],bg=cols[3],col=1,cex=1,pch=21)
text(qth.1[selected],qth2.1[selected]+0.14,syms[selected])
text(qth.1[selected2],qth2.1[selected2]-0.14,syms[selected2])
text(qth.1[7],qth2.1[7]-0.17,syms[7])
text(qth.1[12],qth2.1[12]-0.14,syms[12])
dev.off()



# B. CYCLINS

non_sel <- c('CCNK','CCNA2','CCNH')
selective <- c('CCND1','CCND2','CCND3','CCNE1','CCNE2','CCNF','CCNL1','CCNB1','SPDYE1','CCNC')
others <- c('CCNA1','CCNB2','CCNB3','CCNG1','CCNG2','CCNI','CCNI2','CCNJ','CCNJL','CCNL2','CCNO','CCNT1','CCNT2','CCNY','CCNYL1','SPDYA','SPDYE1','SPDYE3','SPDYE4','CABLES1','CABLES2','CDK5R1','CDK5R2','CDKL1','CDKL2','CDKL3','CDKL4','CDKL5')

pos_non_sel <- c(28,34,39)
pos_sel <- c(20,21,23,25)
pos_othres <- c(27,30,31,32,33,35,36,37,38,41,42,43,44,45,46,47,49,50,51,52,53,54,55,56,57,58,59)

pdf(paste0(PLOTS,"qth_qth2_CCNs.pdf"),width=8,height=8)
# png("qth_qth2_v2.png",width=800,height=800)
par(mar=c(5,5,3,3))
cols <- densCols(data.frame(qth,qth2),colramp = colorRampPalette(blues9[-(1:3)]))
plot(qth,qth2,pch=20,col=cols,type="n",
     xlab="5th percentile score",ylab="95th percentile score")
#cex.lab=2,cex.axis=2)
abline(h=0,v=0,lty=1,col="grey80")
abline(a=0,b=1,lty=2)
points(qth,qth2,pch=20,col=cols,cex=.5)
abline(lf1$coefficients,col=2)
cols <- brewer.pal(11,"Spectral")[c(1,4,7)]
points(qth.1[non_sel],qth2.1[non_sel],bg=cols[1],col=1,cex=1,pch=21)
points(qth.1[selective],qth2.1[selective],bg=cols[2],col=1,cex=1,pch=21)
points(qth.1[others],qth2.1[others],bg=cols[3],col=1,cex=1,pch=21)
text(qth.1[pos_non_sel],qth2.1[pos_non_sel]-0.14,syms[pos_non_sel])
text(qth.1[pos_sel],qth2.1[pos_sel]+0.12,syms[pos_sel])
text(qth.1[24]+0.34,qth2.1[24]+0.07,syms[24])
text(qth.1[48]+0.42,qth2.1[48]-0.05,syms[48])
text(qth.1[29]+0.30,qth2.1[29]+0.07,syms[29])
text(qth.1[22]-0.26,qth2.1[22]-0.11,syms[22])
text(qth.1[40]-0.33,qth2.1[40]+0.03,syms[40])
text(qth.1[26]+0.18,qth2.1[26]+0.11,syms[26])
dev.off()


#############################################################################################################
# 7. CORRELATIONS SELECTIVITY VS EFFICACY
#############################################################################################################
eff <- ef.sel$eff
sel <- ef.sel$sel

ef.sel.1 <- ef.sel[syms,]

eff.1 <- ef.sel.1$eff
sel.1 <- ef.sel.1$sel

ef.sel.2 <- ef.sel[top10_hits,]

eff.2 <- ef.sel.2$eff
sel.2 <- ef.sel.2$sel

#sel.th<- 0.11

# Histograms for efficacy and selectivity
xhist <- hist(eff,breaks=1000,col="pink")
yhist <- hist(sel,breaks=1000,col="pink")
top <- max(xhist$counts,yhist$counts)

top <- max(yhist$counts)


#gids.all <- rownames(sc)

# Selectivity threshold
sel.th <- find.threshold(sel,
                        xlab="selectivity score",plot=FALSE,
                        main="", ps=1e-3,
                        side.fit="lower",side.thres="upper")

effi <- which(ef.sel.1$eff<sc.ths[[1]][1] & ef.sel.1$sel<sel.th )
selec <- which(ef.sel.1$sel>sel.th)
others <- which(ef.sel.1$eff>sc.ths[[1]][1] & ef.sel.1$sel<sel.th )

pdf(paste0(PLOTS,"eff_sel_limit_without_points.pdf"),width=8,height=8)
# Plot base vacío
plot(eff, sel,
     pch = 20,
     type = "n",
     xlab = "Efficacy",
     ylab = "Selectivity",
     xlim = c(-3.5,0.5),
     ylim = c(-0.5,2.6))

abline(h=0,v=0,lty=1,col="grey80")
abline(v=sc.ths[[1]][1],col="grey20",lty=2)
abline(h=sel.th,col="grey20",lty=2)

# 🔹 Todos los puntos en gris transparente
points(eff, sel,
       pch = 20,
       col = rgb(0.6,0.6,0.6,0.4),  # gris con alpha
       cex = 0.5)

# 🔹 Labels solo para seleccionados
text(eff.1 - 0.25,
     sel.1,
     syms,
     col = "blue")

# 🔹 Labels solo para seleccionados
text(eff.2 - 0.25,
     sel.2,
     top10_hits,
     col = "blue")

# 🔹 Puntos seleccionados en azul
points(eff.1,sel.1,bg='#87CEFA',col=1,cex=1,pch=21)
points(eff.2,sel.2,bg='#FFA07A',col=1,cex=1,pch=21)

dev.off()

list_dfs <- list(
  A = A.ef.sel,
  B = B.ef.sel,
  C = C.ef.sel,
  D = D.ef.sel,
  E = E.ef.sel,
  F1 = F1.ef.sel,
  F2 = F2.ef.sel,
  F3 = F3.ef.sel,
  F4 = F4.ef.sel,
  F5 = F5.ef.sel,
  F6 = F6.ef.sel
)

list_dfs <- list(
  BS = Breast_sen.ef.sel,
  BR = Breast_res.ef.sel,
  LS = Lung_sen.ef.sel,
  LR = Lung_res.ef.sel
)

renamed_list <- lapply(names(list_dfs), function(n) {
  df <- list_dfs[[n]]
  colnames(df)[colnames(df) != "Gene"] <- 
    paste0(colnames(df)[colnames(df) != "Gene"], "_", n)
  df
})

merged_eff_2 <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE),
                 renamed_list)
write.table(merged_eff_2, '25Q3_RESULTS/merged_eff_lung_breast.tsv', sep = '\t')

#############################################################################################################
# 8. PLOT INDIVIDUAL DEPENDENCIES
#############################################################################################################


dependency_data_sub$Group <- ifelse(rownames(dependency_data_sub) %in% Resistant_cell_lines, 'Resistant','')
dependency_data_sub$Group <- ifelse(rownames(dependency_data_sub) %in% Sensitive_cell_lines, 'Sensitive',dependency_data_sub$Group)

df_plot <- dependency_data_sub %>%
  filter(Group %in% c("Sensitive", "Resistant")) %>%
  mutate(Group = factor(Group, levels = c("Sensitive", "Resistant")))

cols <- c(
  "Sensitive" = "#E99A98",
  "Resistant" = "#6A81A0"
)


genes <- c("CCND1","CCNE1","CCNE2",
           "CDK2","CDK4","CDK6",
           "E2F1","E2F3","KIF18A")

df_long <- df_plot %>%
  select(Group, all_of(genes)) %>%
  pivot_longer(cols = all_of(genes),
               names_to = "Gene",
               values_to = "Dependency")
df_long$Gene_num <- as.numeric(factor(df_long$Gene))




png('/Users/asanchezb/Desktop/test.png',
    width = 1400, height = 800, res = 250)

ggplot(df_long,
       aes(x = Gene,
           y = Dependency,
           fill = Group)) +
  
  geom_boxplot(position = position_dodge(width = 0.7),
               width = 0.6,
               outlier.shape = NA,
               alpha = 0.8) +
  
  geom_jitter(aes(color = Group),
              position = position_jitterdodge(
                jitter.width = 0.15,
                dodge.width = 0.7),
              size = 0.8,
              alpha = 0.6) +
 
  
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  
  xlab("") +
  ylab("Dependency score")

dev.off()


# Medianas por grupo
# Calcular tabla resumen
results_table <- df_long %>%
  group_by(Gene) %>%
  summarise(
    n_Sensitive = sum(Group == "Sensitive"),
    n_Resistant = sum(Group == "Resistant"),
    
    median_Sensitive = median(Dependency[Group == "Sensitive"], na.rm = TRUE),
    median_Resistant = median(Dependency[Group == "Resistant"], na.rm = TRUE),
    
    p_value = wilcox.test(Dependency ~ Group)$p.value
  ) %>%
  arrange(p_value)

results_table



#############################################################################################################
# 9. COMPUTE SPEARMAN CORRELATION COEFFICIENTS
#############################################################################################################

# This part has to be done in cluster

x <- load(file="eff_and_sel/find_essential_genes.rda") #is.non,is.ess,is.ts,sc.ths,n.ess.gns,n.ts.gns

gns <- rownames(ef)

ess.genes <- gns[is.ess$CRISPR$q1] #select Quartile Q5, Q25

system.time(
  coef <- cor(t(ef[ess.genes,]),method="spearman",use="pairwise.complete.obs")
  # coef <- cor(t(ef),method="spearman",use="pairwise.complete.obs")
)

#saveRDS(coef,file="eff_and_sel/all_coef.rds")


#############################################################################################################
# 10. DIMENSIONAL REDUCTION
#############################################################################################################

# This part has to be done in cluster
library(Rtsne) # 0.15
library(umap)
library(ggplot2)

# Load data
coef <- readRDS(file="../25Q3_RESULTS/ess_coef.rds")
dim(coef)

d <- 1-abs(coef)


# remove mitochondrial genes and genes with NA values
#rem <- read_excel('/Users/asanchezb/Desktop/remove_genes.xlsx', sheet = 1)
#d <- d[!rownames(d) %in% rem$NA_VALUES,!colnames(d) %in% rem$NA_VALUES]
#d <- d[!rownames(d) %in% rem$MIT_GENES,!colnames(d) %in% rem$MIT_GENES]

# set seed
set.seed(123)

# tsne
rt.outs <- Rtsne(X=d,
                 is_distance=TRUE,
                 perplexity=3,
                 max_iter=3000, #20000
                 pca=FALSE,
                 theta=0.8) #0.5

rt <- rt.outs$Y
rownames(rt) <- rownames(d)

# umap
umap_result <- umap(
  d,
  n_neighbors = 15,
  n_components = 2,  
  init = "pca"
)
umap_coords <- umap_result$layout
rownames(umap_coords) <- rownames(d)

saveRDS(rt,'tsne.rds')
saveRDS(umap_coords, 'umap.rds')



#############################################################################################################
# 11. PLOTTING
#############################################################################################################

library(plotly) 

# From here we can do in local
syms <- c('CDK1','CDK2','CDK3','CDK4','CDK5','CDK6','CDK7','CDK8','CDK9','CDK10','CDK12','CDK13','CDK14','CDK15','CDK16','CDK17','CDK18','CDK19','CDK20','CCND1','CCND2','CCND3','CCNC','CCNE1','CCNE2','CCNF','CCNA1','CCNA2','CCNB1','CCNB2','CCNB3','CCNG1','CCNG2','CCNH','CCNI','CCNI2','CCNJ','CCNJL','CCNK','CCNL1','CCNL2','CCNO','CCNQ','CCNP','CCNT1','CCNT2','CCNY','CCNYL1','CNPPD1','CNTD1','SPDYA','SPDYE1','SPDYE3','SPDYE4','CABLES1','CABLES2','CDK5R1','CDK5R2') 


# Read results from cluster
umap_coords <- readRDS('CODEPENDENCIES_CDKS_CYCLINS/25Q3_RESULTS/umap.rds')
rt <- readRDS('tsne.rds')

# ---------------- UMAP ---------------- 

umap_df <- data.frame(UMAP1 = umap_coords[, 1], UMAP2 = umap_coords[, 2])
umap_df$label <- ifelse(rownames(umap_coords) %in% syms, rownames(umap_coords),'')


png(filename = '/Users/asanchezb/Desktop/test.png', width = 700, height = 1000, res = 250)#1300#2500
ggplot(umap_df, aes(x = UMAP1, y = UMAP2,label=label)) +
  geom_point() +
  labs(x = "UMAP1", y = "UMAP2") +
  geom_text_repel(box.padding = 1,
                  max.overlaps = Inf,   # Adjust this value to control the nudge towards the left
                  segment.size = 0.2,  # Adjust the segment size to control the line segment's length
                  segment.color = "red", 
                  colour='red') +
  theme_void()# Adjust segment color to match the theme
dev.off()

# Interactive plot

umap_plot <- plot_ly(data = umap_df, x = ~UMAP1, y = ~UMAP2, type = 'scatter', mode = 'markers')

# Set the hovertext to the row names
umap_plot <- umap_plot %>% add_trace(
  text = ~rownames(umap_df),
  hoverinfo = "text+x+y"
)

# To display the plot, you can use the following code:
umap_plot


# ---------------- tSNE ---------------- 


tsne_df <- data.frame(TSNE1 = rt[, 1], TSNE2 = rt[, 2])
tsne_df$label <- ifelse(rownames(rt) %in% syms, rownames(rt),'')


tsne_df$has_label <- tsne_df$label != ""

png(filename = '/Users/asanchezb/Desktop/test.png', width = 1000, height = 1000, res = 200)

ggplot(tsne_df, aes(x = TSNE1, y = TSNE2)) +
  
  # SIN label (vacío) → pequeños
  geom_point(
    data = subset(tsne_df, label == ""),
    color = "grey70",
    size = 0.3
  ) +
  
  # Empieza por "cdk"
  geom_point(
    data = subset(tsne_df, label != "" & grepl("^CDK", label)),
    color = "#CD5555",
    size = 1
  ) +
  
  # Con label pero NO empieza por "cdk"
  geom_point(
    data = subset(tsne_df, label != "" & !grepl("^CDK", label)),
    color = "#36648B",
    size = 1
  ) +
  
  geom_text_repel(
    data = subset(tsne_df, label != ""),
    aes(label = label),
    box.padding = 1,
    max.overlaps = Inf,
    segment.size = 0.2,
    segment.color = "red",
    colour = "red"
  ) +
  
  theme_void()

dev.off()



# Interactive plot
tsne_plot <- plot_ly(data = tsne_df, x = ~TSNE1, y = ~TSNE2, type = 'scatter', mode = 'markers')

# Set the hovertext to the row names
tsne_plot <- tsne_plot %>% add_trace(
  text = ~rownames(tsne_df),
  hoverinfo = "text+x+y"
)

# To display the plot, you can use the following code:
tsne_plot



