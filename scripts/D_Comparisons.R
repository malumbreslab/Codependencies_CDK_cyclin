# Manuscript:   A reference map of CDK-cyclin dependencies across human cancers
# Project:      Codependencies
# Part:         D. Comparisons
# Author:       Agustín Sánchez-Belmonte
# Date:         04/09/23
# Last update:  08/02/26


# ---------------- INDEX ---------------- 


# 0. SET UP
# 1. LOAD DATA
# 2. PREPARE DATA
# 3. CALCULATE THRESHOLDS
# 4. NINE-SQUEARE
# 5. RANK PLOT
# 6. BAR PLOT




#############################################################################################################
# 0. SET UP
#############################################################################################################

# ---------------- PATHS ---------------- 

setwd('CODEPENDENCIES_CDKS_CYCLINS')


# ---------------- PACKAGES ---------------- 

library(readxl)
library(ggplot2)
library(ggrepel)

figure_font_size=12

#############################################################################################################
# 1. LOAD DATA
#############################################################################################################
#BREAST
sen_data <- read_excel('Supplementary_tables/TableS7_Efficacy_selectivity_PanCancer_Breast.xlsx', sheet = 1)
res_data <- read_excel('Supplementary_tables/TableS7_Efficacy_selectivity_PanCancer_Breast.xlsx', sheet = 2)

#############################################################################################################
# 2. PREPARE DATA
#############################################################################################################

# merge data
df_merged <- merge(sen_data, res_data, by = 'Gene')
names <- df_merged$Gene

# Calculate diff score
df_merged$diff <- df_merged$Efficacy_BS - df_merged$Efficacy_BR

# Calculate Rank
df_merged <- df_merged[order(df_merged$diff, decreasing = T),]
df_merged$Rank <- rev(seq(1,length(df_merged$diff)))

#############################################################################################################
# 3. CALCULATE THRESHOLDS
#############################################################################################################

d <- (df_merged$Efficacy_BR - df_merged$Efficacy_BS) / sqrt(2)

sd_d <- sd(d, na.rm = TRUE)

k <- 3   # o 1.5 / 3 según rigor
upper <-  k * sd_d
lower <- -k * sd_d

offset <- k * sd_d * sqrt(2)


df_merged <- df_merged %>%
  mutate(
    dist_diag = (Efficacy_BR - Efficacy_BS) / sqrt(2),
    class = case_when(
      dist_diag >  upper ~ "Sensitive",
      dist_diag <  lower ~ "Resistant",
      TRUE               ~ "Shared"
    )
  )


#############################################################################################################
# 4. NINE SQUARE
#############################################################################################################
df_labels <- df_merged %>%
  arrange(dist_diag) %>%
  slice_head(n = 50) %>%
  bind_rows(
    df_merged %>%
      arrange(desc(dist_diag)) %>%
      slice_head(n = 50)
  )


cols <- c(
  "Sensitive" = "#E99A98",
  "Shared"         = "#E0E0E0",
  "Resistant" = "#6A81A0"
)

png(filename = '/Users/asanchezb/Desktop/test.png', width = 1300, height = 1200, res = 250)#1300#2500
ggplot(df_merged, aes(Efficacy_BR, Efficacy_BS, color = class)) + #ggplot(df_merged, aes(Efficacy_BS, Efficacy_BR, color = class)) +
  geom_point(alpha = 0.7, size = 1.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid") +
  geom_abline(slope = 1, intercept =  offset, linetype = "dashed") +
  geom_abline(slope = 1, intercept = -offset, linetype = "dashed") +
  scale_color_manual(values = cols) +
  ## Labels solo para top/bottom
  geom_text_repel(
    data = df_labels,
    aes(label = Gene),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50"
  ) +

  theme_bw() +
  xlab("Sensitive") +
  ylab("Resistant") +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank()
  )
dev.off()

write.table(df_merged,paste0(PLOTS,"TEST_ef.sel.tsv"), sep = '\t')

#############################################################################################################
# 5. RANK PLOT
#############################################################################################################

png(filename = '/Users/asanchezb/Desktop/test.png', width = 1300, height = 1200, res = 250)#1300#2500
ggplot(data=df_merged, aes(x=Rank, y=diff, col=class)) + 
  geom_point(alpha = 0.7, size = 1.3) +
  theme_minimal() +
 # box.padding = 0.5,max.overlaps = Inf
  scale_color_manual(values = cols) +
  geom_text_repel(
    data = df_labels,
    aes(label = Gene),
    size = 2,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50"
  ) +
  geom_hline(yintercept=c(0), col="black", size = 0.3) +
  #geom_hline(yintercept=-log10(0.05), col="black") +
  theme(panel.grid = element_blank(),
        plot.margin = margin(15, 15, 15, 15, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))+
  ylim(c(-2.3,2.3))+
  guides(label = FALSE, color=FALSE) 
dev.off()



#############################################################################################################
# 6. BAR PLOT
#############################################################################################################

q <- read_excel('Supplementary_tables/TableS8_GSEA_resistant_vs_sensitive.xlsx', sheet = 1) 
q <- read_excel('Supplementary_tables/TableS6_GSEA_ ClusterD_vs_rest.xlsx', sheet = 1) 



#q <- subset(q, q$fdr < 0.1)
q <- subset(q, q$fdr < 0.25)
q <- q[order(q$nes),]
q$yax <- ifelse(q$nes > 0, -0.02, 0.02)
q$col <- ifelse(q$nes > 0, "#6A81A0","#E99A98")

q <- q[q$nes!= 'Inf',]
q <- q[q$Sample == 'LUNG_H_Cluster6_vs_Rest',]

#q <- rbind(head(q, 8), tail(q,8))
q$Term <- gsub("_"," ",q$Term)
q$Term <- gsub("HALLMARK","",q$Term)
q$nes <- as.numeric(q$nes)
GSEA_figure <- 
  ggplot(q,aes(y=nes,x=reorder(Term,nes),label = Term))+
  #geom_text(aes(y = yax,hjust = nes > 0),size=(figure_font_size)/(14/5))+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-3,3))+
  labs(y = paste0("Normalised enrichment score" ), x = "",title="")+
  scale_fill_manual(values = c("#7A378B","#7F7F7F"),name=paste0("Group: "),labels=c("Cluster 6","Rest"))+ # red, blue, "#2E8B57","#1874CD"
  #scale_fill_manual(values = c("#6A81A0","#E99A98"),name=paste0("Group: "),labels=c("Resistant","Sensitive"))+ # red, blue, "#2E8B57","#1874CD"
  scale_x_discrete(breaks = NULL)+
  theme_manuscript(base_size = figure_font_size)+
  theme(panel.grid.major = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(0.1,0,0.05,0), "lines"))
png('/Users/asanchezb/Desktop/test2.png',width=1400, height=1400, res = 280)
plot(GSEA_figure)
dev.off() 

#scale_fill_manual(values = c("#6A81A0","#E99A98"),
#scale_fill_manual(values=c("#7F7F7F", '#E0E0E0',"#7A378B")) +


#############################################################################################################
# 7. BOXPLOT
#############################################################################################################

png(filename = '/Users/asanchezb/Desktop/test.png', width = 1300, height = 1200, res = 250)#1300#2500
ggplot(data=df_merged, aes(x=class, y=diff, col=class)) + 
  geom_point(alpha = 0.7, size = 1.3) +
  theme_minimal() +
  # box.padding = 0.5,max.overlaps = Inf
  scale_color_manual(values = cols) +
  geom_text_repel(
    data = df_labels,
    aes(label = Gene),
    size = 2,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50"
  ) +
  geom_hline(yintercept=c(0), col="black", size = 0.3) +
  #geom_hline(yintercept=-log10(0.05), col="black") +
  theme(panel.grid = element_blank(),
        plot.margin = margin(15, 15, 15, 15, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))+
  ylim(c(-2.3,2.3))+
  guides(label = FALSE, color=FALSE) 
dev.off()

