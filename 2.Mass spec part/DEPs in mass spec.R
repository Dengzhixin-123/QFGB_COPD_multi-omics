
# loading packages
library(edgeR) 
library(limma)
library(ggplot2)
library(ggrepel)
library(DMwR2)
library(mice)
library(clusterProfiler)
library(org.Mm.eg.db)
library(patchwork)
library(RColorBrewer)
library(ggrepel)
library(ggpubr)
library(VennDiagram)
library(writexl)

#### DEPs analysis (Control vs COPD) ####
raw_data <- read.csv(file.choose(), header = TRUE, check.names = FALSE)
head(raw_data)
trim_data <- raw_data[rowSums(is.na(raw_data)) <= 1,]
md.pattern(raw_data)
colnames(trim_data)
trim_data <- trim_data[,c(1,5:10)]
head(trim_data)

list <- which(duplicated(trim_data$"Gene Name"))
trim_data = trim_data[-list,]
rownames(trim_data) <- trim_data$"Gene Name"

colnames(trim_data)

DMSO <- trim_data[,c(5:7)]
DMSO <- knnImputation(DMSO)
head(DMSO)

celp <- trim_data[,c(2:4)]
celp <- knnImputation(celp)
head(celp)

trim_data <- cbind(DMSO,celp)
md.pattern(trim_data)

head(trim_data)

write.csv(trim_data,"Protein expression in Label free.csv")

expr_data  <- log2(trim_data + 1)

write.csv(data, file = "E:/Proteomics/log2+1(CON vs COPD).csv")

group <- factor(c(rep("Control", 3), rep("COPD", 3)), levels = c("Control", "COPD"))

design <- model.matrix(~group)
colnames(design) <- levels(group)        
rownames(design) <- colnames(data)    

fit <- lmFit(expr_data , design)
fit <- eBayes(fit, trend = TRUE)
deg <- topTable(fit, coef = 2, n = Inf)
cut_off_pvalue <- 0.05
cut_off_logFC <- 0.263 

deg$change <- ifelse(deg$P.Value < cut_off_pvalue & abs(deg$logFC) >= cut_off_logFC,
                     ifelse(deg$logFC > 0, 'Up(423)','Down(185)'),'Stable')

table(deg$change)
#Down(185)    Stable   Up(423) 
#     185      7546       423 

write.csv(deg, "DEP_result_DEP_COPD_vs_Control.csv", quote = FALSE)

p1 <- ggplot(data=deg, aes(x=logFC, y =-log10(P.Value),colour=change)) +
  geom_point(alpha=0.5, size=1.5)+
  scale_color_manual(values=c("#527eba","grey","#dc2a1c"))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.8,alpha=0.8)+
  geom_vline(xintercept = c(0.263,-0.263),lty=4,lwd=0.8,alpha=0.8)+
  labs(x="log2 (fold change)",y="-log10 (P)")+
  ggtitle("Control vs COPD")+ 
  theme_bw() + xlim(-3,3) + ylim(0,5) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) 
p1


#### DEPs analysis (QFGB vs COPD) ####
raw_data <- read.csv(file.choose(), header = TRUE, check.names = FALSE)
head(raw_data)
trim_data <- raw_data[rowSums(is.na(raw_data)) <= 1,]
md.pattern(raw_data) 
colnames(trim_data)
trim_data <- trim_data[,c(1:7)]
head(trim_data)

rownames(trim_data) <- make.unique(as.character(trim_data$`Gene Name`))

colnames(trim_data)

DMSO <- trim_data[,c(2:4)]
DMSO <- knnImputation(DMSO)  
head(DMSO)

celp <- trim_data[,c(5:7)]
celp <- knnImputation(celp)  
head(celp)

trim_data <- cbind(DMSO,celp)
md.pattern(trim_data)

#trim_data <- trim_data[,-1]
head(trim_data)

write.csv(trim_data,"QFGBvsCOPD_Protein expression in Label free.csv")

expr_data  <- log2(trim_data + 1)

write.csv(data, file = "E:/Proteomics/log2+1(QFGB vs COPD).csv")

group <- factor(c(rep("QFGB", 3), rep("COPD", 3)), levels = c("QFGB", "COPD"))

design <- model.matrix(~group)
colnames(design) <- levels(group)       
rownames(design) <- colnames(data)       
fit <- lmFit(expr_data , design)
fit <- eBayes(fit, trend = TRUE)
deg <- topTable(fit, coef = 2, n = Inf)
cut_off_pvalue <- 0.05
cut_off_logFC <- 0.263  

deg$change <- ifelse(deg$P.Value < cut_off_pvalue & abs(deg$logFC) >= cut_off_logFC,
                     ifelse(deg$logFC > 0, 'Up(264)','Down(158)'),'Stable')

table(deg$change)
#Down(158)    Stable   Up(264) 
#     158      7732       264 

write.csv(deg, "DEP_result_DEP_COPD_vs_QFGB.csv", quote = FALSE)

p2 <- ggplot(data=deg, aes(x=logFC, y =-log10(P.Value),colour=change)) +
  geom_point(alpha=0.5, size=1.5)+
  scale_color_manual(values=c("#527eba","grey","#dc2a1c"))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.8,alpha=0.8)+
  geom_vline(xintercept = c(0.263,-0.263),lty=4,lwd=0.8,alpha=0.8)+
  labs(x="log2 (fold change)",y="-log10 (P)")+
  ggtitle("QFGB vs COPD")+ 
  theme_bw() + xlim(-3,3) + ylim(0,4.5) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) 
p2

####UP in Protein ###
DEP_COPD_vs_Control <- read.csv("E:/Proteomics/DEP_result_DEP_COPD_vs_Control.csv",row.names = 1)
table(DEP_COPD_vs_Control$change)

DEP_COPD_vs_Control_up <- rownames(DEP_COPD_vs_Control[DEP_COPD_vs_Control$change == "Up",])
DEP_COPD_vs_Control_down <- rownames(DEP_COPD_vs_Control[DEP_COPD_vs_Control$change == "Down",])

DEP_COPD_vs_QFGB <- read.csv("E:/Proteomics/DEP_result_DEP_COPD_vs_QFGB.csv",row.names = 1)
table(DEP_COPD_vs_QFGB$change)

DEP_COPD_vs_Control_up <- rownames(DEP_COPD_vs_Control[grepl("^Up", DEP_COPD_vs_Control$change), ])
DEP_COPD_vs_QFGB_up    <- rownames(DEP_COPD_vs_QFGB[grepl("^Up", DEP_COPD_vs_QFGB$change), ])

venn.plot <- 
  venn.diagram(
    x = list(
      'DEP_COPD_vs_Control_up' = DEP_COPD_vs_Control_up,
      'DEP_COPD_vs_QFGB_up' = DEP_COPD_vs_QFGB_up
    ),
    filename = NULL,
    col = "black",
    fill = c("#4CBBD5","#F39B7F"),
    alpha = 0.8,
    cex = 0.8,
    cat.col = 'black',
    cat.cex = 0.8,
    cat.fontface = "bold",
    margin = 0.05,
    main = "Overlap up-regulated DEP",
    main.cex = 1.2)

pdf(file="Overlap up-regulated DEP.pdf")
grid.draw(venn.plot)
dev.off()

Overlap_up_regulated_DEP_list <- intersect(DEP_COPD_vs_Control_up, DEP_COPD_vs_QFGB_up)
Overlap_up_regulated_DEP_df <- data.frame(Compound = Overlap_up_regulated_DEP_list)
write_xlsx(list("Overlap_up_regulated_DEP" = Overlap_up_regulated_DEP_df), 
           path = "Overlap_up_regulated_DEP.xlsx")

heatmap_up <- Overlap_up_regulated_DEP_df
head(heatmap_up)

### Down in Protein ###
DEP_COPD_vs_Control <- read.csv("E:/Proteomics/DEP_result_DEP_COPD_vs_Control.csv",row.names = 1)
table(DEP_COPD_vs_Control$change)

DEP_COPD_vs_Control_up <- rownames(DEP_COPD_vs_Control[DEP_COPD_vs_Control$change == "Up",])
DEP_COPD_vs_Control_down <- rownames(DEP_COPD_vs_Control[DEP_COPD_vs_Control$change == "Down",])

DEP_COPD_vs_QFGB <- read.csv("E:/Proteomics/DEP_result_DEP_COPD_vs_QFGB.csv",row.names = 1)
table(DEP_COPD_vs_QFGB$change)

DEP_COPD_vs_Control_down <- rownames(DEP_COPD_vs_Control[grepl("^Down", DEP_COPD_vs_Control$change), ])
DEP_COPD_vs_QFGB_down    <- rownames(DEP_COPD_vs_QFGB[grepl("^Down", DEP_COPD_vs_QFGB$change), ])

venn.plot <- 
  venn.diagram(
    x = list(
      'DEP_COPD_vs_Control_down' = DEP_COPD_vs_Control_down,
      'DEP_COPD_vs_QFGB_down' = DEP_COPD_vs_QFGB_down
    ),
    filename = NULL,
    col = "black",
    fill = c("#4CBBD5","#F39B7F"),
    alpha = 0.8,
    cex = 0.8,
    cat.col = 'black',
    cat.cex = 0.8,
    cat.fontface = "bold",
    margin = 0.05,
    main = "Overlap down-regulated DEP",
    main.cex = 1.2)

pdf(file="Overlap down-regulated DEP.pdf")
grid.draw(venn.plot)
dev.off()

Overlap_down_regulated_DEP_list <- intersect(DEP_COPD_vs_Control_down, DEP_COPD_vs_QFGB_down)
Overlap_down_regulated_DEP_df <- data.frame(Compound = Overlap_up_regulated_DEP_list)
write_xlsx(list("Overlap_down_regulated_DEP" = Overlap_up_regulated_DEP_df), 
           path = "Overlap_down_regulated_DEP.xlsx")

### Down in Protein ###
DEP_COPD_vs_Control <- read.csv("E:/Proteomics/DEP_result_DEP_COPD_vs_Control.csv",row.names = 1)
table(DEP_COPD_vs_Control$change)

DEP_COPD_vs_Control_down <- rownames(DEP_COPD_vs_Control[grepl("^Down", DEP_COPD_vs_Control$change), ])
DEP_COPD_vs_QFGB_down    <- rownames(DEP_COPD_vs_QFGB[grepl("^Down", DEP_COPD_vs_QFGB$change), ])

DEP_COPD_vs_QFGB <- read.csv("E:/Proteomics/DEP_result_DEP_COPD_vs_QFGB.csv",row.names = 1)
table(DEP_COPD_vs_QFGB$change)

DEP_COPD_vs_Control_down <- rownames(DEP_COPD_vs_Control[grepl("^Down", DEP_COPD_vs_Control$change), ])
DEP_COPD_vs_QFGB_down    <- rownames(DEP_COPD_vs_QFGB[grepl("^Down", DEP_COPD_vs_QFGB$change), ])

venn.plot <- 
  venn.diagram(
    x = list(
      'DEP_COPD_vs_Control_down' = DEP_COPD_vs_Control_down,
      'DEP_COPD_vs_QFGB_down' = DEP_COPD_vs_QFGB_down
    ),
    filename = NULL,
    col = "black",
    fill = c("#4CBBD5","#F39B7F"),
    alpha = 0.8,
    cex = 0.8,
    cat.col = 'black',
    cat.cex = 0.8, 
    cat.fontface = "bold",
    margin = 0.05,
    main = "Overlap down-regulated DEP",
    main.cex = 1.2)

pdf(file="Overlap down-regulated DEP.pdf")
grid.draw(venn.plot)
dev.off()

Overlap_down_regulated_DEP_list <- intersect(DEP_COPD_vs_Control_down, DEP_COPD_vs_QFGB_down)
Overlap_down_regulated_DEP_df <- data.frame(Compound = Overlap_down_regulated_DEP_list)
write_xlsx(list("Overlap_down_regulated_DEP" = Overlap_down_regulated_DEP_df), 
           path = "Overlap_down_regulated_DEP.xlsx")

# UP-Pheatmap
data_up <- read.csv(file.choose(),header = T,row.names = 1,check.names = F)
pheatmap::pheatmap(data_up,
                   scale = "row",
                   cluster_cols = F,
                   treeheight_col = 15,
                   treeheight_row = 15,
                   cellheight = 5,
                   cellwidth = 23,
                   fontsize_row = 5,
                   fontsize_col = 7,
                   angle_col = 45,
                   show_rownames = F,
                   border_color = "grey90")

# Down-Pheatmap
data_down <- read.csv(file.choose(),header = T,row.names = 1,check.names = F)
pheatmap::pheatmap(data_down,
                   scale = "row",
                   cluster_cols = F,
                   treeheight_col = 15,
                   treeheight_row = 15,
                   cellheight = 16,
                   cellwidth = 15,
                   fontsize_row = 10,
                   fontsize_col = 10,
                   angle_col = 45,
                   show_rownames = T,
                   border_color = "grey90")

