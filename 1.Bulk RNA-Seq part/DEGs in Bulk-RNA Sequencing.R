
# loading packages
library(BiocManager)
library(limma)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(ggsci)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggalt) # install.packages("ggalt")
library(ggpubr)
library(VennDiagram)
library(viridis)
library(writexl)

data <- read.csv("E:/Bulk-RNA_Seq/data.csv",row.names = 1)
colnames(data)

#### DEGs in Bulk-RNA (Control vs COPD) ####

counts <- data[,c(1:6)] # Control vs COPD
head(counts)

Controldition = factor(c(rep("Control",3),rep("COPD",3)),levels = c("Control","COPD"))
Controldition

genelist = DGEList(counts = counts, group = Controldition)
design <- model.matrix(~Controldition)
colnames(design) <- levels(Controldition)
rownames(design) <- colnames(counts)
design

keep <- rowSums(cpm(genelist)>1)>=2 
genelist.filted <- genelist[keep,,keep.lib.sizes=FALSE] 
genelist.filted

genelist.norm <- calcNormFactors(genelist.filted)  
logCPM <- cpm(genelist.norm, log=TRUE, prior.count=3)

write.csv(logCPM, file = "E:/Bulk-RNA_Seq/logCPM_(Control vs COPD).csv")


fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
DEG_result <- topTable(fit,coef=2,n=Inf)
DEG_result
DEG_result$type = ifelse(DEG_result$P.Value < 0.05 & abs(DEG_result$logFC) >=1,
                         ifelse(DEG_result$logFC > 1 ,'Up(647)','Down(150)'),'Stable')
table(DEG_result$type)
#Down(150)    Stable   Up(647) 
#     150     14603       647 

write.csv(DEG_result,"DEGs in Bulk-RNA Seq(COPD vs Control).csv" )

DEG_up <- row.names(DEG_result[DEG_result$type == "Up(647)",])  
head(DEG_up)
DEG_down <- row.names(DEG_result[DEG_result$type == "Down(150)",]) 
head(DEG_down)

p1 <- ggplot(data=DEG_result, aes(x=logFC, y =-log10(P.Value),colour=type)) +
  geom_point(alpha=0.5, size=1)+
  scale_color_manual(values=c("#527eba","grey","#dc2a1c"))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.8,alpha=0.8)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.8,alpha=0.8)+
  labs(x="log2 (fold change)",y="-log10 (FDR)")+
  ggtitle("Signifianct DEGs in Bulk RNA-seq \n (COPD vs Control)")+ 
  theme_bw() + xlim(-5,5) + ylim(0,5) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) 
p1

#### DEGs in Bulk-RNA (QFGB vs COPD) ####
data <- read.csv("E:/Bulk-RNA_Seq/data.csv",row.names = 1)
colnames(data)
counts <- data[,c(1:6)] # QFGB vs COPD
head(counts)

Controldition = factor(c(rep("QFGB",3),rep("COPD",3)),levels = c("QFGB","COPD"))
Controldition

genelist = DGEList(counts = counts, group = Controldition)
design <- model.matrix(~Controldition)
colnames(design) <- levels(Controldition)
rownames(design) <- colnames(counts)
design

keep <- rowSums(cpm(genelist)>1)>=2 
genelist.filted <- genelist[keep,,keep.lib.sizes=FALSE] 
genelist.filted

genelist.norm <- calcNormFactors(genelist.filted)  
logCPM <- cpm(genelist.norm, log=TRUE, prior.count=3)

write.csv(logCPM, file = "E:/Bulk-RNA_Seq/logCPM_(QFGB vs COPD).csv")

fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
DEG_result <- topTable(fit,coef=2,n=Inf)
DEG_result
DEG_result$type = ifelse(DEG_result$P.Value < 0.05 & abs(DEG_result$logFC) >= 1,
                         ifelse(DEG_result$logFC > 1 ,'Up(305)','Down(75)'),'Stable')
table(DEG_result$type)
# Down(75)   Stable  Up(305) 
#      75    14973      305 

write.csv(DEG_result,"DEGs in Bulk-RNA Seq(COPD vs QFGB).csv" )

DEG_up <- row.names(DEG_result[DEG_result$type == "Up(305)",])  
head(DEG_up)
DEG_down <- row.names(DEG_result[DEG_result$type == "Down(75)",]) 
head(DEG_down)

#DEG_result$label <- ifelse(DEG_result$adj.P.Val < 0.05 & abs(DEG_result$logFC) >= 1,rownames(DEG_result),"")

p2 <- ggplot(data=DEG_result, aes(x=logFC, y =-log10(P.Value),colour=type)) +
  geom_point(alpha=0.5, size=1)+
  scale_color_manual(values=c("#527eba","grey","#dc2a1c"))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.8,alpha=0.8)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.8,alpha=0.8)+
  labs(x="log2 (fold change)",y="-log10 (FDR)")+
  ggtitle("Signifianct DEGs in Bulk RNA-seq \n (COPD vs QFGB)")+ 
  theme_bw() + xlim(-3,3) + ylim(0, 6) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) 
p2


####UP Controlcordance of DEGs in Bulk-RNA ####
DEGs_COPD_vs_Control <- read.csv("E:/Bulk-RNA_Seq/DEGs in Bulk-RNA Seq(COPD vs Control).csv",row.names = 1)
table(DEGs_COPD_vs_Control$type)

DEGs_COPD_vs_Control_up <- rownames(DEGs_COPD_vs_Control[DEGs_COPD_vs_Control$type == "Up(647)",])
DEGs_COPD_vs_Control_down <- rownames(DEGs_COPD_vs_Control[DEGs_COPD_vs_Control$type == "Down(150)",])

DEGs_COPD_vs_QFGB <- read.csv("E:/Bulk-RNA_Seq/DEGs in Bulk-RNA Seq(COPD vs QFGB).csv",row.names = 1)
table(DEGs_COPD_vs_QFGB$type)

DEGs_COPD_vs_QFGB_up <- rownames(DEGs_COPD_vs_QFGB[DEGs_COPD_vs_QFGB$type == "Up(305)",])
DEGs_COPD_vs_QFGB_down <- rownames(DEGs_COPD_vs_QFGB[DEGs_COPD_vs_QFGB$type == "Down(75)",])

venn.plot <- 
  venn.diagram(
    x = list(
      'DEGs_COPD_vs_Control_up (647)' = DEGs_COPD_vs_Control_up,
      'DEGs_COPD_vs_QFGB_up (305)' = DEGs_COPD_vs_QFGB_up
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
    main = "Overlap up-regulated DEGs",
    main.cex = 1.2)

pdf(file="Overlap up-regulated DEGs.pdf")
grid.draw(venn.plot)
dev.off()

Overlap_up_regulated_DEGs_list <- intersect(DEGs_COPD_vs_Control_up, DEGs_COPD_vs_QFGB_up)
Overlap_up_regulated_DEGs_df <- data.frame(Compound = Overlap_up_regulated_DEGs_list)
write_xlsx(list("Overlap_up_regulated_DEGs" = Overlap_up_regulated_DEGs_df), 
           path = "Overlap_up_regulated_DEGs.xlsx")
Overlap_up_regulated_DEGs_list <- intersect(DEGs_COPD_vs_Control_up,DEGs_COPD_vs_QFGB_up)
write.csv(Overlap_up_regulated_DEGs_list, "E:/Bulk-RNA_Seq/overlap_up_genes.csv", row.names = FALSE)

# UP-Pheatmap
data <- read.csv("E:/Bulk-RNA_Seq/overlap_up_genes.csv",row.names = 1)
colnames(data)
counts <- data[,c(1:9)] # Control vs COPD
head(counts)

Controldition = factor(c(rep("Control",3),rep("COPD",3),rep("QFGB",3)),levels = c("Control","COPD","QFGB"))
Controldition

genelist = DGEList(counts = counts, group = Controldition)
design <- model.matrix(~Controldition)
colnames(design) <- levels(Controldition)
rownames(design) <- colnames(counts)
design

keep <- rowSums(cpm(genelist)>1)>=2 
genelist.filted <- genelist[keep,,keep.lib.sizes=FALSE] 
genelist.filted

genelist.norm <- calcNormFactors(genelist.filted)  
logCPM <- cpm(genelist.norm, log=TRUE, prior.count=3)

pheatmap::pheatmap(data[Overlap_up_regulated_DEGs_list,],
                   scale = "row",
                   cluster_cols = F,
                   treeheight_col = 15,
                   treeheight_row = 15,
                   cellheight = 0.5,
                   cellwidth = 10,
                   fontsize_row = 5,
                   fontsize_col = 5,
                   angle_col = 45,
                   show_rownames = F,
                   border_color = "grey90")

####Down  Controlcordance of DEGs in Bulk-RNA ####

DEGs_COPD_vs_Control <- read.csv("E:/Bulk-RNA_Seq/DEGs in Bulk-RNA Seq(COPD vs Control).csv",row.names = 1)
table(DEGs_COPD_vs_Control$type)

DEGs_COPD_vs_Control_up <- rownames(DEGs_COPD_vs_Control[DEGs_COPD_vs_Control$type == "Up(647)",])
DEGs_COPD_vs_Control_down <- rownames(DEGs_COPD_vs_Control[DEGs_COPD_vs_Control$type == "Down(150)",])

DEGs_COPD_vs_QFGB <- read.csv("E:/Bulk-RNA_Seq/DEGs in Bulk-RNA Seq(COPD vs QFGB).csv",row.names = 1)
table(DEGs_COPD_vs_QFGB$type)

DEGs_COPD_vs_QFGB_up <- rownames(DEGs_COPD_vs_QFGB[DEGs_COPD_vs_QFGB$type == "Up(305)",])
DEGs_COPD_vs_QFGB_down <- rownames(DEGs_COPD_vs_QFGB[DEGs_COPD_vs_QFGB$type == "Down(75)",])

venn.plot <- venn.diagram(
  x = list(
    'COPD vs Control Down' = DEGs_COPD_vs_Control_down,
    'COPD vs QFGB Down' = DEGs_COPD_vs_QFGB_down
  ),
  filename = NULL,
  col = "black",
  fill = c("#4CBBD5", "#F39B7F"),
  alpha = 0.8,
  cex = 0.8,
  cat.col = 'black',
  cat.cex = 0.8,
  cat.fontface = "bold",
  margin = 0.05,
  main = "Overlap down-regulated DEGs",
  main.cex = 1.2
)


Overlap_down_regulated_DEGs_list <- intersect(DEGs_COPD_vs_Control_down,DEGs_COPD_vs_QFGB_down)
write.csv(Overlap_down_regulated_DEGs_list, "E:/Bulk-RNA_Seq/overlap_down_genes.csv", row.names = FALSE)


# Down-Pheatmap
data <- read.csv("E:/Bulk-RNA_Seq/overlap_down_genes.csv",row.names = 1)
counts <- data[,c(1:9)] # Control vs COPD
head(counts)

Controldition = factor(c(rep("Control",3),rep("COPD",3),rep("QFGB",3)),levels = c("Control","COPD","QFGB"))
Controldition

genelist = DGEList(counts = counts, group = Controldition)
design <- model.matrix(~Controldition)
colnames(design) <- levels(Controldition)
rownames(design) <- colnames(counts)
design

keep <- rowSums(cpm(genelist)>1)>=2 
genelist.filted <- genelist[keep,,keep.lib.sizes=FALSE] 
genelist.filted

genelist.norm <- calcNormFactors(genelist.filted)  
logCPM <- cpm(genelist.norm, log=TRUE, prior.count=3)

pheatmap::pheatmap(data[Overlap_down_regulated_DEGs_list,],
                   scale = "row",
                   cluster_cols = F,
                   treeheight_col = 15,
                   treeheight_row = 15,
                   cellheight = 9,
                   cellwidth = 13,
                   fontsize_row = 5,
                   fontsize_col = 5,
                   angle_col = 45,
                   show_rownames = F,
                   border_color = "grey90")
# pheatmap::pheatmap(data[Overlap_down_regulated_DEGs_list,],scale = "row",cluster_cols = F)


