library(tidyverse)
library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(sctransform)
library(hdf5r)
library(ggrepel)
require(toolboxH)
library(ComplexHeatmap)
library(clusterProfiler)
library(enrichplot)
require(DOSE)
library(pathview)
library(org.Mm.eg.db)
library(writexl)
library(ggthemes)

ScSeuMerged <- readRDS("Data/s565_1_seurat_scRNA_zenodo.rds")

ScSeuMerged@active.ident <-factor(ScSeuMerged$celltype_v3)
reihenfolge = c('AM', 'Macrophages+DC', 'Neutrophils', 'TNK_cells', 'B_cells', 'AT2', 'Endothelial', 'Fibroblasts')
UMAPPlot(ScSeuMerged, label = T)

NogpaletteReihe <-  c("#CB769E", "#A85C85", "#0081AF", "#368F8B", "#62C370", "#F97E44", "#0D3B66", "#B2675E")

my_levels <- c(reihenfolge)
ScSeuMerged@active.ident <- factor(x = ScSeuMerged@active.ident, levels = my_levels)

UMAPPlot(ScSeuMerged, label = T, cols = c(NogpaletteReihe)) + NoLegend() + (coord_fixed(ratio=1))

plot <- DotPlot(ScSeuMerged, features = c("Adgre1", "Marco", "Flt3", "S100a8", "Cd3d", "Nkg7", "Cd79a", "Epcam", "Lamp3", "Pecam1", "Col1a1", "Inmt"))
plot2 <- plot + theme(axis.text.x = element_text(angle = 45, hjust=1))
plot2

Subsets <- subset(ScSeuMerged, idents = c("AT2", "Endothelial"))
UMAPPlot(Subsets, group.by = "Gruppe", cols = c("darkgrey", "black", "#0099a1", "#285D5E")) + NoLegend()

####Differential expression analysis
N8_AT2_edgr <- read_excel2("Data/Supplementary data CXCL5.xlsx", 3)
N8_endo_edgr <- read_excel2("Data/Supplementary data CXCL5.xlsx", 5)

N8_AT2_edgrPBS <- read_excel2("Data/Supplementary data CXCL5.xlsx", 2)
N8_endo_edgrPBS <- read_excel2("Data/Supplementary data CXCL5.xlsx", 4)


N8_AT2_edgr$Significant <- "FALSE"
N8_AT2_edgr$Significant[N8_AT2_edgr$PValue < 0.05] <- "TRUE"

N8_endo_edgr$Significant <- "FALSE"
N8_endo_edgr$Significant[N8_endo_edgr$PValue < 0.05] <- "TRUE"

N8_AT2_edgrPBS$Significant <- "FALSE"
N8_AT2_edgrPBS$Significant[N8_AT2_edgrPBS$PValue < 0.05] <- "TRUE"

N8_endo_edgrPBS$Significant <- "FALSE"
N8_endo_edgrPBS$Significant[N8_endo_edgrPBS$PValue < 0.05] <- "TRUE"

N8_endo_edgr$diffexpressed <- "FALSE"
N8_endo_edgr$diffexpressed[N8_endo_edgr$logFC > 0 & N8_endo_edgr$PValue < 0.05] <- "UP"
N8_endo_edgr$diffexpressed[N8_endo_edgr$logFC < 0 & N8_endo_edgr$PValue < 0.05] <- "DOWN"

N8_AT2_edgr$diffexpressed <- "FALSE"
N8_AT2_edgr$diffexpressed[N8_AT2_edgr$logFC > 0 & N8_AT2_edgr$PValue < 0.05] <- "UP"
N8_AT2_edgr$diffexpressed[N8_AT2_edgr$logFC < 0 & N8_AT2_edgr$PValue < 0.05] <- "DOWN"

N8_endo_edgrPBS$diffexpressed <- "FALSE"
N8_endo_edgrPBS$diffexpressed[N8_endo_edgrPBS$logFC > 0 & N8_endo_edgrPBS$PValue < 0.05] <- "UP"
N8_endo_edgrPBS$diffexpressed[N8_endo_edgrPBS$logFC < 0 & N8_endo_edgrPBS$PValue < 0.05] <- "DOWN"

N8_AT2_edgrPBS$diffexpressed <- "FALSE"
N8_AT2_edgrPBS$diffexpressed[N8_AT2_edgrPBS$logFC > 0 & N8_AT2_edgrPBS$PValue < 0.05] <- "UP"
N8_AT2_edgrPBS$diffexpressed[N8_AT2_edgrPBS$logFC < 0 & N8_AT2_edgrPBS$PValue < 0.05] <- "DOWN"

AT2 <- ggplot(data=N8_AT2_edgr, aes(x=logFC, y=-log10(PValue),col=diffexpressed, label=symbol)) + geom_point(size = 4) + theme_classic() + geom_text_repel(size = 10, max.overlaps = 4) + geom_vline(xintercept=c(-1, 1), linetype = "dashed") + geom_hline(yintercept=-log10(0.05), linetype = "dashed") + scale_colour_manual(values=c("blue", "black", "red")) + NoLegend()
AT2

N8_AT2_edgr1 <- N8_AT2_edgr %>% 
  filter(symbol %in% c("Egfr", "Vcl", "Stxbp6", "Egfr", "Lrg1", "Lrp4"))

p1 <- ggplot(data=N8_AT2_edgr, aes(x=logFC, y=-log10(PValue),col=diffexpressed, label=symbol)) + geom_point(size = 4) + geom_point(data=N8_AT2_edgr1, aes(x=logFC,y=-log10(PValue)), col="pink", size = 4) + theme_classic() + geom_text_repel(size = 10, max.overlaps = 4) + geom_vline(xintercept=c(-1, 1), linetype = "dashed") + geom_hline(yintercept=-log10(0.05), linetype = "dashed") + scale_colour_manual(values=c("blue", "black", "red"))  + geom_text_repel(data = N8_AT2_edgr1 %>%   filter(symbol %in% c("Egfr", "Vcl", "Stxbp6", "Egfr", "Lrg1", "Lrp4")), size = 4) + NoLegend()
p1

Endo <- ggplot(data=N8_endo_edgr, aes(x=logFC, y=-log10(PValue),col=diffexpressed, label=symbol)) + geom_point(size = 4) + theme_classic() + geom_text_repel(size = 10, max.overlaps = 4) + geom_vline(xintercept=c(-1, 1), linetype = "dashed") + geom_hline(yintercept=-log10(0.05), linetype = "dashed") + scale_colour_manual(values=c("blue", "black", "red")) + NoLegend()
Endo

N8_endo_edgr1 <- N8_endo_edgr %>% 
  filter(symbol %in% c("Tjp1", "Vegfa", "Lasp1"))

p2 <- ggplot(data=N8_endo_edgr, aes(x=logFC, y=-log10(PValue),col=diffexpressed, label=symbol)) + geom_point(size = 4) + geom_point(data=N8_endo_edgr1, aes(x=logFC,y=-log10(PValue)), col="pink", size = 4) + theme_classic() + geom_text_repel(size = 10, max.overlaps = 4) + geom_vline(xintercept=c(-1, 1), linetype = "dashed") + geom_hline(yintercept=-log10(0.05), linetype = "dashed") + scale_colour_manual(values=c("blue", "black", "red"))  + geom_text_repel(data = N8_endo_edgr1 %>%   filter(symbol %in% c("Tjp1", "Vegfa", "Lasp1")), size = 8) + NoLegend()
p2

AT2PBS <- ggplot(data=N8_AT2_edgrPBS, aes(x=logFC, y=-log10(PValue),col=diffexpressed, label=symbol)) + geom_point(size = 4) + theme_classic() + geom_text_repel(size = 10, max.overlaps = 4) + geom_vline(xintercept=c(-1, 1), linetype = "dashed") + geom_hline(yintercept=-log10(0.05), linetype = "dashed") + scale_colour_manual(values=c("blue", "black", "red")) + NoLegend()
AT2PBS

EndoPBS <- ggplot(data=N8_endo_edgrPBS, aes(x=logFC, y=-log10(PValue),col=diffexpressed, label=symbol)) + geom_point(size = 4) + theme_classic() + geom_text_repel(size = 8, max.overlaps = 3) + geom_vline(xintercept=c(-1, 1), linetype = "dashed") + geom_hline(yintercept=-log10(0.05), linetype = "dashed") + scale_colour_manual(values=c("blue", "black", "red")) + NoLegend()
EndoPBS



################Complex Heatmap############################################################################################

HeatmapAT2 <-filter(N8_AT2_edgr, symbol %in% c("Egf", "Egfr", "Eif2ak3", "Lamp1", "Scnn1a", "Scnn1b", "Scnn1g", "Cldn3", "Cldn4", "Cldn8", "Cldn18", "Clcn2", "Slc26a9", "Slc12a2", "Cftr", "Clcn1", "Clcn2", "Clcn3b", "Clcn4a", "Clcn4b", "Clcnka", "Clcnkb", "Bsnd", "Clca1", "Clca2", "Clca3b", "Clca4a", "Ano1", "Ano2", "Best1", "Best2", "Best3", "Slc12a1", "Slc12a4", "Slc12a5", "Slc4a8", "Slc4a10", "Kcnk2", "Mapk8", "Mapk9", "Mapk10", "Nrg1", "Erbb3", "Erbb2", "Tjp1",  "F11r", "Gja3", "Gja8", "Gja1", "Kirrel", "Afdn", "Tjp3", "Tjp2", "Ocln", "Itgb4", "Ctnnb1", "Cdh1",  "Kcnk2", "Kcnk3", "Kcnk4", "Kcnk5", "Kcnk6", "Kcnk7", "Kcnk9", "Kcnk10", "Kcnk12", "Kcnk13", "Kcnk15", "Kcnk16", "Kcnk18", "Jam2", "Jam3", "Map2k1", "Pard3", "Pard6b", "Pten"))

HeatmapAT2select <-filter(N8_AT2_edgr, symbol %in% c("Cldn18", "Cldn3", "Lamp1", "Slc12a2", "Scnn1a", "Ctnnb1", "Pard3", "Afdn", "Cdh1", "Tjp2", "Tjp1", "F11r", "Pten", "Gja1", "Slc26a9", "Map2k1", "Erbb3", "Cftr", "Ocln", "Pard6b", "Eif2ak3", "Tjp3", "Mapk8", "Mapk9", "Scnn1g", "Egfr", "Cldn4", "Kcnk2", "Ano1", "Erbb2")) 

HeatmapAT2 <- HeatmapAT2select %>% column_to_rownames(var="symbol")
my_AT2_matrix <- as.matrix(HeatmapAT2[,c(13:18)])

Heatmap(my_AT2_matrix, column_order = c("AT2 A1", "AT2 A2", "AT2 A3", "AT2 B1", "AT2 B2", "AT2 B3"))

HeatmapEndo <-filter(N8_endo_edgr, symbol %in% c("Egf", "Egfr", "Eif2ak3", "Lamp1", "Scnn1a", "Scnn1b", "Scnn1g", "Cldn3", "Cldn4", "Cldn8", "Cldn18", "Clcn2", "Slc26a9", "Slc12a2", "Cftr", "Clcn1", "Clcn2", "Clcn3b", "Clcn4a", "Clcn4b", "Clcnka", "Clcnkb", "Bsnd", "Clca1", "Clca2", "Clca3b", "Clca4a", "Ano1", "Ano2", "Best1", "Best2", "Best3", "Slc12a1", "Slc12a4", "Slc12a5", "Slc4a8", "Slc4a10", "Kcnk2", "Mapk8", "Mapk9", "Mapk10", "Nrg1", "Erbb3", "Erbb2", "Tjp1",  "F11r", "Gja3", "Gja8", "Gja1", "Kirrel", "Afdn", "Tjp3", "Tjp2", "Ocln", "Itgb4", "Ctnnb1", "Cdh1",  "Kcnk2", "Kcnk3", "Kcnk4", "Kcnk5", "Kcnk6", "Kcnk7", "Kcnk9", "Kcnk10", "Kcnk12", "Kcnk13", "Kcnk15", "Kcnk16", "Kcnk18", "Jam2", "Jam3", "Map2k1", "Pard3", "Pard6b", "Pten"))

HeatmapEndoselect <-filter(N8_endo_edgr, symbol %in% c("Tjp1", "Lamp1", "Pard3", "Afdn", "Ctnnb1", "Pten", "Tjp2", "F11r", "Mapk9", "Jam2", "Slc12a2", "Map2k1", "Slc12a4", "Eif2ak3", "Ocln", "Mapk8", "Cldn18", "Erbb2", "Kcnk6", "Kcnk5", "Scnn1a", "Jam3", "Itgb4", "Ano2", "Gja1", "Slc4a8", "Scnn1b", "Kirrel", "Nrg1", "Pard6b")) 

colnames(HeatmapEndoselect)[14] <- "ECs A1"
colnames(HeatmapEndoselect)[15] <- "ECs A2"
colnames(HeatmapEndoselect)[16] <- "ECs A3"
colnames(HeatmapEndoselect)[17] <- "ECs B1"
colnames(HeatmapEndoselect)[18] <- "ECs B2"
colnames(HeatmapEndoselect)[19] <- "ECs B3"

HeatmapEndo <- HeatmapEndoselect %>% column_to_rownames(var="symbol")
my_Endo_matrix <- as.matrix(HeatmapEndo[,c(13:18)])

Heatmap(my_Endo_matrix, column_order = c("ECs A1", "ECs A2", "ECs A3", "ECs B1", "ECs B2", "ECs B3"))


################Complex PBS Heatmap############################################################################################


HeatmapAT2PBS <-filter(N8_AT2_edgrPBS, symbol %in% c("Egf", "Egfr", "Eif2ak3", "Lamp1", "Scnn1a", "Scnn1b", "Scnn1g", "Cldn3", "Cldn4", "Cldn8", "Cldn18", "Clcn2", "Slc26a9", "Slc12a2", "Cftr", "Clcn1", "Clcn2", "Clcn3b", "Clcn4a", "Clcn4b", "Clcnka", "Clcnkb", "Bsnd", "Clca1", "Clca2", "Clca3b", "Clca4a", "Ano1", "Ano2", "Best1", "Best2", "Best3", "Slc12a1", "Slc12a4", "Slc12a5", "Slc4a8", "Slc4a10", "Kcnk2", "Mapk8", "Mapk9", "Mapk10", "Nrg1", "Erbb3", "Erbb2", "Tjp1",  "F11r", "Gja3", "Gja8", "Gja1", "Kirrel", "Afdn", "Tjp3", "Tjp2", "Ocln", "Itgb4", "Ctnnb1", "Cdh1",  "Kcnk2", "Kcnk3", "Kcnk4", "Kcnk5", "Kcnk6", "Kcnk7", "Kcnk9", "Kcnk10", "Kcnk12", "Kcnk13", "Kcnk15", "Kcnk16", "Kcnk18", "Jam2", "Jam3", "Map2k1", "Pard3", "Pard6b", "Pten"))

HeatmapAT2selectPBS <-filter(N8_AT2_edgrPBS, symbol %in% c("Cldn18", "Cldn3", "Lamp1", "Slc12a2", "Scnn1a", "Ctnnb1", "Pard3", "Afdn", "Cdh1", "Tjp2", "Tjp1", "F11r", "Pten", "Gja1", "Slc26a9", "Map2k1", "Erbb3", "Cftr", "Ocln", "Pard6b", "Eif2ak3", "Tjp3", "Mapk8", "Mapk9", "Scnn1g", "Egfr", "Cldn4", "Kcnk2", "Ano1", "Erbb2")) 

HeatmapAT2PBS <- HeatmapAT2selectPBS %>% column_to_rownames(var="symbol")
my_AT2_matrixPBS <- as.matrix(HeatmapAT2PBS[,c(13:18)])

Heatmap(my_AT2_matrixPBS, column_order = c("AT2 X1", "AT2 X2", "AT2 X3", "AT2 Y1", "AT2 Y2", "AT2 Y3"))

HeatmapEndoPBS <-filter(N8_endo_edgrPBS, symbol %in% c("Egf", "Egfr", "Eif2ak3", "Lamp1", "Scnn1a", "Scnn1b", "Scnn1g", "Cldn3", "Cldn4", "Cldn8", "Cldn18", "Clcn2", "Slc26a9", "Slc12a2", "Cftr", "Clcn1", "Clcn2", "Clcn3b", "Clcn4a", "Clcn4b", "Clcnka", "Clcnkb", "Bsnd", "Clca1", "Clca2", "Clca3b", "Clca4a", "Ano1", "Ano2", "Best1", "Best2", "Best3", "Slc12a1", "Slc12a4", "Slc12a5", "Slc4a8", "Slc4a10", "Kcnk2", "Mapk8", "Mapk9", "Mapk10", "Nrg1", "Erbb3", "Erbb2", "Tjp1",  "F11r", "Gja3", "Gja8", "Gja1", "Kirrel", "Afdn", "Tjp3", "Tjp2", "Ocln", "Itgb4", "Ctnnb1", "Cdh1",  "Kcnk2", "Kcnk3", "Kcnk4", "Kcnk5", "Kcnk6", "Kcnk7", "Kcnk9", "Kcnk10", "Kcnk12", "Kcnk13", "Kcnk15", "Kcnk16", "Kcnk18", "Jam2", "Jam3", "Map2k1", "Pard3", "Pard6b", "Pten"))

HeatmapEndoselectPBS <-filter(N8_endo_edgrPBS, symbol %in% c("Tjp1", "Lamp1", "Pard3", "Afdn", "Ctnnb1", "Pten", "Tjp2", "F11r", "Mapk9", "Jam2", "Slc12a2", "Map2k1", "Slc12a4", "Eif2ak3", "Ocln", "Mapk8", "Cldn18", "Erbb2", "Kcnk6", "Kcnk5", "Scnn1a", "Jam3", "Itgb4", "Ano2", "Gja1", "Slc4a8", "Scnn1b", "Kirrel", "Nrg1", "Pard6b")) 

colnames(HeatmapEndoselectPBS)[14] <- "ECs X1"
colnames(HeatmapEndoselectPBS)[15] <- "ECs X2"
colnames(HeatmapEndoselectPBS)[16] <- "ECs X3"
colnames(HeatmapEndoselectPBS)[17] <- "ECs Y1"
colnames(HeatmapEndoselectPBS)[18] <- "ECs Y2"
colnames(HeatmapEndoselectPBS)[19] <- "ECs Y3"

HeatmapEndoPBS <- HeatmapEndoselectPBS %>% column_to_rownames(var="symbol")
my_Endo_matrixPBS <- as.matrix(HeatmapEndoPBS[,c(13:18)])

Heatmap(my_Endo_matrixPBS, column_order = c("ECs X1", "ECs X2", "ECs X3", "ECs Y1", "ECs Y2", "ECs Y3"))

##########################################################################################################################################################################################################################################################

megaEndo <- merge(my_Endo_matrix, my_Endo_matrixPBS, by="row.names") 
megaAT2 <- merge(my_AT2_matrix, my_AT2_matrixPBS, by="row.names") 

HeatmapEndocombine <- megaEndo %>% column_to_rownames(var="Row.names")
my_Endo_matrix_combine <- as.matrix(HeatmapEndocombine)

Heatmap(my_Endo_matrix_combine, column_order = c("ECs A1", "ECs A2", "ECs A3", "ECs B1", "ECs B2", "ECs B3", "ECs X1", "ECs X2", "ECs X3", "ECs Y1", "ECs Y2", "ECs Y3"))

HeatmapAT2combine <- megaAT2 %>% column_to_rownames(var="Row.names")
my_AT2_matrix_combine <- as.matrix(HeatmapAT2combine)

Heatmap(my_AT2_matrix_combine, column_order = c("AT2 A1", "AT2 A2", "AT2 A3", "AT2 B1", "AT2 B2", "AT2 B3", "AT2 X1", "AT2 X2", "AT2 X3", "AT2 Y1", "AT2 Y2", "AT2 Y3"))





#GSEA analysis

N8_endo_edgr$FStats <- sign(N8_endo_edgr$logFC) * N8_endo_edgr$F
original_gene_list_endo <- N8_endo_edgr$FStats
names(original_gene_list_endo) <- N8_endo_edgr$symbol
head(original_gene_list_endo)
gene_list_endo <-na.omit(original_gene_list_endo)
# sort the list in decreasing order (required for clusterProfiler)
gene_list_endo = sort(gene_list_endo, decreasing = TRUE)
head(gene_list_endo)

gse_endo <- gseGO(geneList=gene_list_endo, 
                  ont ="ALL", 
                  keyType = "SYMBOL", 
                  nPerm = 10000, 
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.05, 
                  verbose = TRUE, 
                  OrgDb = org.Mm.eg.db, 
                  pAdjustMethod = "none")

gse_endo_df <- as.data.frame(gse_endo)
dotplot(gse_endo, showCategory=10, split=".sign", font.size = 2, color = "pvalue") + theme_classic() + facet_grid(.~.sign) + theme(panel.border = element_rect(fill = NA))

N8_AT2_edgr$FStats <- sign(N8_AT2_edgr$logFC) * N8_AT2_edgr$F
original_gene_list_AT2 <- N8_AT2_edgr$FStats
names(original_gene_list_AT2) <- N8_AT2_edgr$symbol
head(original_gene_list_AT2)
gene_list_AT2 <-na.omit(original_gene_list_AT2)
gene_list_AT2 = sort(gene_list_AT2, decreasing = TRUE)
head(gene_list_AT2)

gse_AT2 <- gseGO(geneList=gene_list_AT2, 
                 ont ="ALL", 
                 keyType = "SYMBOL", 
                 nPerm = 10000, 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = org.Mm.eg.db, 
                 pAdjustMethod = "none")

gse_AT2_df <- as.data.frame(gse_AT2)
dotplot(gse_AT2, showCategory=10, split=".sign", font.size = 2, color = "pvalue") + theme_classic() + facet_grid(.~.sign) + theme(panel.border = element_rect(fill = NA))

