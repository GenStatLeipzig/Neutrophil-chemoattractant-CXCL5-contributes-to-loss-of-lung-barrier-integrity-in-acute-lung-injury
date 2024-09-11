require(toolboxH)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Rfast2)
library(ggrepel)
library(cowplot)
library(patchwork)
require(ggthemes) 
library(reshape2)
library(ggforce)

StSeuMerged <- readRDS("Data/ST_Seu_Cxcl5_final.rds")

StSeuMerged@active.ident = factor(StSeuMerged$SCT_snn_res.0.4)
colors04 <- c("0" = "#2176ae", "1" = "#FC7536", "2" = "#fbb13c", "3" = "#b66d0d" ,"4" = "#87a878" ,"5" = "#dbf9b8", "6" =  "#EC3241","7" =  "#57b8ff","8" =  "#00998F","9" =  "#B8C0FF","10" =  "#c7ccb9")
SpatialDimPlot(StSeuMerged, crop = F, label.size = 3, cols = colors04)

#heatmap
StSeuMerged.markers04 <- FindAllMarkers(StSeuMerged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
StSeuMerged.markers04 %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- StSeuMerged.markers04 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 <- StSeuMerged.markers04 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(StSeuMerged, features = top5$gene, group.colors = c(colors04)) + NoLegend()
DoHeatmap(StSeuMerged, features = top10$gene, group.colors = c(colors04)) + NoLegend()

#diff expr
StSeuMerged$cluster_group04 <- paste(StSeuMerged$group, StSeuMerged$SCT_snn_res.0.4, sep = "_")
StSeuMerged@active.ident <- factor(StSeuMerged$cluster_group04)
Idents(StSeuMerged) <- "cluster_group04"

#C1
diff.markersC1 <- FindMarkers(StSeuMerged, ident.1 = "LixkoSpn_1", ident.2 = "WtSpn_1") 
diff.markersC1 <- tibble::rownames_to_column(diff.markersC1, "symbol")

diff.markersC1$Significant <- "FALSE"
diff.markersC1$Significant[diff.markersC1$p_val_adj < 0.05] <- "TRUE"
diff.markersC1$diffexpressed <- "FALSE"
diff.markersC1$diffexpressed[diff.markersC1$avg_log2FC > 0 & diff.markersC1$p_val_adj < 0.05] <- "UP"
diff.markersC1$diffexpressed[diff.markersC1$avg_log2FC < 0 & diff.markersC1$p_val_adj < 0.05] <- "DOWN"

ggplot(data=diff.markersC1, aes(x=avg_log2FC, y=-log10(p_val_adj),col=diffexpressed, label=symbol)) + geom_point(size = 4) + theme_classic() + geom_text_repel(size = 8, max.overlaps = 3) + geom_vline(xintercept=c(-1, 1), linetype = "dashed") + geom_hline(yintercept=-log10(0.05), linetype = "dashed") + scale_colour_manual(values=c("blue", "black", "red")) + NoLegend()

#C2
diff.markersC2 <- FindMarkers(StSeuMerged, ident.1 = "LixkoSpn_2", ident.2 = "WtSpn_2") 
diff.markersC2 <- tibble::rownames_to_column(diff.markersC2, "symbol")

diff.markersC2$Significant <- "FALSE"
diff.markersC2$Significant[diff.markersC2$p_val_adj < 0.05] <- "TRUE"

diff.markersC2$diffexpressed <- "FALSE"
diff.markersC2$diffexpressed[diff.markersC2$avg_log2FC > 0 & diff.markersC2$p_val_adj < 0.05] <- "UP"
diff.markersC2$diffexpressed[diff.markersC2$avg_log2FC < 0 & diff.markersC2$p_val_adj < 0.05] <- "DOWN"
ggplot(data=diff.markersC2, aes(x=avg_log2FC, y=-log10(p_val_adj),col=diffexpressed, label=symbol)) + geom_point(size = 4) + theme_classic() + geom_text_repel(size = 8, max.overlaps = 3) + geom_vline(xintercept=c(-1, 1), linetype = "dashed") + geom_hline(yintercept=-log10(0.05), linetype = "dashed") + scale_colour_manual(values=c("blue", "black", "red")) + NoLegend()



#C3
diff.markersC3 <- FindMarkers(StSeuMerged, ident.1 = "LixkoSpn_3", ident.2 = "WtSpn_3") 
diff.markersC3 <- tibble::rownames_to_column(diff.markersC3, "symbol")

diff.markersC3$Significant <- "FALSE"
diff.markersC3$Significant[diff.markersC3$p_val_adj < 0.05] <- "TRUE"

diff.markersC3$diffexpressed <- "FALSE"
diff.markersC3$diffexpressed[diff.markersC3$avg_log2FC > 0 & diff.markersC3$p_val_adj < 0.05] <- "UP"
diff.markersC3$diffexpressed[diff.markersC3$avg_log2FC < 0 & diff.markersC3$p_val_adj < 0.05] <- "DOWN"
ggplot(data=diff.markersC3, aes(x=avg_log2FC, y=-log10(p_val_adj),col=diffexpressed, label=symbol)) + geom_point(size = 4) + theme_classic() + geom_text_repel(size = 8, max.overlaps = 3) + geom_vline(xintercept=c(-1, 1), linetype = "dashed") + geom_hline(yintercept=-log10(0.05), linetype = "dashed") + scale_colour_manual(values=c("blue", "black", "red")) + NoLegend()

Idents(StSeuMerged) <- "group"
StSeuMergedNoctr <- subset(StSeuMerged, idents = c("WtSpn", "LixkoSpn"))
Idents(StSeuMergedNoctr) <- "SCT_snn_res.0.4"
Cluster2 <- subset(StSeuMergedNoctr, idents = "2")
Cluster3 <- subset(StSeuMergedNoctr, idents = "3")
Cluster1 <- subset(StSeuMergedNoctr, idents = "1")

Idents(Cluster3) <- "group"
VlnPlot(Cluster3, features = "Il6", cols = c("black", "#285D5E")) + NoLegend()

Idents(Cluster2) <- "group"
VlnPlot(Cluster2, features = "Il6", cols = c("black", "#285D5E")) + NoLegend()

Idents(Cluster3) <- "group"
VlnPlot(Cluster3, features = "Cxcl1", cols = c("black", "#285D5E")) + NoLegend()

Idents(Cluster2) <- "group"
VlnPlot(Cluster2, features = "Cxcl1", cols = c("black", "#285D5E")) + NoLegend()

Idents(Cluster3) <- "group"
VlnPlot(Cluster3, features = "Cxcl5", cols = c("black", "#285D5E")) + NoLegend()

Idents(Cluster2) <- "group"
VlnPlot(Cluster2, features = "Cxcl5", cols = c("black", "#285D5E")) + NoLegend()

Idents(Cluster3) <- "group"
VlnPlot(Cluster3, features = "Cxcl10", cols = c("black", "#285D5E")) + NoLegend()

Idents(Cluster2) <- "group"
VlnPlot(Cluster2, features = "Cxcl10", cols = c("black", "#285D5E")) + NoLegend()

SpatialFeaturePlot(StSeuMerged, features = "Cxcl1", crop = F, alpha = c(0.1, 1))
SpatialFeaturePlot(StSeuMerged, features = "Cxcl5", crop = F, alpha = c(0.1, 1))
SpatialFeaturePlot(StSeuMerged, features = "Cxcl10", crop = F, alpha = c(0.1, 1))
SpatialFeaturePlot(StSeuMerged, features = "Il6", crop = F, alpha = c(0.1, 1))
