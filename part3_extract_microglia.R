setwd("/home/users/yliu/projects/scRNA/Kenny_v2")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(pheatmap)
library(RColorBrewer)

load("./rdas/part2_itg.combined.sct.rda")

##########################################
#### perform clustering (21 clusters) ####
##########################################
itg.combined.sct <- FindNeighbors(itg.combined.sct, reduction = "integrated.dr", dims = 1:50)
itg.combined.sct <- FindClusters(itg.combined.sct, resolution = 0.2)
itg.combined.sct <- RunUMAP(itg.combined.sct, dims = 1:50, reduction = "integrated.dr")

###############################################
#### Extract microglia cells in cluster 4  ####
###############################################
table(Idents(itg.combined.sct))
combined.mg <- subset(x = itg.combined.sct, idents="4")
save(combined.mg, file="./rdas/part3_combined.mg.rda")

##########################################
#### clustering analysis for MG cells ####
##########################################
combined.mg <- RunPCA(combined.mg, verbose = FALSE)
combined.mg <- FindNeighbors(combined.mg, reduction = "pca", dims = 1:10)
combined.mg <- FindClusters(combined.mg, resolution = 0.15)
combined.mg <- RunUMAP(combined.mg, reduction = "pca", dims = 1:10)
combined.mg <- RunTSNE(combined.mg, reduction = "pca", dims = 1:10)

###################
#### UMAP plot ####
###################
tmp <- combined.mg
Idents(object = tmp) <- "stim"
tmp <- subset(x = tmp, downsample=500)
Idents(object = tmp) <- "seurat_clusters"
pdf(file="./figs/Microglia_4clusters_UMAP.pdf")
	DimPlot(tmp, reduction = "umap", label = TRUE, pt.size = 1.2) + NoLegend()
dev.off()

################################################
#### calculate cell counts for each cluster ####
################################################
cellcounts <- matrix(nrow=4, ncol=11)
sample_ids <- c("s1001", "s1002", "s1003", "p1001", "p1002", "p1003", "p1004", "p1005", "p1006", "p1007", "p1008")
colnames(cellcounts) <- sample_ids
for(i in 1:11) {
	id <- sample_ids[i]
	tmp <- subset(x = combined.mg, orig.ident == id)
	cellcounts[,i] <- table(Idents(tmp))
}
write.csv(cellcounts, file="./tables/Microglia_4cluster_cellcounts.csv")


#####################################
#### plot marker gene expression ####
#####################################
genelist <- c("Mef2a", "Tgfbr1", "Hexb", "Cx3cr1", "P2ry12", "Sall1", "Siglech", "Ctnnd2", "Gpc5", "Gpc4", "Sparcl1", "Sdc2", "Cdh22", "Prrx1", "Parva",  "Pdgfd", "Bcan", "Wwc1", "Egfr", "Frem1", "Adamts9", "Adamts5", "Yap1", "Tgfb2", "Tnc", "Apoe", "Lpl", "Clec7a", "Tyrobp", "Trem2", "Itgax", "Lgals3", "Gpnmb", "Mgll", "Alcam", "Airn", "Pparg", "Mmp12", "Tnf", "Cd74", "B2m", "Ciita", "H2-Aa", "Lrp1", "Cd36", "Nr1h3", "Abca1", "Abcg1", "Soat1", "Plin2", "C4b", "Csf1", "Cd68", "Mbp", "Plp1", "Mobp", "Mag")
pdf(file="./figs/Microglia_4clusters_marker_genes_bubbleplot.pdf", width=20, height=4)
	DotPlot(combined.mg, scale.by = "size", features = genelist, cols = c("lightgrey", "red")) + RotatedAxis()
dev.off()

markers.to.plot <- c("Lgals3", "Gpnmb", "Mgll", "Alcam", "Airn", "Pparg")
pdf(file="./figs/Microglia_4cluster_marker_genes.pdf", width=5, height=50)
	plots <- VlnPlot(combined.mg, features = markers.to.plot, split.by =  "seurat_clusters", pt.size = 0, combine = FALSE, log=T, adjust = 0.8)
	wrap_plots(plots = plots, ncol = 1)
dev.off()

###############################################
#### identify marker gene for each cluster ####
###############################################
all.markers <- FindAllMarkers(object = combined.mg, recorrect_umi = F)
write.csv(all.markers,file="./tables/Microglia_4clusters_allmarker_genes.csv")
