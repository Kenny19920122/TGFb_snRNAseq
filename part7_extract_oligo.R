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

############################################################
#### Extract oligodendrocyte cells in cluster 0, 8, 11  ####
############################################################
table(Idents(itg.combined.sct))
combined.oligo <- subset(x = itg.combined.sct, idents=c("0", "8", "11"))
save(combined.oligo, file="./rdas/part7_combined.oligo.rda")


#######################################################
#### clustering analysis for oligodendrocyte cells ####
#######################################################
combined.oligo <- RunPCA(combined.oligo, verbose = FALSE)
combined.oligo <- FindNeighbors(combined.oligo, reduction = "pca", dims = 1:20)
combined.oligo <- FindClusters(combined.oligo, resolution = 0.5)
combined.oligo <- RunUMAP(combined.oligo, reduction = "pca", dims = 1:20)
combined.oligo <- RunTSNE(combined.oligo, reduction = "pca", dims = 1:20)

###################
#### tSNE plot ####
###################
tmp <- combined.oligo
Idents(object = tmp) <- "stim"
tmp <- subset(x = tmp, downsample=5000)
Idents(object = tmp) <- "seurat_clusters"
pdf(file="./figs/Oligo_15clusters_tSNE_labels.pdf")
	DimPlot(tmp, reduction = "tsne", label = TRUE, pt.size = 0.7) + NoLegend()
	DimPlot(tmp, reduction = "tsne", group.by = "stim", pt.size = 0.7)
dev.off()

################################################
#### calculate cell counts for each cluster ####
################################################
cellcounts <- matrix(nrow=15, ncol=11)
sample_ids <- c("s1001", "s1002", "s1003", "p1001", "p1002", "p1003", "p1004", "p1005", "p1006", "p1007", "p1008")
colnames(cellcounts) <- sample_ids
for(i in 1:11) {
	id <- sample_ids[i]
	tmp <- subset(x = combined.oligo, orig.ident == id)
	cellcounts[,i] <- table(Idents(tmp))
}
write.csv(cellcounts, file="./tables/Oligo_15clusters_cellcounts.csv")

#####################################
#### plot marker gene expression ####
#####################################
tmp <- Idents(object = combined.oligo)
tmp[tmp %in% c(2, 6, 8, 13)] <- 2
tmp[tmp %in% c(4, 10, 14)] <- 4
tmp[tmp %in% c(0, 5, 7)] <- 0
tmp[tmp %in% c(1, 11, 12)] <- 1
combined.oligo$newgroup <- tmp
Idents(object = combined.oligo) <- "newgroup"

genelist <- c("Pdgfra", "Cspg4", "Tnr", "Bcas1", "Tcf7l2", "Prom1", "Thbs3", "Ctps", "Opalin", "Aspa", "Apod", "Anln", "Pcsk6", "Mob3b", "S100b", "Klk6", "Scd3", "Faah", "Hopx", "Anxa5", "Mgst3", "Ptgds", "Il33", "Grm3", "Car2", "Jph4", "C4b", "Serpina3n", "Serpina3i", "Mapk11", "Il20ra", "Bdnf", "Fgfr4")
pdf(file="./figs/Oligo_15clusters_combined_marker_genes_bubbleplot.pdf", width=10, height=4)
	DotPlot(combined.oligo, scale.by = "size", features = genelist, cols = c("lightgrey", "red")) + RotatedAxis()
dev.off()
