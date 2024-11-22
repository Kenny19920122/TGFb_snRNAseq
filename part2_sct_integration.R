setwd("/home/users/yliu/projects/scRNA/Kenny_v2")

library(Seurat)
library(SeuratObject)

load("./rdas/part1.rda")

# integrate datasets
itg.combined.sct <- merge(x = s1001, y = list(s1002,s1003,p1001,p1002,p1003,p1004,p1005,p1006,p1007,p1008))

library(future)
options(future.globals.maxSize = 100000 * 1024^5)
itg.combined.sct <- SCTransform(itg.combined.sct, vars.to.regress = "percent.mt")

itg.combined.sct <- RunPCA(itg.combined.sct)
itg.combined.sct <- RunUMAP(itg.combined.sct, dims = 1:50)
itg.combined.sct <- IntegrateLayers(object = itg.combined.sct, method = CCAIntegration, normalization.method = "SCT", verbose = F)
save(itg.combined.sct,file="./rdas/part2_itg.combined.sct.rda")

##########################################
#### perform clustering (21 clusters) ####
##########################################
itg.combined.sct <- FindNeighbors(itg.combined.sct, reduction = "integrated.dr", dims = 1:50)
itg.combined.sct <- FindClusters(itg.combined.sct, resolution = 0.2)
itg.combined.sct <- RunUMAP(itg.combined.sct, dims = 1:50, reduction = "integrated.dr")

###################
#### UMAP plot ####
###################
tmp <- itg.combined.sct
Idents(object = tmp) <- "stim"
tmp <- subset(x = tmp, downsample=7500)
Idents(object = tmp) <- "seurat_clusters"
pdf(file="./figs/SCT_UMAP_cluster21_labels.pdf")
	DimPlot(tmp, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

################################################
#### calculate cell counts for each cluster ####
################################################
cellcounts <- matrix(nrow=21, ncol=11)
sample_ids <- c("s1001", "s1002", "s1003", "p1001", "p1002", "p1003", "p1004", "p1005", "p1006", "p1007", "p1008")
colnames(cellcounts) <- sample_ids
for(i in 1:11) {
	id <- sample_ids[i]
	tmp <- subset(x = itg.combined.sct, orig.ident == id)
	cellcounts[,i] <- table(Idents(tmp))
}
write.csv(cellcounts, file="./tables/sct_cluster21_cellcounts.csv")


#####################################
#### plot marker gene expression ####
#####################################
genelist <- c("Aqp4", "Gfap", "Col23a1", "Hsbp1", "Ptprc", "Itgam", "Hexb", "Sall1", "Mertk", "Tgfbr1", "Csf1r", "Pdgfra", "Cspg4", "Olig2", "Ptprz1", "Tnr", "Tcf7l2", "Bcas1", "Fyn", "Mbp", "Plp1", "St18", "Prr5l", "Mobp", "Scd3", "Bsg", "Pecam1", "Flt1", "Vtn", "Abcc9", "Rgs5", "Dnah12", "Nnat", "Cfap43", "Lama1", "Cemip", "Cped1", "Pkd1l2", "Pkd2l1", "Myo3b", "Rbfox3", "Syt1", "Snap25", "Slc17a6", "Slc32a1", "Gad1", "Gad2", "Chat", "Col6a6")
pdf(file="./figs/all_21clusters_marker_genes_bubbleplot.pdf", width=20, height=7)
	DotPlot(itg.combined.sct, scale.by = "size", features = genelist, cols = c("lightgrey", "red")) + RotatedAxis()
dev.off()
