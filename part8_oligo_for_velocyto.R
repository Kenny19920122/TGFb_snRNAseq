setwd("/home/users/yliu/projects/scRNA/Kenny_v2")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

load(file="./rdas/part7_combined.oligo.rda")


#######################################################
#### clustering analysis for oligodendrocyte cells ####
#######################################################
combined.oligo <- RunPCA(combined.oligo, verbose = FALSE)
combined.oligo <- FindNeighbors(combined.oligo, reduction = "pca", dims = 1:20)
combined.oligo <- FindClusters(combined.oligo, resolution = 0.5)
combined.oligo <- RunUMAP(combined.oligo, reduction = "pca", dims = 1:20)
combined.oligo <- RunTSNE(combined.oligo, reduction = "pca", dims = 1:20)


tmp <- Idents(object = combined.oligo)
tmp[tmp %in% c(2, 6, 8, 13)] <- 2
tmp[tmp %in% c(4, 10, 14)] <- 4
tmp[tmp %in% c(0, 5, 7)] <- 0
tmp[tmp %in% c(1, 11, 12)] <- 1
tmp[tmp == 9] <- 5
tmp <- droplevels(tmp)
combined.oligo$newgroup <- tmp
Idents(object = combined.oligo) <- "newgroup"


combined.oligo$barcode <- colnames(combined.oligo)
combined.oligo$UMAP_1 <- combined.oligo@reductions$umap@cell.embeddings[,1]
combined.oligo$UMAP_2 <- combined.oligo@reductions$umap@cell.embeddings[,2]
combined.oligo$tSNE_1 <- combined.oligo@reductions$tsne@cell.embeddings[,1]
combined.oligo$tSNE_2 <- combined.oligo@reductions$tsne@cell.embeddings[,2]
meta.data <- combined.oligo@meta.data
write.csv(meta.data, file='./tables/combined.oligo_15cluster_metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- combined.oligo$SCT@counts
writeMM(counts_matrix, file='./tables/combined.oligo_15cluster_counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(combined.oligo@reductions$pca@cell.embeddings, file='./tables/combined.oligo_15cluster_pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='./tables/combined.oligo_15cluster_gene_names.csv',
  quote=F,row.names=F,col.names=F
)
