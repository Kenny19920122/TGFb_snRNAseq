setwd("/home/users/yliu/projects/scRNA/Kenny_v2")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

load(file="./rdas/part3_combined.mg.rda")
combined.mg <- RunPCA(combined.mg, verbose = FALSE)
combined.mg <- FindNeighbors(combined.mg, reduction = "pca", dims = 1:10)
combined.mg <- FindClusters(combined.mg, resolution = 0.15)
combined.mg <- RunUMAP(combined.mg, reduction = "pca", dims = 1:10)
combined.mg <- RunTSNE(combined.mg, reduction = "pca", dims = 1:10)

combined.mg$barcode <- colnames(combined.mg)
combined.mg$UMAP_1 <- combined.mg@reductions$umap@cell.embeddings[,1]
combined.mg$UMAP_2 <- combined.mg@reductions$umap@cell.embeddings[,2]
combined.mg$tSNE_1 <- combined.mg@reductions$tsne@cell.embeddings[,1]
combined.mg$tSNE_2 <- combined.mg@reductions$tsne@cell.embeddings[,2]
write.csv(combined.mg@meta.data, file='./tables/combined.mg_metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- combined.mg$SCT@counts
writeMM(counts_matrix, file='./tables/combined.mg_counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(combined.mg@reductions$pca@cell.embeddings, file='./tables/combined.mg_pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='./tables/combined.mg_gene_names.csv',
  quote=F,row.names=F,col.names=F
)
