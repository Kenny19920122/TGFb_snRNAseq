setwd("/home/users/yliu/projects/scRNA/Kenny_v2")

library(FastCAR)
library(Matrix)
library(Seurat)
library(qlcMatrix)
library(pheatmap)
library(ggplot2)
library(gridExtra)
library(stringr)
library(DoubletFinder)

contaminationChanceCutoff = 0.05


# Define the sample IDs
sample_ids <- c("1001", "1002", "1003")

for (id in sample_ids) {
  print(paste0("s", id))	
  # Read the data
  data_dir <- paste0("./cellranger/count/P26770_", id, "/run_count_P26770_", id, "/outs/")
  data <- Read10X(data.dir = paste0(data_dir, "filtered_feature_bc_matrix/"))
  full <- Read10X(data.dir = paste0(data_dir, "raw_feature_bc_matrix/"))
  
  # Describe ambient RNA sequence
  ambProfile <- describe.ambient.RNA.sequence(fullCellMatrix = full, 
		start = 10, stop = 500, by = 10, contaminationChanceCutoff = contaminationChanceCutoff)
  
  # Save the ambient profile plot
  pdf(paste0("./figs/revision_ambient_s", id, ".pdf"))
  plot.ambient.profile(ambProfile)
  dev.off()
  
  # Determine background to remove
  emptyDropletCutoff <- recommend.empty.cutoff(ambProfile)
  print(emptyDropletCutoff)
  ambientProfile <- determine.background.to.remove(full, emptyDropletCutoff, contaminationChanceCutoff)
  data <- remove.background(data, ambientProfile)
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = data, project = paste0("s", id), min.cells = 3, min.features = 500)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 800 & nFeature_RNA < 3000 & percent.mt < 5)
  
  # Assign the Seurat object to a variable
  assign(paste0("s", id), seurat_obj)
}


# Define the sample IDs
sample_ids <- c("1001", "1002", "1003", "1004", "1005", "1006", "1007", "1008")

for (id in sample_ids) {
  print(paste0("p", id))	
  # Read the data
  data_dir <- paste0("./cellranger/count/P32113_", id, "/run_count_P32113_", id, "/outs/")
  data <- Read10X(data.dir = paste0(data_dir, "filtered_feature_bc_matrix/"))
  full <- Read10X(data.dir = paste0(data_dir, "raw_feature_bc_matrix/"))
  
  # Describe ambient RNA sequence
  ambProfile <- describe.ambient.RNA.sequence(fullCellMatrix = full, 
		start = 10, stop = 500, by = 10, contaminationChanceCutoff = contaminationChanceCutoff)
  
  # Save the ambient profile plot
  pdf(paste0("./figs/revision_ambient_p", id, ".pdf"))
  plot.ambient.profile(ambProfile)
  dev.off()
  
  # Determine background to remove
  emptyDropletCutoff <- recommend.empty.cutoff(ambProfile)
  print(emptyDropletCutoff)
  ambientProfile <- determine.background.to.remove(full, emptyDropletCutoff, contaminationChanceCutoff)
  data <- remove.background(data, ambientProfile)
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = data, project = paste0("p", id), min.cells = 3, min.features = 500)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 800 & nFeature_RNA < 3000 & percent.mt < 5)
  
  # Assign the Seurat object to a variable
  assign(paste0("p", id), seurat_obj)
}

s1001$stim <- p1001$stim <- p1002$stim <- "ctrl"
p1003$stim <- p1004$stim <- "d10"
s1002$stim <- p1005$stim <- p1006$stim <- "d20"
s1003$stim <- p1007$stim <- p1008$stim <- "d30"


sample_ids <- c("s1001", "s1002", "s1003", "p1001", "p1002", "p1003", "p1004", "p1005", "p1006", "p1007", "p1008")

for (id in sample_ids) {
  # SCTransform, PCA and UMAP
  seurat_obj <- get(id)
  seurat_obj <- SCTransform(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:50)
  
  # ParamSweep and summarizeSweep
  sweep_res_list <- paramSweep(seurat_obj, PCs = 1:50, sct = TRUE, num.cores = 10)
  sweep_stats <- summarizeSweep(sweep_res_list, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  
  # DoubletFinder
  nExp_poi <- round(0.075 * nrow(seurat_obj@meta.data))
  print(paste(id, pK_bcmvn, nExp_poi, sep=" "))
  seurat_obj <- doubletFinder(seurat_obj, PCs = 1:50, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
    
  # Assign the Seurat object to a variable
  assign(id, seurat_obj)
}


s1001 <- subset(s1001, subset = DF.classifications_0.25_0.28_435 == "Singlet")
s1002 <- subset(s1002, subset = DF.classifications_0.25_0.06_951 == "Singlet")
s1003 <- subset(s1003, subset = DF.classifications_0.25_0.17_115 == "Singlet")
p1001 <- subset(p1001, subset = DF.classifications_0.25_0.1_1020 == "Singlet")
p1002 <- subset(p1002, subset = DF.classifications_0.25_0.28_370 == "Singlet")
p1003 <- subset(p1003, subset = DF.classifications_0.25_0.3_1509 == "Singlet")
p1004 <- subset(p1004, subset = DF.classifications_0.25_0.15_453 == "Singlet")
p1005 <- subset(p1005, subset = DF.classifications_0.25_0.11_3483 == "Singlet")
p1006 <- subset(p1006, subset = DF.classifications_0.25_0.11_2648 == "Singlet")
p1007 <- subset(p1007, subset = DF.classifications_0.25_0.09_2370 == "Singlet")
p1008 <- subset(p1008, subset = DF.classifications_0.25_0.29_1406 == "Singlet")

save(s1001,s1002,s1003,p1001,p1002,p1003,p1004,p1005,p1006,p1007,p1008, file="./rdas/part1.rda")

