import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd

# load sparse matrix:
X = io.mmread("/home/users/yliu/projects/scRNA/Kenny_v2/tables/combined.oligo_15cluster_counts.mtx")

# create anndata object
adata = anndata.AnnData(
    X=X.transpose().tocsr()
)

# load cell metadata:
cell_meta = pd.read_csv("/home/users/yliu/projects/scRNA/Kenny_v2/tables/combined.oligo_15cluster_metadata.csv")

# load gene names:
with open("/home/users/yliu/projects/scRNA/Kenny_v2/tables/combined.oligo_15cluster_gene_names.csv", 'r') as f:
    gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

# load dimensional reduction:
pca = pd.read_csv("/home/users/yliu/projects/scRNA/Kenny_v2/tables/combined.oligo_15cluster_pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T
adata.obsm['X_tsne'] = np.vstack((adata.obs['tSNE_1'].to_numpy(), adata.obs['tSNE_2'].to_numpy())).T

# save dataset as anndata format
adata.write('/home/users/yliu/projects/scRNA/Kenny_v2/tables/combined.oligo_15cluster_my_data.h5ad')
