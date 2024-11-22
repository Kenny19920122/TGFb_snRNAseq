import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2

# reload dataset
adata = sc.read_h5ad('/home/users/yliu/projects/scRNA/Kenny_v2/tables/combined.oligo_15cluster_my_data.h5ad')

# load loom files for spliced/unspliced matrices for each sample:
ldata1 = scv.read('/home/users/yliu/projects/scRNA/Kenny_v2/cellranger/count/P26770_1001/run_count_P26770_1001/velocyto/run_count_P26770_1001.loom', cache=True)
ldata2 = scv.read('/home/users/yliu/projects/scRNA/Kenny_v2/cellranger/count/P26770_1002/run_count_P26770_1002/velocyto/run_count_P26770_1002.loom', cache=True)
ldata3 = scv.read('/home/users/yliu/projects/scRNA/Kenny_v2/cellranger/count/P26770_1003/run_count_P26770_1003/velocyto/run_count_P26770_1003.loom', cache=True)
ldata4 = scv.read('/home/users/yliu/projects/scRNA/Kenny_v2/cellranger/count/P32113_1001/run_count_P32113_1001/velocyto/run_count_P32113_1001.loom', cache=True)
ldata5 = scv.read('/home/users/yliu/projects/scRNA/Kenny_v2/cellranger/count/P32113_1002/run_count_P32113_1002/velocyto/run_count_P32113_1002.loom', cache=True)
ldata6 = scv.read('/home/users/yliu/projects/scRNA/Kenny_v2/cellranger/count/P32113_1003/run_count_P32113_1003/velocyto/run_count_P32113_1003.loom', cache=True)
ldata7 = scv.read('/home/users/yliu/projects/scRNA/Kenny_v2/cellranger/count/P32113_1004/run_count_P32113_1004/velocyto/run_count_P32113_1004.loom', cache=True)
ldata8 = scv.read('/home/users/yliu/projects/scRNA/Kenny_v2/cellranger/count/P32113_1005/run_count_P32113_1005/velocyto/run_count_P32113_1005.loom', cache=True)
ldata9 = scv.read('/home/users/yliu/projects/scRNA/Kenny_v2/cellranger/count/P32113_1006/run_count_P32113_1006/velocyto/run_count_P32113_1006.loom', cache=True)
ldata10 = scv.read('/home/users/yliu/projects/scRNA/Kenny_v2/cellranger/count/P32113_1007/run_count_P32113_1007/velocyto/run_count_P32113_1007.loom', cache=True)
ldata11 = scv.read('/home/users/yliu/projects/scRNA/Kenny_v2/cellranger/count/P32113_1008/run_count_P32113_1008/velocyto/run_count_P32113_1008.loom', cache=True)


# rename barcodes in order to merge:
barcodes = [bc.split(':')[1] for bc in ldata1.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_1' for bc in barcodes]
ldata1.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata2.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_2' for bc in barcodes]
ldata2.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata3.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_3' for bc in barcodes]
ldata3.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata4.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_4' for bc in barcodes]
ldata4.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata5.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_5' for bc in barcodes]
ldata5.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata6.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_6' for bc in barcodes]
ldata6.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata7.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_7' for bc in barcodes]
ldata7.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata8.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_8' for bc in barcodes]
ldata8.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata9.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_9' for bc in barcodes]
ldata9.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata10.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_10' for bc in barcodes]
ldata10.obs.index = barcodes

barcodes = [bc.split(':')[1] for bc in ldata11.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '-1_11' for bc in barcodes]
ldata11.obs.index = barcodes


ldata1 = ldata1[np.isin(ldata1.obs.index, adata.obs.index)]
ldata2 = ldata2[np.isin(ldata2.obs.index, adata.obs.index)]
ldata3 = ldata3[np.isin(ldata3.obs.index, adata.obs.index)]
ldata4 = ldata4[np.isin(ldata4.obs.index, adata.obs.index)]
ldata5 = ldata5[np.isin(ldata5.obs.index, adata.obs.index)]
ldata6 = ldata6[np.isin(ldata6.obs.index, adata.obs.index)]
ldata7 = ldata7[np.isin(ldata7.obs.index, adata.obs.index)]
ldata8 = ldata8[np.isin(ldata8.obs.index, adata.obs.index)]
ldata9 = ldata9[np.isin(ldata9.obs.index, adata.obs.index)]
ldata10 = ldata10[np.isin(ldata10.obs.index, adata.obs.index)]
ldata11 = ldata11[np.isin(ldata11.obs.index, adata.obs.index)]

# make variable names unique
ldata1.var_names_make_unique()
ldata2.var_names_make_unique()
ldata3.var_names_make_unique()
ldata4.var_names_make_unique()
ldata5.var_names_make_unique()
ldata6.var_names_make_unique()
ldata7.var_names_make_unique()
ldata8.var_names_make_unique()
ldata9.var_names_make_unique()
ldata10.var_names_make_unique()
ldata11.var_names_make_unique()

# concatenate the three loom
ldata = ldata1.concatenate([ldata2, ldata3, ldata4, ldata5, ldata6, ldata7, ldata8, ldata9, ldata10, ldata11])

# merge matrices into the original adata object
adata = scv.utils.merge(adata, ldata)
scv.pl.proportions(adata, groupby='celltype_full', save='_oligo.pdf')

# pre-process
scv.pp.filter_and_normalize(adata)

sc.pp.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)

# compute velocity
scv.tl.recover_dynamics(adata, n_jobs = 20)
scv.tl.velocity(adata, mode = 'dynamical')
scv.tl.velocity_graph(adata, n_jobs = 20)
scv.pl.velocity_embedding(adata, basis='tsne', color='newgroup', frameon=False, save='oligo_15cluster_combined_dynamical_embedding_tsne.pdf')
scv.pl.velocity_embedding_grid(adata, basis='tsne', color='newgroup', save='oligo_15cluster_combined_dynamical_embedding_grid_tsne.pdf', title='', scale=0.25)
scv.pl.velocity_embedding_stream(adata, basis='tsne', color='newgroup', save='oligo_15cluster_combined_dynamical_embedding_stream_tsne.pdf', title='')

scv.tl.latent_time(adata)
scv.pl.scatter(adata, basis='tsne', color = 'latent_time', save = 'oligo_15cluster_combined_dynamical_latent_time_tsne.pdf')