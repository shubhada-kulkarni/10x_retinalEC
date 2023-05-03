# runnning scvelo on Denise's samples by selecting same number of cells as processed by Seurat

import scanpy as sc
import numpy as np
import pandas as pd
import scvelo as scv


adata_A = scv.read("velocyto_sampleA_newassignment/out_sorted_PGD6G.loom", cache=True)
adata_A.obs.index = [word.replace('out_sorted_PGD6G:', '') for word in adata_A.obs.index]
adata_A.var_names_make_unique()
adata_B = scv.read("velocyto_sampleB/out_merged_NZDOP.loom", cache=True)
adata_B.obs.index = [word.replace('out_merged_NZDOP:', '') for word in adata_B.obs.index]
adata_B.var_names_make_unique()

# filtering on cells and genes from Denise Seurat data
cells_A = pd.read_csv("cellnames_SampleA.csv", header = None, index_col = 0, dtype=str)
cells_B = pd.read_csv("cellnames_SampleB.csv", header = None, index_col = 0)
genes = pd.read_csv("genenames.csv", header = None, index_col = 0)

# 3 cells in sample A were not present in the scvelo object. So remove these
cells_to_remove = ['sampleA_hash6_AGCGTCGTCTCAACTT', 'sampleA_hash2_GTGCAGCCAAGCGAGT', 'sampleA_hash2_TAGTTGGTCGTTGCCT']
cells_A_updated = list(set(cells_A.index.to_list()) - set(cells_to_remove))

# checking for duplicate genes
xs = adata_A.var.index.to_list()
xs = adata_B.var.index.to_list()
import collections
[item for item, count in collections.Counter(xs).items() if count > 1]

# 13 genes from seurat were also not present in the scanpy genes list. Also remove these.
genes_to_remove = ['Gm35558.1', 'Ptp4a1.1', 'Gm16499.1', 'Gm16701.1', 'Gm5089.1', 'Snhg4.1', 'St6galnac2.1', 'Pcdha11.1', '4930594M22Rik.1', 'Fam220a.1', 'Pakap.1', 'Aldoa.1', 'Gm16364.1']
genes_updated = list(set(genes.index.to_list()) - set(genes_to_remove))

# subset
adata_A = adata_A[cells_A_updated , genes_updated]
adata_B = adata_B[cells_B.index , genes_updated]

# merge the two samples
adata = adata_A.concatenate(adata_B, batch_categories=['A', 'B'])

## Here onwards, followed the tutorial from https://scvelo.readthedocs.io/en/stable/VelocityBasics/ 
# process the merged object
#scv.pp.normalize_per_cell(adata)
sc.pp.normalize_total(adata, target_sum=1e4)
scv.pp.log1p(adata)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.umap(adata)
scv.tl.louvain(adata, resolution = 2)
scv.pl.umap(adata, color='louvain')   # plotting the umap with clusters
scv.pl.umap(adata, color=['louvain'])

# Estimating the RNA velocity
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

# Project the velocities
scv.pl.velocity_embedding_stream(adata, basis='umap')
scv.pl.velocity_embedding(adata, arrow_length=3, arrow_size=2, dpi=120)

# add the cre and the cluster information from Seurat to the anndata object
adata.obs['sample_hash'] = ['_'.join(i.split('_')[:2]) for i in adata.obs_names.to_list()]
hash_dict = {"sampleA_hash1": "neg", "sampleA_hash3": "neg", "sampleA_hash4": "neg", "sampleA_hash6":"neg", "sampleA_hash9":"neg",
"sampleB_hash3": "neg", "sampleB_hash4": "neg", "sampleB_hash6":"neg", "sampleB_hash8" :"neg", "sampleA_hash2":"pos", "sampleA_hash5":"pos", "sampleA_hash7": "pos", "sampleA_hash8": "pos", "sampleA_hash10":"pos", "sampleB_hash1":"pos","sampleB_hash2":"pos", "sampleB_hash5":"pos", "sampleB_hash7":"pos", "sampleB_hash9":"pos", "sampleB_hash10":"pos"}
adata.obs['cre'] = (adata.obs['sample_hash'].map(hash_dict).astype('category'))

# UMAP with all informations
scv.pl.umap(adata, color=['louvain', 'sample_hash', 'cre']) # 'new_clusters' add these later

# Finding marker genes
sc.tl.rank_genes_groups(adata, 'louvain', method='t-test')
sc.pl.rank_genes_groups(adata)

## saving marker genes into a file
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
df = pd.DataFrame({group + '_' + key[:1]: result[key][group]
for group in groups for key in ['names','logfoldchanges','pvals','pvals_adj']})
df.to_csv('markergenes_scanpy.csv')

# checking marker expressions in clusters
scv.pl.umap(adata, color=['louvain', "Col4a1", "Apol7d", "Cd93", "Nid1"])
