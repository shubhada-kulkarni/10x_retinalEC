# file with all commands used to create dynamo for Denies's 10x data

# import required modules
import dynamo as dyn
import anndata   # reading with anndata didnt work so tried scanpy
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# from gseapy.plot import barplot, dotplot
import warnings
warnings.filterwarnings('ignore')


# reading data with scanpy
adata_A = sc.read_10x_mtx("output_SampleA/outs/multi/count/raw_feature_bc_matrix/")
adata_B = sc.read_10x_mtx("output_SampleB/outs/multi/count/raw_feature_bc_matrix/")

# read cells that needs to be filtered. And replace the sample names in the cellnames as original data does not have these
# also read subset genes
cells_A = [line.rstrip() for line in open('cellnames_SampleA.csv')]
cells_B = [line.rstrip() for line in open('cellnames_SampleB.csv')]
genes = [line.rstrip() for line in open('genenames.csv')]

# subset the above adata only for these cells
adata_A = adata_A[cells_A , genes]
adata_B = adata_B[cells_B , genes]

########################################################################################
#### seurat to anndata from Theis lab (also see the code in analysis_Seurat.Rmd file)
# read data into pandas dataframe
import pandas as pd
exprs = pd.read_csv('exprs_denise_10x_annotated.csv', index_col=0)
meta = pd.read_csv('meta_denise_10x_annotated.csv', index_col=0)
feature_meta = pd.read_csv('feature_meta_denise_10x_annotated.csv', index_col=0)
embedding = pd.read_csv('embedding_denise_10x_annotated.csv', index_col=0)

import numpy as np
exprs = np.genfromtxt('exprs_denise_10x_annotated.csv', delimiter=",", names = True)

# create anndata object
adata = sc.AnnData(X = exprs)
adata.obs_names = meta
adata.var_names = feature_meta
adata.obsm['umap'] = embedding
sc.pl.umap(adata, color='louvain')

#########################################################################################
# creating the unspliced and spliced labels using velocyto before running dynamo.
# in terminal
# velocyto run -b cellbarcodes_SampleA.txt -o velocyto_sampleA output_SampleA/outs/multi/count/unassigned_alignments.bam /omics/odcf/analysis/OE0228_projects/endothelial_cells/s467i/scRNA-seq/try_cellranger/refdata-gex-mm10-2020-A/genes/genes.gtf
# velocyto run -b cellbarcodes_SampleB.txt -o velocyto_sampleB output_SampleB/outs/multi/count/unassigned_alignments.bam /omics/odcf/analysis/OE0228_projects/endothelial_cells/s467i/scRNA-seq/try_cellranger/refdata-gex-mm10-2020-A/genes/genes.gtf

# 09-02-2023
# because above velocyto command was run for unassigned cells, reran it for assigned bam files. For this, first merged the bamfiles per sample using mergeBams
mergeBams -i DG22-1/count/sample_alignments.bam,DG22-2/count/sample_alignments.bam,DG22-3/count/sample_alignments.bam,DG22-4/count/sample_alignments.bam,DG22-5/count/sample_alignments.bam,DG22-6/count/sample_alignments.bam,DG22-7/count/sample_alignments.bam,DG22-8/count/sample_alignments.bam,DG22-9/count/sample_alignments.bam,DG22-10/count/sample_alignments.bam -l sampleA_hash1_,sampleA_hash2_,sampleA_hash3_,sampleA_hash4_,sampleA_hash5_,sampleA_hash6_,sampleA_hash7_,sampleA_hash8_,sampleA_hash9_,sampleA_hash10_ -o merged_bam
mergeBams -i DG22-11/count/sample_alignments.bam,DG22-12/count/sample_alignments.bam,DG22-13/count/sample_alignments.bam,DG22-14/count/sample_alignments.bam,DG22-15/count/sample_alignments.bam,DG22-16/count/sample_alignments.bam,DG22-17/count/sample_alignments.bam,DG22-18/count/sample_alignments.bam,DG22-19/count/sample_alignments.bam,DG22-20/count/sample_alignments.bam -l sampleB_hash1_,sampleB_hash2_,sampleB_hash3_,sampleB_hash4_,sampleB_hash5_,sampleB_hash6_,sampleB_hash7_,sampleB_hash8_,sampleB_hash9_,sampleB_hash10_ -o merged_bam_sampleB

# then sort the bam files
samtools sort output_SampleA/outs/per_sample_outs/merged_bam/out.bam -o output_SampleA/outs/per_sample_outs/merged_bam/out_sorted.bam
samtools sort output_SampleB/outs/per_sample_outs/merged_bam_sampleB/out.bam -o output_SampleB/outs/per_sample_outs/merged_bam_sampleB/out_merged.bam

# then run velocyto again
time velocyto run -o velocyto_sampleB output_SampleB/outs/per_sample_outs/merged_bam_sampleB/out_merged.bam /omics/odcf/analysis/OE0228_projects/endothelial_cells/s467i/scRNA-seq/try_cellranger/refdata-gex-mm10-2020-A/genes/genes.gtf
time velocyto run -o velocyto_sampleA output_SampleA/outs/per_sample_outs/merged_bam/out_sorted.bam /omics/odcf/analysis/OE0228_projects/endothelial_cells/s467i/scRNA-seq/try_cellranger/refdata-gex-mm10-2020-A/genes/genes.gtf
#########################################################################################
# 13-02-2023

adata_A = dyn.read_loom("velocyto_sampleA/out_sorted_2WDZZ.loom")
adata_A.obs.index = [word.replace('out_sorted_2WDZZ:', '') for word in adata_A.obs.index]
adata_B = dyn.read_loom("velocyto_sampleB/out_merged_NZDOP.loom")
adata_B.obs.index = [word.replace('out_merged_NZDOP:', '') for word in adata_B.obs.index]

# # look for common genes (removed because they were the same)
# var_names = adata_B.var_names.intersection(adata_A.var_names)
# adata_A = adata_A[:, var_names]
# adata_B = adata_B[:, var_names]

#########################################################################################################
# 17-02-2023
def single_cell_processing(adata):
  # processing
  adata.var_names_make_unique() 
  #sc.pl.highest_expr_genes(adata, n_top=20, )
  sc.pp.filter_cells(adata, min_genes=200)
  sc.pp.filter_genes(adata, min_cells=3)
  ## quality control
  adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
  sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
  #sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)
  adata = adata[adata.obs.pct_counts_mt < 5, :]  # subset for the cells without high mt genes
  ## Normalizing and scaling and PCA
  sc.pp.normalize_total(adata, target_sum=1e4)
  sc.pp.log1p(adata)
  sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
  #sc.pl.highly_variable_genes(adata)
  sc.pp.regress_out(adata, ['total_counts'])  #Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed
  sc.pp.scale(adata, max_value=5)
  sc.tl.pca(adata, svd_solver='arpack', use_highly_variable=False)
  #sc.pl.pca(adata, color='CST3')
  #sc.pl.pca_variance_ratio(adata, log=True)
  return adata

# call this function on both sample objects
adata_A = single_cell_processing(adata_A)
adata_B = single_cell_processing(adata_B)

## INTEGRATION PART
# find the common genes and subset both objects
var_names = adata_A.var_names.intersection(adata_B.var_names)
adata_A = adata_A[:, var_names]
adata_B = adata_B[:, var_names]

adata = adata_A.concatenate(adata_B, batch_categories=['A', 'B'])

######################################################################

### AFTER MERGING THE TWO OBJECTS
# neighborhood and clustering
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata)
sc.tl.paga(adata)
sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
# Featureplot
sc.pl.umap(adata, color=['batch', 'leiden', 'Pecam1', 'Terf1'])

# Finding marker genes
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

######################################################################,   ,  
# tentative cluster annotation
new_cluster_names = ['Cap', 'Prol', 'Cap', 'Artery', 'Venous', "Venous", "Stalk", "Venous", "Venous", "Cap", "Cap", "ND"]
old_to_new = {'0':'Cap', '1':'Prol', '2':'Cap', '3':'Artery', '4':'Venous', '5':"Venous", '6':"Stalk", '7':"Venous", '8':"Venous", '9':"Cap", '10':"Cap", '11':"ND"}
adata.obs['new_clusters'] = (adata.obs['leiden'].map(old_to_new).astype('category'))
# old and new UMAP
sc.pl.umap(adata, color=['batch', 'leiden','new_clusters'])

######################################################################
### add cre information to the .obs and separate cre neg and pos

# adata.obs_names.to_list()
#['_'.join(i.split('_')[:2]) for i in adata.obs_names.to_list()[:5]] 
adata.obs['sample_hash'] = ['_'.join(i.split('_')[:2]) for i in adata.obs_names.to_list()]
hash_dict = {"sampleA_hash1": "neg", "sampleA_hash3": "neg", "sampleA_hash4": "neg", "sampleA_hash6":"neg", "sampleA_hash9":"neg",
"sampleB_hash3": "neg", "sampleB_hash4": "neg", "sampleB_hash6":"neg", "sampleB_hash8" :"neg", "sampleA_hash2":"pos", "sampleA_hash5":"pos", "sampleA_hash7": "pos", "sampleA_hash8": "pos", "sampleA_hash10":"pos", "sampleB_hash1":"pos","sampleB_hash2":"pos", "sampleB_hash5":"pos", "sampleB_hash7":"pos", "sampleB_hash9":"pos", "sampleB_hash10":"pos"}
adata.obs['cre'] = (adata.obs['sample_hash'].map(hash_dict).astype('category'))
sc.pl.umap(adata, color=['leiden','new_clusters', 'sample_hash', 'cre'])

## separate
adata_neg = adata[adata.obs['cre'].isin(["neg"])]
adata_pos = adata[adata.obs['cre'].isin(["pos"])]

######################################################################
## DYNAMO
# Initial processsing + velocity calculations
def calculate_velocity(adata):
  dyn.pp.recipe_monocle(adata)
  dyn.tl.dynamics(adata, model='stochastic', cores=3)
  # dyn.pl.umap(adata, color='Cell_type')	
  dyn.tl.cell_velocities(adata)
  # confidence of gene-wise velocity
  dyn.tl.gene_wise_confidence(adata, group='leiden', lineage_dict={'1': ['0']})
  # dyn.pl.streamline_plot(adata, color=['leiden']) 			# velocity vector of each cell
  # Differential geometry analysis
  # 1. vector field learning in velocity space
  dyn.tl.cell_velocities(adata, basis='pca'); 				# first project the RNA velocities into PCA space
  dyn.vf.VectorField(adata, basis='pca', M=100)				# learns the vector field function in PCA space
  # 2. Velocity, acceleration and curvature ranking + plotting
  dyn.vf.rank_velocity_genes(adata, groups='leiden', vkey="velocity_S");	# rank genes based on their velocity matrix.
  rank_speed = adata.uns['rank_velocity_S'];
  rank_abs_speed = adata.uns['rank_abs_velocity_S'];			# save the speed ranking information to rank_speed or rank_abs_speed for future usages if needed
  dyn.vf.acceleration(adata, basis='pca')					# to compute acceleration for each cell with the learned vector field in adata2
  dyn.vf.rank_acceleration_genes(adata, groups='leiden', akey="acceleration", prefix_store="rank");	# rank acceleration same as velocity 
  rank_acceleration = adata.uns['rank_acceleration'];
  rank_abs_acceleration = adata.uns['rank_abs_acceleration'];
  dyn.vf.curvature(adata, basis='pca');					# ranks genes based on their raw or absolute curvature values in different cell groups.
  dyn.vf.rank_curvature_genes(adata, groups='leiden');
  # dyn.pl.umap(adata, color=['Filip1l', 'Fbxl7'], layer='velocity_S', frontier=True)			# velocity plot
  # dyn.pl.umap(adata, color=['Filip1l', 'Fbxl7'], layer='acceleration', frontier=True)		# acceleration plot
  # dyn.pl.umap(adata, color=['Filip1l', 'Fbxl7'], layer='curvature', frontier=True)			# curvature plot
  return adata

adata_pos = calculate_velocity(adata_pos)
adata_neg = calculate_velocity(adata_neg)

adata_pos.var.use_for_pca[["Slc7a1","Slc7a5","Slc38a3","Mtor","Rps6","Bax"]] = ["True", "True","True", "True","True", "True"]
adata_neg.var.use_for_pca[["Slc7a1","Slc7a5","Slc38a3","Mtor","Rps6","Bax"]] = ["True", "True","True", "True","True", "True"]

for eachgene in adata_pos.var_names:
  if eachgene in ["Slc7a1","Slc7a5","Slc38a3","Mtor","Rps6","Bax"]:
    adata_pos.var.use_for_pca[eachgene] = True


# Jacobian Calculation and Ranking
# 1. selection of top genes
dyn.pp.top_pca_genes(adata_pos, n_top_genes=1000);				# for Jacobian matrix calculation, to reduce the computational power used
dyn.pp.top_pca_genes(adata_neg, n_top_genes=1000);
# take unique of all pca genes
top_pca_genes = list(set(adata_pos.var.index[adata_pos.var.top_pca_genes].to_list() + adata_neg.var.index[adata_neg.var.top_pca_genes].to_list() + ["Slc7a1", "Slc7a5", "Slc38a3", "mTOR", "S6", "Bax"]))
#top_pca_genes = ["erbb3b", "col6a3", "vwa1", "slc35c2", "col6a2", "col6a1"] + list(top_pca_genes)

dyn.vf.jacobian(adata_pos, regulators=top_pca_genes, effectors=top_pca_genes);
dyn.vf.jacobian(adata_neg, regulators=top_pca_genes, effectors=top_pca_genes);

# recalculating PCA because required genes were not present in PCA genes
sc.tl.pca(adata_pos, svd_solver='arpack', use_highly_variable=False)
sc.tl.pca(adata_neg, svd_solver='arpack', use_highly_variable=False)

################################################################
# in-silico perturbations
gene = "Slc38a3"
dyn.pd.perturbation(adata_pos, gene, [100], emb_basis="umap")
dyn.pl.streamline_plot(adata_pos, color=["new_clusters", gene], basis="umap_perturbation")

################################################################

# writing output files using cellbrowser module (didn't work)
import cellbrowser
sc.external.exporting.cellbrowser(adata, "scanpy_export", "denise_10x")

################################################################
# commented part

# # read seurat cellnames to subset
# cells_A = [line.rstrip() for line in open('cellnames_SampleA.csv')]
# cells_B = [line.rstrip() for line in open('cellnames_SampleB.csv')]
# genes = [line.rstrip() for line in open('genenames.csv')]
# 
# # to check if unassigned reads are somehow getting stored somewhere
# b = dyn.read_loom("velocyto_sampleB/unassigned_alignments_LZRTD.loom")
# b.obs.index = [word.replace('unassigned_alignments_LZRTD:', '') for word in b.obs.index]
