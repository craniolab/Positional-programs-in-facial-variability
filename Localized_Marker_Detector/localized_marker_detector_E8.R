library(Seurat)
library(SeuratWrappers)
library(anndata)
library(LocalizedMarkerDetector)
library(ggplot2)
library(patchwork)
library(RhpcBLASctl)
options(renv.config.cache.symlinks=FALSE)
blas_set_num_threads(16)
setwd('/mnt/mpistaff/Cranio_Lab/Louk_Seton/mesenchyme_project_2023/')

ad <- read_h5ad('anndata_objects/interactive/mesen_E8_new.h5ad')
cell_cycle_genes = read.delim("cell_cycle_genes_mesen.txt",header = FALSE)
cell_cycle_genes = cell_cycle_genes$V1
mesen.seurat <- CreateSeuratObject(counts = t(as.matrix(ad$layers['original_counts'])), meta.data = ad$obs,)
mesen.seurat[['RNA']]@meta.data = ad$var

X_umap = ad$obsm$X_umap
X_pca = ad$obsm$X_pca
dimnames(X_umap)[[2]]=c('UMAP_1','UMAP_2')
dimnames(X_pca)[[2]]=sprintf("PCA_%d", 1:dim(X_pca)[[2]])
dimnames(X_umap)[[1]]=rownames(ad$obs)
dimnames(X_pca)[[1]]=rownames(ad$obs)
mesen.seurat@reductions[["umap"]] = CreateDimReducObject(embeddings = as.matrix(X_umap), key = "UMAP_", assay = DefaultAssay(mesen.seurat)) 
mesen.seurat@reductions[["pca"]] = CreateDimReducObject(embeddings = as.matrix(X_pca), key = "PCA_", assay = DefaultAssay(mesen.seurat)) 

#mesen.seurat = subset(x = mesen.seurat, subset = sample %in% c('10','11','12','13','14'))
mesen.seurat = subset(mesen.seurat, features = rownames(mesen.seurat)[!rownames(mesen.seurat) %in% cell_cycle_genes]) #remove cell cycle genes


mesen.seurat <- NormalizeData(mesen.seurat, normalization.method = "LogNormalize", scale.factor = 10000)

DimPlot(object = mesen.seurat, reduction = "umap",group.by='prominences_merged')
FeaturePlot(object = mesen.seurat, features = "Alx4")

DefaultAssay(mesen.seurat) <- "RNA"
n_dim = dim(mesen.seurat@reductions$pca@cell.embeddings)[2]
feature_space = as.matrix(mesen.seurat@reductions$pca@cell.embeddings[,1:n_dim])
visual_space = data.frame(mesen.seurat@reductions$umap@cell.embeddings)
dat = as.matrix(mesen.seurat[[DefaultAssay(mesen.seurat)]]$data)
cell_label = mesen.seurat$prominences_merged

Gene_detected_count <- apply(dat > apply(dat,2,median),1,sum)
selected_genes = (Gene_detected_count >= 30) & (Gene_detected_count <= ncol(dat) * 0.66)
dat = dat[selected_genes,,drop = FALSE]

# Construct knn graph
knn_result = ConstructKnnGraph(knn = 5, feature_space = feature_space)
saveRDS(knn_result, 'r_env/knn_result_E8.rds')
#> Constructing KNN graph
#knn_result = readRDS('r_env/knn_result.rds')
A = knn_result$adj_matrix # Adjacency Matrix
W = knn_result$graph # Symmetrized Graph ((A + AT) / 2)

# saveRDS(A, 'r_env/A.rds')
# saveRDS(W, 'r_env/W.rds')


# Plot knn graph
#VisualizeGraph(affinity_m = W, label = cell_label, layout = visual_space) + 
#  guides(color = guide_legend(ncol = 1, byrow = TRUE)) + 
#  theme(
#    legend.title = element_text(size = rel(0.7)),
#    legend.text = element_text(size = rel(0.7)),
#    legend.key.height = unit(1, "null"))

#had to create a custom function because diag doesn't work on such a large input if it's not a matrix.

Doubly_stochastic <- LocalizedMarkerDetector:::Doubly_stochastic
sinkhorn_knopp <- LocalizedMarkerDetector:::sinkhorn_knopp
sinkhorn_knopp_largedata <- LocalizedMarkerDetector:::sinkhorn_knopp_largedata

ConstructDiffusionOperators_diag <- function(W, max_time){
  cat("Create a list of diffusion operators...\n")
  pb <- txtProgressBar(min = 0, max = 100, style = 3)
  P = Doubly_stochastic(W)
  P = as.matrix(P)
  setTxtProgressBar(pb, 10)
  # P = RowwiseNormalize(W)
  # eig_res = RSpectra::eigs_sym(P, k = 1, which = "LM")
  max_step = max_time
  # P_ls = list(P)
  P_ls = NULL
  if(max_step < 1){
    stop("Incorrect diffusion time, no propogation")
  }else{
    t = 1
    # automatic decide max_step by checking whether the diag of P -> 1/n
    while(t <= floor(log(max_step,2)) & 
          max(abs((ncol(P) * diag(P)) * ncol(P) - ncol(P))) >= 1e-2 * ncol(P)){  ######## I added the "as.matrix"
      P = crossprod(P); t = t + 1
      P_ls = c(P_ls,list(as(P,"sparseMatrix")))
      setTxtProgressBar(pb, 10 * (t + 1))
    }
    setTxtProgressBar(pb, 100)
  }
  # Add t = 0
  W.diag = as(diag(nrow(W)),"sparseMatrix")
  P_ls = c(list(W.diag),P_ls)
  # make it sparse
  cat("\nConverting diffusion operators to sparse matrices...\n")
  # P_ls = lapply(P_ls,function(x) as(x,"sparseMatrix"))
  names(P_ls) = c(0,2^seq(1,t-1))
  cat("\nMax diffusion time:",2^(t-1),"\n")
  return(P_ls)
}
P_ls = ConstructDiffusionOperators_diag(W = W, max_time = 2^6)

rho = RowwiseNormalize(dat[,colnames(W)])

gene = "Alx4"
VisualizeDiffusion(coord = visual_space,init_state = rho[gene,],P_ls = P_ls,check_time = c(0,2,2^4,2^6),gene_name = gene) & 
  theme(
    plot.title = element_text(size = rel(0.7)),
    plot.subtitle = element_text(size = rel(0.7)),
    legend.title = element_text(size = rel(0.7)),
    legend.text = element_text(size = rel(0.7)) )

score_result = fast_get_lmds(W = W, init_state = rho, max_time = 2^6,P_ls = P_ls)

lmd.cumulative = as.data.frame(score_result$cumulative_score)

res = show_result_lmd(score_result, kneeplot = TRUE)
head(res$gene_table,10)

res.cutoff = res$gene_table[1:length(res$cut_off_gene),]
res.cutoff[c('Alx4'),]
FeaturePlot(object = mesen.seurat, features = rownames(head(res$gene_table,10)))

top_lmd_genes = names(res$cut_off_gene)

# tic()
# ALRA imputation
if(!"alra" %in% names(mesen.seurat@assays)){
  mesen.seurat <- RunALRA(mesen.seurat, assay = "RNA")
  DefaultAssay(mesen.seurat) <- "RNA"
}
dat_alra = as.matrix(mesen.seurat[["alra"]]@data)[top_lmd_genes,]

# Compute the gene-gene pairwise distance
dist = CalculateGeneDistance(dat_alra, method = "jaccard")
# toc()

gene_cl_res = ClusterGenes(dist, accu=.55,min_gene = 10,clustering_method = "complete", return_tree = TRUE, deepSplit = 1)

VisualizeGeneHeatmap(dist, gene_partition = gene_cl_res$gene_partition, gene_hree = gene_cl_res$gene_hree)

gene_partition = gene_cl_res$gene_partition
levels(gene_partition) = 1:nlevels(gene_partition) # rename gene modules

AddModuleActivityScore_v5 <- function(srat, gene_partition, assay = "RNA", 
                                      do_local_smooth = FALSE, knn = 10, major_vote = 6, 
                                      nloop = 100, module_names = NULL){
  # adjust the format of gene_partition
  gene_partition = setNames(as.character(gene_partition),names(gene_partition))
  gene_partition = gene_partition[names(gene_partition) %in% rownames(srat)]
  gene_partition = gene_partition[!is.na(gene_partition)]
  gene_partition = as.factor(gene_partition)
  dat = as.matrix(srat[[assay]]$data[names(gene_partition),])
  if(do_local_smooth){
    if ("pca" %in% names(srat@reductions)){
      ndims = FindPC(srat)
      feature_space <- Embeddings(srat, reduction = "pca")[,1:ndims]
      cat(ncol(feature_space),"PC used for building graph for majority vote\n")
    }else if ("lsi" %in% names(srat@reductions)){
      ndims = 50
      feature_space <- Embeddings(srat, reduction = "lsi")[,2:ndims]
      cat(ncol(feature_space),"LSI used for building graph for majority vote\n")
    }
    else{
      stop("Please RunPCA")
    }
    A = ConstructKnnGraph(knn = knn, feature_space = feature_space)$'adj_matrix'
    cat(sprintf("Do Major Vote: at least %d out of %d neighbors(self-include) expressed\n", major_vote, knn + 1))
  }else{
    A = NULL
  }
  cat("Start Loop\n")
  pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                   total = nloop,
                                   complete = "=",   # Completion bar character
                                   incomplete = "-", # Incomplete bar character
                                   current = ">",    # Current bar character
                                   clear = FALSE,    # If TRUE, clears the bar when finish
                                   width = 100)
  cell_block_loop = NULL
  for(loop in 1:nloop){
    cell_block_loop = c(cell_block_loop,list(GMM_subsampling(seed = loop, 
                                                             gene_partition, 
                                                             expr_dat = dat, 
                                                             cell_kNN_graph = A, major_vote)))
    pb$tick()
  }
  cell_block_prop = Reduce(`+`, cell_block_loop) / length(cell_block_loop)
  cell_block = cell_block_prop
  if(is.null(module_names) | length(module_names) != ncol(cell_block)){
    colnames(cell_block) = paste0("Module",colnames(cell_block))
  }else{
    colnames(cell_block) = module_names
  }
  srat <- Seurat::AddMetaData(srat, cell_block, col.name = colnames(cell_block))
  return(srat)
}
GMM_subsampling <- LocalizedMarkerDetector:::GMM_subsampling

mesen.seurat = AddModuleActivityScore_v5(mesen.seurat, gene_partition = gene_partition)

# Visualize patterns
(FeaturePlot(mesen.seurat, features = colnames(mesen.seurat@meta.data)[grepl("Module",colnames(mesen.seurat@meta.data))], order = TRUE, reduction = "umap", ncol = 5) & NoAxes() & 
    scale_color_gradient(low = "lightgrey", high = "blue", limits = c(0,1)) & labs(color = "ModuleScore") & theme(
      plot.title = element_text(size = rel(0.5)),
      legend.title = element_text(size = rel(0.5)),
      legend.text = element_text(size = rel(0.5)) )) + plot_layout(guides = "collect")

module_id = 9
p1 = FeaturePlot(mesen.seurat, features = paste0("Module",module_id), order = TRUE, reduction = "umap") + ggtitle(sprintf("Module%s (%d genes)",module_id,sum(gene_partition == module_id))) + labs(color = "ModuleScore") + NoAxes()
pl = FeaturePlot(mesen.seurat, features = names(gene_partition)[gene_partition == module_id][1:12],ncol = 3, order = TRUE) & NoAxes()
p = (p1 + pl + plot_layout(design = c("#BBB\nABBB\nABBB\n#BBB"))) & theme(
  plot.title = element_text(size = rel(0.5)),
  legend.title = element_text(size = rel(0.5)),
  legend.text = element_text(size = rel(0.5)) )
p

write.csv(res$gene_table, 'r_lmd_result_E8.txt')
