#
# Covid19 Healthy Lung Atlas - integrate ---- 
# this script is used to integrate sn/scRNA-seq datasets after processing..
# Benchmarking - Seurat V4 RPCA approach, Seurat V4 RPCA + Harmony, Harmony alone.. ----
# https://satijalab.org/seurat/articles/integration_rpca.html

# load libs ----
suppressPackageStartupMessages({
  library(ggplot2)
  library(SingleCellExperiment)
  library(scater)
  library(Seurat)
  library(dplyr)
  library(png)
  library(cowplot)
  library(parallel)
  library(harmony)
  library(cetcolor)
  library(gridExtra)
  library(mixOmics)
  #library(plotly)
  #library(SeuratDisk)
})

# set some R options
# visualize huge matrix, standard notation vs scientific one, & digits after comma
options(max.print = 6e5, scipen = 500, digits = 6) 

# parallelization in Seurat with future ----
library(future)
# set maximum for object size
# https://satijalab.org/seurat/articles/future_vignette.html
options(future.globals.maxSize = 100000 * 1024^2) # ~104G RAM
parallel::detectCores() # how many cores in total
#plan("multiprocess", workers = 22)
#nbrOfWorkers() # how many cores are set..
#plan() # check 

# load filtered (doublets removed) Lung COVID19-Healthy Atlas obj ----
dir_working <- "/rds/general/user/ahuseyno/projects/cardiac_single_cell_biology/live/sanger_covid19_snRNAseq/lung/lung_v2/"
setwd(dir_working)
object <- readRDS("./lung_atlas_filtered_scrublet.rds")
# clean up
object <- DietSeurat(object, counts = T, assays = "RNA", scale.data = F)
object

# split by "scn_source_study.batch" ie "scn_source" per "study batch" ----
object_split_list <- SplitObject(object, split.by = "scn_source_study.batch")
object_split_list %>% names() %>% length()
# clean up
rm(object)
gc()

##
# Seurat v4 (SCT+RPCA) ----
# correct for sc-sn per study batch c("scn_source_study.batch"), regressing out c("percent_mito","percent_ribo")
# run SCTransform on each obj in the list ----
gc()                                        
#plan("multiprocess", workers = 20)
#plan()
for (i in seq(names(object_split_list))) {
  object_split_list[[i]] <- SCTransform(object_split_list[[i]], 
                                        variable.features.n = 3000,
                                        do.correct.umi = T, 
                                        method = "glmGamPoi", residual_type = "pearson",
                                        vars.to.regress = c("percent_mito","percent_ribo"),
                                        conserve.memory = T, verbose = T)
}
# prepare for integration ----
gc()
object_split_list.features <- SelectIntegrationFeatures(object.list = object_split_list, nfeatures = 3000)
object_split_list <- PrepSCTIntegration(object.list = object_split_list, anchor.features = object_split_list.features)
# run PCA on each obj in the list ----
for (i in seq(names(object_split_list))) {
  object_split_list[[i]] <- RunPCA(object_split_list[[i]], npcs = 50, verbose = F, features = object_split_list.features)
}

# reference is Melms21 dataset ----
# covid lung atlas paper - https://www.nature.com/articles/s41586-021-03569-1#Sec1
reference_dataset <- grep(pattern = "Melms_2021", x = names(object_split_list))

# find integration anchors ----
# basically MNNs between batches in shared PCA space..
gc()                                        
#plan("multiprocess", workers = 22)
#plan()                                       
object_split_list.anchors <- FindIntegrationAnchors(object.list = object_split_list, 
                                                    normalization.method = "SCT",
                                                    anchor.features = object_split_list.features, 
                                                    dims = 1:30, reduction = "rpca", k.anchor = 20,
                                                    reference = reference_dataset)

# integrate ----
gc()
object <- IntegrateData(anchorset = object_split_list.anchors, normalization.method = "SCT", dims = 1:30)
DefaultAssay(object)
# keep only RNA and integrated assays ----
object <- DietSeurat(object, counts = T, assays = c("RNA","integrated"), data = T, scale.data = T)
object
rm(list = c("object_split_list.anchors","object_split_list.features","object_split_list"))
gc()

# run PCA ----
DefaultAssay(object) <- "integrated"
object <- RunPCA(object, npcs = 50, verbose = F)
# calculate the explained variance per PC
# add PC vars to @misc list.
object@misc <- c(c(object@reductions$pca@stdev)^2 / sum(matrixStats::rowVars(GetAssayData(object, assay = "integrated", slot = "scale.data")))) %>% list()
object@misc[[1]] <- (object@misc[[1]]/sum(object@misc[[1]]))*100 # get variance in %
names(object@misc[[1]]) <- paste0("PC",seq(length(object@reductions$pca@stdev)),rep("_var",length(object@reductions$pca@stdev)))
# add names to the list
names(object@misc)[1] <- "PCA_var_integrated"

# run IPCA from mixomics ----
#library(mixOmics)
DefaultAssay(object) <- "integrated"
#t(GetAssayData(object, assay = "integrated", slot = "data")) %>% str()
#t(GetAssayData(object, assay = "integrated", slot = "scale.data")) %>% str()
gc()
ipca_dim <- ipca(t(GetAssayData(object, assay = "integrated", slot = "scale.data")),
                 ncomp = 15, scale = F, fun = "logcosh", mode = "deflation", max.iter = 100)
#ipca_dim %>% str()
# add IPCA to object ----
object[["ipca"]] <- CreateDimReducObject(embeddings = ipca_dim$x, loadings = ipca_dim$loadings$X,
                                         stdev = ipca_dim$kurtosis,
                                         key = "IPC_", assay = "integrated")

# run UMAP on PCA ----
gc()
object <- RunUMAP(object, reduction.name = "umap.seurat", n.epochs = 400, spread = 1, min.dist = 0.3, n.neighbors = 30L, umap.method = "uwot", dims = 1:30, reduction = "pca", n.components = 3L, verbose = F)

# run UMAP on IPCA ---- 
DefaultAssay(object) <- "integrated"
gc()
object <- RunUMAP(object, reduction.name = "umap.seurat.ipca",  umap.method = "uwot",
                  n.epochs = 400, spread = 1, min.dist = 0.3, dims = 1:15, n.neighbors = 30L,
                  reduction = "ipca", n.components = 3L,
                  verbose = F)

# run NormalizeData on RNA assay ----
gc()                                        
#plan("multiprocess", workers = 20)
#plan()
DefaultAssay(object) <- "RNA"
object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000, verbose = T) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 10000) #%>% ScaleData(verbose = T)

# run SCTransform on RNA assay ----
# regress out some vars as needed
gc()                                        
#plan("multiprocess", workers = 20)
#plan()
DefaultAssay(object) <- "RNA"                                       
object <- SCTransform(object, new.assay.name = "SCT", method = "glmGamPoi", residual_type = "pearson", variable.features.n = 10000, 
                      vars.to.regress = c("percent_mito","percent_ribo"),
                      do.correct.umi = T, conserve.memory = T, verbose = T)
object

# save temp integrated object ----
# check if dir exist
dir_path_temp <- "/rds/general/user/ahuseyno/projects/cardiac_single_cell_biology/ephemeral/ah_covid/"
if (!dir.exists(dir_path_temp)) { dir.create(path = dir_path_temp) }
gc()
# save on RDS ./ephemeral
saveRDS(object, file = paste0(dir_path_temp,"lung_atlas_integrated_temp_test.rds"))
# save on RDS ./live
#saveRDS(object, file = paste0(dir_working,"lung_atlas_integrated.rds"))

## 
# Seurat v4 (SCT+RPCA) + Harmony ----
# use integrated assay, ie corrected for sc-sn per study batch c("scn_source_study.batch"), regressing out c("percent_mito","percent_ribo")
# 2nd step - Harmony - correct PCA for c("donor_study.batch","version_10X")

# run Harmony integration - as 2nd correction ----
DefaultAssay(object) <- "integrated"
group.by.vars <- c("donor_study.batch","version_10X")
gc()
object <- object %>% harmony::RunHarmony(
  group.by.vars <- group.by.vars, assay.use = "integrated",
  theta = rep(2,length(group.by.vars)), lambda = rep(1,length(group.by.vars)),
  reduction = "pca", reduction.save = "harmony", dims.use = 1:30, max.iter.harmony = 30, plot_convergence = F)

# Run Harmony integration - 2nd correction for IPCA ----
#DefaultAssay(object) <- "integrated"
#group.by.vars <- c("donor_study.batch","version_10X")
gc()
object <- object %>% harmony::RunHarmony(
  group.by.vars <- group.by.vars, assay.use = "integrated",
  theta = rep(2,length(group.by.vars)), lambda = rep(1,length(group.by.vars)),
  reduction = "ipca", reduction.save = "harmony.ipca", dims.use = 1:15, max.iter.harmony = 30, plot_convergence = F)
object

# run UMAP on Harmony-corrected PCA ----
gc()
object <- RunUMAP(object, reduction.name = "umap", n.epochs = 400, spread = 1, min.dist = 0.3, 
                  n.neighbors = 30L, 
                  umap.method = "uwot", dims = 1:30, 
                  reduction = "harmony", n.components = 3L, verbose = F)

# run UMAP on Harmony-corrected IPCA ----
gc()
object <- RunUMAP(object, reduction.name = "umap.harmony.ipca",  umap.method = "uwot",
                  n.epochs = 300, spread = 1, min.dist = 0.3, dims = 1:15, n.neighbors = 30L,
                  reduction = "harmony.ipca", n.components = 3L,
                  verbose = F)
object

# save temp obj on RDS ./ephemeral ----
saveRDS(object, file = paste0(dir_path_temp,"lung_atlas_integrated_temp2_test.rds"))

##
# run Harmony solo integration (SCT assay) ----
# correct PCA for c("scn_source_study.batch","donor_study.batch","version_10X")

# run PCA on SCT assay ----
DefaultAssay(object) <- "SCT"
object <- RunPCA(object, npcs = 50, assay = "SCT", verbose = F, reduction.name = "pca.uncorrected")

# run IPCA on SCT assay ----
gc()
ipca_dim <- ipca(t(GetAssayData(object, assay = "SCT", slot = "scale.data")),
                 ncomp = 15, scale = F, fun = "logcosh", mode = "deflation", max.iter = 100)
#ipca_dim %>% str()
# add IPCA to object ----
object[["ipca.uncorrected"]] <- CreateDimReducObject(embeddings = ipca_dim$x, loadings = ipca_dim$loadings$X,
                                         stdev = ipca_dim$kurtosis,
                                         key = "IPC_uncorrected_", assay = "SCT")
object
# variables to correct
group.by.vars <- c("scn_source_study.batch","donor_study.batch","version_10X")
gc()
object <- object %>% harmony::RunHarmony(
  group.by.vars <- group.by.vars, assay.use = "SCT",
  theta = rep(2,length(group.by.vars)), lambda = rep(1,length(group.by.vars)),
  reduction = "pca.uncorrected", reduction.save = "harmony.solo.pca", dims.use = 1:30, max.iter.harmony = 30, plot_convergence = F)

#..on ipca
object <- object %>% harmony::RunHarmony(
  group.by.vars <- group.by.vars, assay.use = "SCT",
  theta = rep(2,length(group.by.vars)), lambda = rep(1,length(group.by.vars)),
  reduction = "ipca.uncorrected", reduction.save = "harmony.solo.ipca", dims.use = 1:15, max.iter.harmony = 30, plot_convergence = F)

# run UMAP on Harmony-corrected PCA ----
gc()
object <- RunUMAP(object, reduction.name = "umap.harmony.solo.pca", n.epochs = 400, spread = 1, min.dist = 0.3, n.neighbors = 30L, 
                  umap.method = "uwot", dims = 1:30, reduction = "harmony.solo.pca", n.components = 3L, verbose = F)

# run UMAP on Harmony-corrected IPCA ----
gc()
object <- RunUMAP(object, reduction.name = "umap.harmony.solo.ipca", n.epochs = 400, spread = 1, min.dist = 0.3, n.neighbors = 30L, umap.method = "uwot", dims = 1:15, reduction = "harmony.solo.ipca", n.components = 3L, verbose = F)
object

# save integrated object ----
# check if dir exist
dir_path_temp <- "/rds/general/user/ahuseyno/projects/cardiac_single_cell_biology/ephemeral/ah_covid/"
if (!dir.exists(dir_path_temp)) { dir.create(path = dir_path_temp) }
gc()
# save on RDS ./ephemeral
saveRDS(object, file = paste0(dir_path_temp,"lung_atlas_integrated_test.rds"))
# save on RDS ./home
saveRDS(object, file = paste0("/rds/general/user/ahuseyno/home/projects_imperial/covid19_project/lung_atlas_integrated_test.rds"))
# save on RDS ./live
saveRDS(object, file = paste0(dir_working,"lung_atlas_integrated_test.rds"))
