# Analysis of scRNA seq data for supporting cells in 3D (MSC, HUVEC) - Magdalini Kanari

## load packages ----
library(Seurat)
library(cowplot)
library(patchwork)
library(ggplot2)
library(dplyr)
library(ragg)
library(qusage)

## create seurat object of all supporting cells - samples: MSC, HUVEC, MSC_HUVEC, CO1, CO2, CO3, CO4, CO5, CO6 ----
dirs1 <- list.dirs("~/kispi/Research/12_Sequencing_raw/MagdaliniKanari/ECTICA_scRNAseq/Controls/", recursive = FALSE, full.names = FALSE)
dirs2 <- list.dirs("~/kispi/Research/12_Sequencing_raw/MagdaliniKanari/ECTICA_scRNAseq/Leukemic/co/", recursive = FALSE, full.names = FALSE)

# controls
for(name in dirs1){
  
  cts1 <- Read10X(data.dir = paste0("~/kispi/Research/12_Sequencing_raw/MagdaliniKanari/ECTICA_scRNAseq/Controls/", name, "/"))
  
  assign(name, CreateSeuratObject(counts = cts1, 
                                  min.cells = 3, 
                                  min.features = 200))}
rm("2d", "cts1", dirs1)

# cocultures
for(name in dirs2){
  
  cts2 <- Read10X(data.dir = paste0("~/kispi/Research/12_Sequencing_raw/MagdaliniKanari/ECTICA_scRNAseq/Leukemic/co/", name, "/"))
  
  assign(name, CreateSeuratObject(counts = cts2, 
                                  min.cells = 3, 
                                  min.features = 200))}

rm("cts2", dirs2, name)

# set correct identities
new_identity <- "msc"
msc$orig.ident <- new_identity

new_identity <- "mschuvec"
mschuvec$orig.ident <- new_identity

new_identity <- "huvec"
huvec$orig.ident <- new_identity

for (i in 1:6) {
  new_identity <- paste0("co", i)
  obj <- get(new_identity)  
  obj$orig.ident <- new_identity  
  assign(new_identity, obj)}

rm(i, "obj", new_identity)

## filter low quality cells (mitochondria genes, nFeature, nCount) ----
## check every sample independently, thresholds below were used 
sample_names <- c(paste0("co", 1:6), "msc", "huvec", "mschuvec")
for (sample in sample_names) {
  obj <- get(sample)
  obj[["percent.mito"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  upper_lim <- quantile(obj$nFeature_RNA, 0.98)
  obj <- subset(obj, subset = nFeature_RNA < upper_lim & nCount_RNA < 50000 & percent.mito < 5)
  assign(sample, obj)}

rm("obj")

## merge and process the cocultures ----
co <- merge(co1, y = c(co2, co3, co4, co5, co6), add.cell.ids = c("CO1", "CO2", "CO3", "CO4", "CO5", "CO6"), project = "co")
co <- NormalizeData(co)
co <- FindVariableFeatures(co, selection.method = "vst")
co <- ScaleData(co)
co <- RunPCA(co, features = VariableFeatures(object = co))
co <- FindNeighbors(co, dims = 1:20)
co <- FindClusters(co, resolution = 0.1)
co <- RunUMAP(co, dims = 1:20)
DimPlot(co, reduction = "umap", label = TRUE)

## merge and process all seurats ----
sup_cells_MSCHUVEC <- merge(msc, y = c(mschuvec, huvec, co), add.cell.ids = c("MSC", "MSC_HUVEC", "HUVEC", "CO"), project = "sup_cells_MSCHUVEC")
sup_cells_MSCHUVEC <- NormalizeData(sup_cells_MSCHUVEC)
sup_cells_MSCHUVEC <- FindVariableFeatures(sup_cells_MSCHUVEC, selection.method = "vst")
sup_cells_MSCHUVEC <- ScaleData(sup_cells_MSCHUVEC)
sup_cells_MSCHUVEC <- RunPCA(sup_cells_MSCHUVEC, features = VariableFeatures(object = sup_cells_MSCHUVEC))
sup_cells_MSCHUVEC <- FindNeighbors(sup_cells_MSCHUVEC, dims = 1:20)
sup_cells_MSCHUVEC <- FindClusters(sup_cells_MSCHUVEC, resolution = 0.1)
sup_cells_MSCHUVEC <- RunUMAP(sup_cells_MSCHUVEC, dims = 1:20)
DimPlot(sup_cells_MSCHUVEC, reduction = "umap", label = TRUE)

## remove leukemia cells ----
FeaturePlot(sup_cells_MSCHUVEC, features = c("CD7","CD5", "CD22", "CD19"))
sup_cells_MSCHUVEC <- subset(x = sup_cells_MSCHUVEC, idents = c(4,7,3,2), invert = TRUE)

# re-run pipeline without leukemic cells
sup_cells_MSCHUVEC <- NormalizeData(sup_cells_MSCHUVEC)
sup_cells_MSCHUVEC <- FindVariableFeatures(sup_cells_MSCHUVEC, selection.method = "vst")
sup_cells_MSCHUVEC <- ScaleData(sup_cells_MSCHUVEC)
sup_cells_MSCHUVEC <- RunPCA(sup_cells_MSCHUVEC, features = VariableFeatures(object = sup_cells_MSCHUVEC))
sup_cells_MSCHUVEC <- FindNeighbors(sup_cells_MSCHUVEC, dims = 1:20)
sup_cells_MSCHUVEC <- FindClusters(sup_cells_MSCHUVEC, resolution = 0.3)
sup_cells_MSCHUVEC <- RunUMAP(sup_cells_MSCHUVEC, dims = 1:20)
DimPlot(sup_cells_MSCHUVEC, reduction = "umap", label = TRUE)

# saveRDS(sup_cells_MSCHUVEC, "~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Analysis/Supporting_Cells/Combined_MSCHUVEC/sup_cells_MSCHUVECg.rds")

## divide into MSCs and HUVECs ----
FeaturePlot(sup_cells_MSCHUVEC, features = c("LEPR", "NES", "CXCL12", "PDGFRB", "CDH5", "PECAM1"))
DimPlot(sup_cells_MSCHUVEC, group.by = 'RNA_snn_res.0.1', label = TRUE)
Idents(sup_cells_MSCHUVEC) <- "RNA_snn_res.0.1"
# Select MSCS
sup_cells_MSCs <- subset(x = sup_cells_MSCHUVEC, idents = c(5, 8), invert = TRUE)
sup_cells_MSCs <- NormalizeData(sup_cells_MSCs)
sup_cells_MSCs <- FindVariableFeatures(sup_cells_MSCs, selection.method = "vst")
sup_cells_MSCs <- ScaleData(sup_cells_MSCs)
sup_cells_MSCs <- RunPCA(sup_cells_MSCs, features = VariableFeatures(object = sup_cells_MSCs))
sup_cells_MSCs <- FindNeighbors(sup_cells_MSCs, dims = 1:20)
sup_cells_MSCs <- FindClusters(sup_cells_MSCs, resolution = 0.3)
sup_cells_MSCs <- RunUMAP(sup_cells_MSCs, dims = 1:20)
DimPlot(sup_cells_MSCs, reduction = "umap", label = TRUE, group.by = "orig.ident")

# Further MSC decontamination
# CDH5 and CD19 detected in downstream analysis, so remove them
sup_cells_MSC_CDH5neg <-subset(x = sup_cells_MSCHUVEC, subset = CDH5 < 1)
sup_cells_MSC_CDH5neg <- NormalizeData(sup_cells_MSC_CDH5neg)
sup_cells_MSC_CDH5neg <- FindVariableFeatures(sup_cells_MSC_CDH5neg, selection.method = "vst")
sup_cells_MSC_CDH5neg <- ScaleData(sup_cells_MSC_CDH5neg)
sup_cells_MSC_CDH5neg <- RunPCA(sup_cells_MSC_CDH5neg, features = VariableFeatures(object = sup_cells_MSC_CDH5neg))
sup_cells_MSC_CDH5neg <- FindNeighbors(sup_cells_MSC_CDH5neg, dims = 1:20)
sup_cells_MSC_CDH5neg <- FindClusters(sup_cells_MSC_CDH5neg, resolution = 0.3)
sup_cells_MSC_CDH5neg <- RunUMAP(sup_cells_MSC_CDH5neg, dims = 1:20)
DimPlot(sup_cells_MSC_CDH5neg, reduction = 'umap', label = TRUE, group.by = "RNA_snn_res.0.3")

sup_cells_MSC_CDH5negCD19neg <-subset(x = sup_cells_MSC_CDH5neg, subset = CD19 < 0.1)
sup_cells_MSC_CDH5negCD19neg <- NormalizeData(sup_cells_MSC_CDH5negCD19neg)
sup_cells_MSC_CDH5negCD19neg <- FindVariableFeatures(sup_cells_MSC_CDH5negCD19neg, selection.method = "vst")
sup_cells_MSC_CDH5negCD19neg <- ScaleData(sup_cells_MSC_CDH5negCD19neg)
sup_cells_MSC_CDH5negCD19neg <- RunPCA(sup_cells_MSC_CDH5negCD19neg, features = VariableFeatures(object = sup_cells_MSC_CDH5negCD19neg))
sup_cells_MSC_CDH5negCD19neg <- FindNeighbors(sup_cells_MSC_CDH5negCD19neg, dims = 1:20)
sup_cells_MSC_CDH5negCD19neg <- FindClusters(sup_cells_MSC_CDH5negCD19neg, resolution = 0.3)
sup_cells_MSC_CDH5negCD19neg <- RunUMAP(sup_cells_MSC_CDH5negCD19neg, dims = 1:20)
DimPlot(sup_cells_MSC_CDH5negCD19neg, reduction = 'umap', label = TRUE, group.by = "orig.ident")

# saveRDS(sup_cells_MSC_CDH5negCD19neg, "~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Supporting_Cells/Combined_MSCHUVEC/Subset_from_MSCHUVEC/CDH5negCD19neg_MSCs_FINAL_USED/sup_cells_MSC_CDH5negCD19neg.rds")

# Select HUVECs
sup_cells_HUVEC <- subset(x = sup_cells_MSCHUVEC, idents = c(5, 8))
sup_cells_HUVEC <- NormalizeData(sup_cells_HUVEC)
sup_cells_HUVEC <- FindVariableFeatures(sup_cells_HUVEC, selection.method = "vst")
sup_cells_HUVEC <- ScaleData(sup_cells_HUVEC)
sup_cells_HUVEC <- RunPCA(sup_cells_HUVEC, features = VariableFeatures(object = sup_cells_HUVEC))
sup_cells_HUVEC <- FindNeighbors(sup_cells_HUVEC, dims = 1:20)
sup_cells_HUVEC <- FindClusters(sup_cells_HUVEC, resolution = 0.3)
sup_cells_HUVEC <- RunUMAP(sup_cells_HUVEC, dims = 1:20)
DimPlot(sup_cells_HUVEC, reduction = "umap", label = TRUE, group.by = "orig.ident")

# saveRDS(sup_cells_HUVEC, "~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Analysis/Supporting_Cells/Combined_MSCHUVEC/Subset_from_MSCHUVEC/HUVEC_condition_comparison_FINAL_USED/sup_cells_HUVEC.rds")

## Find markers for HUVECs ----
Idents(sup_cells_HUVEC) <- "orig.ident"
sup_cells_HUVEC_1 <-FindMarkers(sup_cells_HUVEC, ident.1 = c("HUVEC"), ident.2 = c("MSC_HUVEC"), subset.ident = 'orig.ident')
sup_cells_HUVEC_3 <-FindMarkers(sup_cells_HUVEC, ident.1 = c("MSC_HUVEC"), ident.2 = c("CO1", "CO2", "CO3", "CO4", "CO5", "CO6"), subset.ident = 'orig.ident')

# saveRDS (...)

## Find markers for MSCs -----
Idents(sup_cells_MSC_CDH5negCD19neg) <- "orig.ident"
sup_cells_MSC_CDH5negCD19neg_markers_all <- FindMarkers(sup_cells_MSC_CDH5negCD19neg, ident.1 = c("CO1", "CO2", "CO3", "CO4", "CO5", "CO6", "MSC_HUVEC"), 
                                                        ident.2 = "MSC" , subset.ident = 'orig.ident')

# saveRDS (...)

