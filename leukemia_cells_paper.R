# Analysis of scRNA seq data for leukemia cells in 3D and in vivo - Magdalini Kanari

## load packages ----
library(Seurat)
library(cowplot)
library(patchwork)
library(ggplot2)
library(dplyr)
library(ragg)
library(qusage)
library(paleteer)
library(DoubletFinder)
library(fgsea)

## create seurat objects for all samples ----
dirs1 <- list.dirs("~/kispi/Research/12_Sequencing_raw/MagdaliniKanari/ECTICA_scRNAseq/Leukemic/co/", recursive = FALSE, full.names = FALSE)
dirs2 <- list.dirs("~/kispi/Research/12_Sequencing_raw/MagdaliniKanari/ECTICA_scRNAseq/Leukemic/mono/", recursive = FALSE, full.names = FALSE)
dirs3 <- list.dirs("~/kispi/Research/12_Sequencing_raw/MagdaliniKanari/in_vivo/", recursive = FALSE, full.names = FALSE)

for(name in dirs1){
  
  cts1 <- Read10X(data.dir = paste0("~/kispi/Research/12_Sequencing_raw/MagdaliniKanari/ECTICA_scRNAseq/Leukemic/co/", name, "/"))
  
  assign(name, CreateSeuratObject(counts = cts1, 
                                  min.cells = 3, 
                                  min.features = 200))}

for(name in dirs2){
  
  cts2 <- Read10X(data.dir = paste0("~/kispi/Research/12_Sequencing_raw/MagdaliniKanari/ECTICA_scRNAseq/Leukemic/mono/", name, "/"))
  
  assign(name, CreateSeuratObject(counts = cts2, 
                                  min.cells = 3, 
                                  min.features = 200))}

for(name in dirs3){
  
  cts3 <- Read10X(data.dir = paste0("~/kispi/Research/12_Sequencing_raw/MagdaliniKanari/in_vivo/", name, "/"))
  
  assign(name, CreateSeuratObject(counts = cts3[["Gene Expression"]], 
                                  min.cells = 3, 
                                  min.features = 200))}

# assign correct identities
seurat_objects <- c(paste0("co", 1:6), paste0("mono", 1:6), paste0("IN_VIVO", 1:6))
for (object_name in seurat_objects) {
  obj <- get(object_name)
  Idents(obj) <- object_name
  assign(object_name, obj)}

# pre-processing with standard pipeline for each sample ----
for (object_name in seurat_objects) {
  obj <- get(object_name)
  
  # 1. QC filtering
  obj[["percent.mito"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  upper_lim <- quantile(obj$nFeature_RNA, 0.98)
  obj <- subset(obj, subset = nFeature_RNA < upper_lim)
  obj <- subset(obj, subset = nCount_RNA < 50000)
  obj <- subset(obj, subset = percent.mito < 5)
  
  # 2. standard workflow 
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst")
  all.genes <- rownames(obj)
  obj <- ScaleData(obj, features = all.genes)
  obj <- RunPCA(obj, features = VariableFeatures(object = obj))
  obj <- FindNeighbors(obj, dims = 1:20)
  obj <- FindClusters(obj)
  set.seed(42)
  obj <- RunUMAP(obj, dims = 1:20)
  
  assign(object_name, obj)}

# doublet finder for each sample  ----
# calculate doublet rate based on number of cells per sample, information online from doubletFinder
doublet_rates <- c(CO1 = 0.054, CO2 = 0.046, CO3 = 0.039, CO4 = 0.046, CO5 = 0.031, CO6 = 0.023,
                   MONO1 = 0.039, MONO2 = 0.004, MONO3 = 0.016, MONO4 = 0.031, MONO5 = 0.069, MONO6 = 0.039,
                   IN_VIVO_1 = 0.004, IN_VIVO_2 = 0.008, IN_VIVO_3 = 0.004, IN_VIVO_4 = 0.008, IN_VIVO_5 = 0.008, IN_VIVO_6 = 0.023)

# Calculate independently for each sample - CHANGE SEURAT OBJECT and  DOUBLET RATE
# check correct sample name (eg. CO2)
sweep.res.list <- paramSweep_v3(CO2, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn_CO2 <- find.pK(sweep.stats)

# Select pK that corresponds to max bcmvn to optimize doublet detection
pK <- bcmvn_CO2 %>% 
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

# Homotypic Doublet Proportion Estimate
annotations <- CO2@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
# check doublet rate (eg. 0.048)
nExp_poi <- round(0.048*nrow(CO2@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run doubletFinder adn remove doublets
CO2 <- doubletFinder_v3(CO2, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
CO2 <- subset(CO2, cells=rownames(CO2@meta.data)[which(CO2@meta.data$DF.classification == "Singlet")])

# assign cell cycle for each sample ----
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

for (object_name in seurat_objects) {
  obj <- get(object_name)
  obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  assign(object_name, obj)}

 
## merge seurat objects based on condition ----
mono <- merge(MONO1, y = c(MONO2, MONO3, MONO4, MONO5, MONO6), 
              add.cell.ids = c("MONO1", "MONO2", "MONO3", "MONO4", "MONO5", "MONO6"), project = "mono")
# saveRDS(mono, "~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Doublets_excluded/rds.files/mono.rds")

co <- merge(CO1, y = c(CO2, CO3, CO4, CO5, CO6), 
            add.cell.ids = c("CO1", "CO2", "CO3", "CO4", "CO5", "CO6"), project = "co")
# saveRDS(co, "~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Doublets_excluded/rds.files/co.rds")

invivo<- merge(IN_VIVO1, y = c(IN_VIVO2, IN_VIVO3, IN_VIVO4, IN_VIVO5, IN_VIVO6), 
               add.cell.ids = c("invivo1", "invivo2", "invivo3", 'invivo4', "invivo5", "invivo6"), project = "invivo")
# saveRDS(invivo, "~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Leukemia_Doublets_excluded/BASIC_PIPELINE_merged_co_mono/co_mono_invivo(cmi)/invivo.rds")

# merge ex vivo data together to remove any supporting cells from co-culture----
co_mono<- merge(co, y = mono, add.cell.ids = c("co", "mono"), project = "co_mono")
# saveRDS(co_mono, "~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Doublets_excluded/rds.files/co_mono.rds")

# run standard pipeline 
co_mono <- NormalizeData(co_mono, normalization.method = "LogNormalize", scale.factor = 10000)
co_mono <- FindVariableFeatures(co_mono, selection.method = "vst")
co_mono <- ScaleData(co_mono)
co_mono <- RunPCA(co_mono, features = VariableFeatures(object = co_mono))
co_mono <- FindNeighbors(co_mono, dims = 1:20)
co_mono <- FindClusters(co_mono, resolution = 0.5)
set.seed(42)
co_mono <- RunUMAP(co_mono, dims = 1:20)

# remove supporting cells
FeaturePlot(co_mono, features = c("CD19", "CD22", "CD7", "CD5", "CDH5", "ENG", "PECAM1", "MCAM", "CXCL12", "PDGFRB"))
DimPlot(co_mono, reduction = "umap", label = TRUE)
co_mono_Leu <- subset(x = co_mono, idents = c(4, 8, 13, 3, 10, 12), invert = TRUE)

# run  pipeline for  leukemic cells
co_mono_Leu <- NormalizeData(co_mono_Leu)
co_mono_Leu <- FindVariableFeatures(co_mono_Leu, selection.method = "vst")
co_mono_Leu <- ScaleData(co_mono_Leu)
co_mono_Leu <- RunPCA(co_mono_Leu, features = VariableFeatures(object = co_mono_Leu))
co_mono_Leu <- FindNeighbors(co_mono_Leu, dims = 1:20)
co_mono_Leu <- FindClusters(co_mono_Leu, resolution = 0.3)
set.seed(42)
co_mono_Leu <- RunUMAP(co_mono_Leu, dims = 1:20)
DimPlot(co_mono_Leu, reduction = "umap", group.by = "orig.ident", label = TRUE)
# saveRDS(co_mono_Leu, "~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Doublets_excluded/rds.files/co_mono_Leu.rds")

# subset datasets based on subtype ----
# this step was performed to look into subtype-based signatures, and to account for subtype variability downstream
Idents(co_mono_Leu) <- "orig.ident"
co_mono_Leu_HLF <- subset(x = co_mono_Leu, idents = c("CO1", "CO2", "MONO1", "MONO2"))
co_mono_Leu_PBX1 <- subset(x = co_mono_Leu, idents = c("CO3", "CO4", "MONO3", "MONO4"))
co_mono_Leu_T<- subset(x = co_mono_Leu, idents = c("CO5", "CO6", "MONO5", "MONO6"))

# saveRDS(co_mono_Leu_HLF, "~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Leukemia_Doublets_excluded/BASIC_PIPELINE_merged_co_mono/co_mono_subtype/co_mono_Leu_HLF.rds")
# saveRDS(co_mono_Leu_PBX1, "~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Leukemia_Doublets_excluded/BASIC_PIPELINE_merged_co_mono/co_mono_subtype/co_mono_Leu_PBX1.rds")
# saveRDS(co_mono_Leu_T, "~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Leukemia_Doublets_excluded/BASIC_PIPELINE_merged_co_mono/co_mono_subtype/co_mono_Leu_T.rds")

Idents(invivo)<- "orig.ident"
invivo_HLF <- subset(x = invivo, idents = c("IN_VIVO1", "IN_VIVO2"))
invivo_PBX1 <- subset(x = invivo, idents = c("IN_VIVO3", "IN_VIVO4"))
invivo_T<- subset(x = invivo, idents = c("IN_VIVO5", "IN_VIVO6"))

# saveRDS(invivo_HLF, "~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Analysis/Leukemia_Doublets_excluded/BASIC_PIPELINE_merged_co_mono_invivo/7.co_mono_invivo(cmi)/invivo_d0/invivo_HLF.rds")
# saveRDS(invivo_PBX1, "~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Analysis/Leukemia_Doublets_excluded/BASIC_PIPELINE_merged_co_mono_invivo/7.co_mono_invivo(cmi)/invivo_d0/invivo_PBX1.rds")
# saveRDS(invivo_T, "~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Analysis/Leukemia_Doublets_excluded/BASIC_PIPELINE_merged_co_mono_invivo/7.co_mono_invivo(cmi)/invivo_d0/invivo_T.rds")

# process all new datasets with standard pipeline
seurat_objects_new <- c("co_mono_Leu_HLF","co_mono_Leu_PBX1", "co_mono_Leu_T", "invivo_HLF", "invivo_PBX1", "invivo_T")
for (object_name in seurat_objects_new) {
  obj <- get(object_name)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst")
  all.genes <- rownames(obj)
  obj <- ScaleData(obj, features = all.genes)
  obj <- RunPCA(obj, features = VariableFeatures(object = obj))
  obj <- FindNeighbors(obj, dims = 1:20)
  obj <- FindClusters(obj)
  set.seed(42)
  obj <- RunUMAP(obj, dims = 1:20)
  assign(object_name, obj)}

## merge ex vivo and in vivo 
cmi<- merge(co_mono_Leu_HLF, y = c(co_mono_Leu_PBX1, co_mono_Leu_T, invivo_HLF, invivo_PBX1, invivo_T), 
            add.cell.ids = c("cm_H", "cm_P", "cm_T", "i_H", "i_P", "i_T"), project = "cm_i_cc")

cmi <- NormalizeData(cmi)
cmi <- FindVariableFeatures(cmi)
cmi <- ScaleData(cmi)
cmi <- RunPCA(cmi)
cmi <- FindNeighbors(cmi)
cmi <- FindClusters(cmi, resolution = 0.1)
set.seed(42)
cmi <- RunUMAP(cmi, dims = 1:20)
DimPlot(cmi, reduction = "umap", group.by = 'orig.ident', label = FALSE)

# saveRDS(cmi, "~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Leukemia_Doublets_excluded/BASIC_PIPELINE_merged_co_mono_invivo/7.co_mono_invivo(cmi)/cmi_int_final_used/cmi.rds")

## data integration based on experiment ---- 
# create column defining ex vivo and in vivo(check if identities co, mono, invivo are correct)
cmi$experiment <- NA
cmi$experiment[grepl("^CO", cmi$orig.ident)] <- "exvivo"
cmi$experiment[grepl("^MONO", cmi$orig.ident)] <- "exvivo"
cmi$experiment[grepl("^invivo", cmi$orig.ident)] <- "invivo"

# integration
# split the dataset into a list of two seurat objects based on experiment
cmi.list <- SplitObject(cmi, split.by = c("experiment"))

# normalize and identify variable features for each dataset independently
cmi.list <- lapply(X = cmi.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)})

# select features that are repeatedly variable across datasets for integration run PCA 
features <- SelectIntegrationFeatures(object.list = cmi.list)
cmi.list <- lapply(X = cmi.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)})

# integrate
cmi.anchors <- FindIntegrationAnchors(object.list = cmi.list, anchor.features = features, reduction = "rpca")
cmi_int <- IntegrateData(anchorset = cmi.anchors)
DefaultAssay(cmi_int) <- "integrated"
cmi_int <- ScaleData(cmi_int, verbose = FALSE)
cmi_int <- RunPCA(cmi_int, npcs = 30, verbose = FALSE)
cmi_int <- RunUMAP(cmi_int, reduction = "pca", dims = 1:30)
cmi_int <- FindNeighbors(cmi_int, reduction = "pca", dims = 1:30)
cmi_int <- FindClusters(cmi_int, resolution = 0.2)
DimPlot(cmi_int, reduction = "umap", group.by = c("Phase", "condition")) + ggtitle("cmi_int")

# saveRDS(cmi_int, "~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Leukemia_Doublets_excluded/BASIC_PIPELINE_merged_co_mono_invivo/7.co_mono_invivo(cmi)/cmi_int_final_used/cmi_int.rds")

## downstream analysis on RNA ----
DefaultAssay(cmi_int) <- "RNA"
Idents(cmi_int) <- "condition"

# up- and down-regulated genes comparing mono- and co-culture ----
co_markers <-FindMarkers(cmi_int, ident.1 = "co", ident.2 = "mono", subset.ident = 'condition')
co_markers$geneSymbol <- rownames(co_markers)
# saveRDS(co_markers, "~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Leukemia_Doublets_excluded/BASIC_PIPELINE_merged_co_mono_invivo/7.co_mono_invivo(cmi)/cmi_final_used/co_markers.rds")

# ...dot plot for representative upregulated markers - figures_paper.R

# run gsea on these markers ----

# load refernce files 
Reactome <- fgsea::gmtPathways("~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Reference_gmt_files/c2.cp.reactome.v2023.1.Hs.symbols.gmt")
hallmark <- fgsea::gmtPathways("~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Reference_gmt_files/h.all.v2023.1.Hs.symbols.gmt")
KEGG <- fgsea::gmtPathways("~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Reference_gmt_files/c2.cp.kegg.v2023.1.Hs.symbols.gmt")
GO <- fgsea::gmtPathways("~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Reference_gmt_files/c5.all.v2023.1.Hs.symbols.gmt")

# run the gsea
co_markers <- co_markers %>% arrange(desc(avg_log2FC))
fold_changes <- co_markers$avg_log2FC
names(fold_changes) <- co_markers$geneSymbol
gsea_co <- fgsea(pathways = KEGG, stats = fold_changes, eps = 0.0, minSize=15, maxSize=500)
# saveRDS(gsea_co, "~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Leukemia_Doublets_excluded/BASIC_PIPELINE_merged_co_mono_invivo/7.co_mono_invivo(cmi)/cmi_final_used/gsea.co.rds")

# ...gsea dot plot - figures_paper.R

# emt scoring ----
EMT_genes <- read.gmt("~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTicA/Reference_gmt_files/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.v2023.2.Hs.gmt")

# subset object to ex vivo only, since emt was found upregulated in co vs mono
Idents(cmi_int) <- "condition"
cmi_int_exvivo <- subset(x = cmi_int, idents = c("co", "mono"))
EMT_score_cmi_int <- AddModuleScore(object = cmi_int_exvivo, features = EMT_genes,  name = 'EMT_Features_cmi')
FeaturePlot(EMT_score_cmi_int, features = "EMT_Features_cmi1", cols = c("lightgray", "#BE2641"), min.cutoff = 0, max.cutoff = 0.3) + ggtitle("EMT signature")

# look into emt based on subtype in co-culture
EMT_score_cmi_int_co <- subset(x = EMT_score_cmi_int, idents = "co")

# percentage of cells expressing emt
EMT_score_cmi_int_co$subtype2 <- NA
EMT_score_cmi_int_co$subtype2[grepl("^CO1|^CO2", EMT_score_cmi_int_co$orig.ident)] <- "B-ALL"
EMT_score_cmi_int_co$subtype2[grepl("^CO3|^CO4", EMT_score_cmi_int_co$orig.ident)] <- "B-ALL"
EMT_score_cmi_int_co$subtype2[grepl("^CO5|^CO6", EMT_score_cmi_int_co$orig.ident)] <- "T-ALL"

Idents(EMT_score_cmi_int_co) <- "subtype2"
EMT_score_cmi_int_co_B <- subset(x = EMT_score_cmi_int_co, idents = "B-ALL")
EMT_score_cmi_int_co_T <- subset(x = EMT_score_cmi_int_co, idents = "T-ALL")

positive_cellsB <- which(EMT_score_cmi_int_co_B@meta.data$EMT_Features_cmi > 0)
percentage_positiveB <- length(positive_cellsB) / nrow(EMT_score_cmi_int_co_B@meta.data) * 100

positive_cellsT <- which(EMT_score_cmi_int_co_T@meta.data$EMT_Features_cmi > 0)
percentage_positiveT <- length(positive_cellsT) / nrow(EMT_score_cmi_int_co_T@meta.data) * 100

# ... plot violin plot - figures_paper.R

# saveRDS(cmi_int_exvivo, "~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Leukemia_Doublets_excluded/BASIC_PIPELINE_merged_co_mono_invivo/7.co_mono_invivo(cmi)/cmi_FINAL_USED/cmi_int_exvivo.rds")
# saveRDS(EMT_score_cmi_int_co, "~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Analysis/Leukemia_Doublets_excluded/BASIC_PIPELINE_merged_co_mono_invivo/7.co_mono_invivo(cmi)/cmi_FINAL_USED/EMT_score_cmi_int_co.rds")

# cycling vs non cycling cells ----
Cycling <- WhichCells(cmi_int, expression = S.Score > 0 | G2M.Score > 0)
nonCycling <- WhichCells(cmi_int, expression = S.Score <= 0 & G2M.Score <= 0)
cmi_int$Cycling_prop <- ifelse(colnames(cmi_int) %in% Cycling, "Cycling", "Non-cycling")
DimPlot(cmi_int, reduction = "umap", group.by = "Cycling_prop", cols = c("#edc48e", "#BE2641")) + ggtitle("Cell cycle state")

# look into subtype differences
cmi_int$subtype2 <- NA
cmi_int$subtype2[grepl("^CO1|^CO2|^MONO1|^MONO2|^invivo1|^invivo2", cmi_int$orig.ident)] <- "B-ALL"
cmi_int$subtype2[grepl("^CO3|^CO4|^MONO3|^MONO4|^invivo3|^invivo4", cmi_int$orig.ident)] <- "B-ALL"
cmi_int$subtype2[grepl("^CO5|^CO6|^MONO5|^MONO6|^invivo5|^invivo6", cmi_int$orig.ident)] <- "T-ALL"

Idents(cmi_int) <- "subtype2"
cmi_int_B <- subset(x = cmi_int, idents = "B-ALL")
cmi_int_T <- subset(x = cmi_int, idents = "T-ALL")

Idents(cmi_int_B) <- "condition"

# ... stacked plots -  figures_paper.R

