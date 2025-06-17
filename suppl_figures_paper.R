# Code for supplementary figures in manuscript - Magdalini Kanari
# RDS objects created in supporting_cells_paper.R & leukemia_cells_paper.R

## load packages----
library(Seurat)
library(cowplot)
library(patchwork)
library(ggplot2)
library(dplyr)
library(CellChat)
library(qusage)

## suppl. figure 2 ----
sup_cells_MSC_CDH5negCD19neg <- readRDS("~/sup_cells_MSC_CDH5negCD19neg.rds")

# marker dot plot
dp <- DotPlot(sup_cells_MSC_CDH5negCD19neg, features = c("PDGFRA", "PDGFRB", "CXCL12", "LEPR", "DCN", "FAP", "CCN2", 
                                                         "SOX5", "SOX9", "BMP2", "S100A6", "BGLAP", "TAGLN", "TAGLN2",
                                                         "ACTA2", "CALD1", "COL4A2", "MCAM", "ANGPT1", "ACAN",
                                                         "VCAM1", "HOXC6", "HOXC8", "MKI67", "CXCR4"), 
              dot.scale = 8, cols = c("whitesmoke","#426AAC"),col.min = 0, col.max = 2)
dp + coord_flip() + theme(axis.text.x = element_text(angle = 30, hjust=1))

# stacked plot for MSC type per original sample 
sup_cells_MSC_CDH5negCD19neg$condition <- NA
sup_cells_MSC_CDH5negCD19neg$condition[grepl("^CO", sup_cells_MSC_CDH5negCD19neg$orig.ident)] <- "MSC leukemic"
sup_cells_MSC_CDH5negCD19neg$condition[grepl("^MSC", sup_cells_MSC_CDH5negCD19neg$orig.ident)] <- "MSC only"
sup_cells_MSC_CDH5negCD19neg$condition[grepl("^MSC_HUVEC", sup_cells_MSC_CDH5negCD19neg$orig.ident)] <- "MSC vascular"
sup_cells_MSC_CDH5negCD19neg$condition[grepl("^HUVEC", sup_cells_MSC_CDH5negCD19neg$orig.ident)] <- "HUVEC"

id_composition <- as.data.frame(table(sup_cells_MSC_CDH5negCD19neg$condition, Idents(sup_cells_MSC_CDH5negCD19neg)))
id_composition$Percentage <- ave(id_composition$Freq, id_composition$Var1, FUN = function(x) x/sum(x) * 100)
id_composition <- id_composition[id_composition$Var1 != "HUVEC", ]
id_composition$Var1 <- factor(id_composition$Var1, levels = c("MSC only", "MSC vascular", "MSC leukemic"))
id_composition$Var2 <- factor(id_composition$Var2, levels = c("pre-fibro", "pre-osteo", "pre-fibro-chondro", "osteo-CAR", "pre-chondro", 
                                                              "smooth-muscle", "perivascular", "osteo-CAR", "pre-fibro-chondro",
                                                              "pre-adipo", "smooth-muscle", "cycling", "CXCR4", "perivascular"))
ggplot(id_composition, aes(x = Var1, y = Percentage, fill = Var2)) +
  geom_bar(stat = "identity", color = "gray80") +
  scale_fill_manual(values = c("#5a91eb", "#798cd9", "#8e87c8", "#9f82b7", "#ac7da6", "#b87895", "#c17384", "#ca6d73", "#d16863", "#d76252")) +
  labs(x = "Sample", y = "Percentage", title = "Identified populations per sample - MSC") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1)) + coord_flip()

# gsea 
Idents(sup_cells_MSC_CDH5negCD19neg)<-"orig.ident"
sup_cells_MSC_CDH5negCD19neg_markers1 <-FindMarkers(sup_cells_MSC_CDH5negCD19neg, ident.1 = c("MSC_HUVEC"), 
                                                        ident.2 = "MSC" , subset.ident = 'orig.ident')
sup_cells_MSC_CDH5negCD19neg_markers2 <- FindMarkers(sup_cells_MSC_CDH5negCD19neg, ident.1 = c("CO1", "CO2", "CO3", "CO4", "CO5", "CO6"), 
                                                         ident.2 = "MSC" , subset.ident = 'orig.ident')
# gsea
sup_cells_MSC_CDH5negCD19neg_markers1 <- sup_cells_MSC_CDH5negCD19neg_markers1 %>% arrange(desc(avg_log2FC))
fold_changes <- sup_cells_MSC_CDH5negCD19neg_markers1$avg_log2FC
sup_cells_MSC_CDH5negCD19neg_markers1$geneSymbol <-rownames(sup_cells_MSC_CDH5negCD19neg_markers1)
names(fold_changes) <- sup_cells_MSC_CDH5negCD19neg_markers1$geneSymbol
gsea1 <- fgsea(pathways = Reactome, stats = fold_changes, eps = 0.0, minSize=15, maxSize=500)

sup_cells_MSC_CDH5negCD19neg_markers2 <- sup_cells_MSC_CDH5negCD19neg_markers2 %>% arrange(desc(avg_log2FC))
fold_changes <- sup_cells_MSC_CDH5negCD19neg_markers2$avg_log2FC
sup_cells_MSC_CDH5negCD19neg_markers2$geneSymbol <-rownames(sup_cells_MSC_CDH5negCD19neg_markers2)
names(fold_changes) <- sup_cells_MSC_CDH5negCD19neg_markers2$geneSymbol
gsea2 <- fgsea(pathways = Reactome, stats = fold_changes, eps = 0.0, minSize=15, maxSize=500)

# dot plots
top_pathways1 <- gsea1[order(-gsea1$NES), ][1:10, ]
top_pathways2 <- gsea2[order(-gsea2$NES), ][1:10, ]

top_pathways1$pathway <- gsub("_", " ", top_pathways1$pathway)
top_pathways2$pathway <- gsub("_", " ", top_pathways2$pathway)
top_pathways1$pathway <- str_wrap(top_pathways1$pathway, width = 30)
top_pathways2$pathway <- str_wrap(top_pathways2$pathway, width = 30)

# repeat top_pathways1 and top_pathways2
ggplot(top_pathways1, aes(x = reorder(pathway, NES), y = NES, fill = -log10(padj))) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.15, binwidth = 0.2) +
  labs(title = "Leukemic MSC to Vascular MSC") +
  theme_minimal() + theme(plot.title = element_text(size = 12, hjust = 0)) +
  scale_fill_gradient(low = "#9ed0ff", high = "#20329e") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0)) + coord_flip()

# huvecs
sup_cells_HUVEC <- readRDS("~/sup_cells_HUVEC.rds")
sup_cells_HUVEC_1 <- readRDS ("~/sup_cells_HUVEC_1.rds")

Reactome <- fgsea::gmtPathways("~/c2.cp.reactome.v2023.1.Hs.symbols.gmt")

sup_cells_HUVEC_1$geneSymbol <- rownames(sup_cells_HUVEC_1)
sup_cells_HUVEC_1 <- sup_cells_HUVEC_1 %>% arrange(desc(avg_log2FC))
fold_changes <- sup_cells_HUVEC_1$avg_log2FC
names(fold_changes) <- sup_cells_HUVEC_1$geneSymbol
gsea1 <- fgsea(pathways = Reactome, stats = fold_changes, eps = 0.0, minSize=15, maxSize=500)
path_1 <- gsea1[order(-gsea1$NES), ][1:5, ]
path_1$source <- "HUVEC vascular to HUVEC only"

ggplot(path_1, aes(x = pathway, y = NES, fill = -log10(padj))) +
  geom_segment(aes(xend = pathway, yend = 0), color = "#009380", linewidth = 7, alpha = 0.3) +
  geom_text(aes(label = pathway, y = 0.1), size = 3, color = "black", angle = 0, hjust = 0) +
  geom_point(size = 8, shape = 21, color = "#009380") +  # Use shape 21 for filled circles
  scale_fill_gradient(low = "whitesmoke", high = "#009380") +  # Set the gradient from white to royalblue
  labs(x = "Pathway", y = "Normalized Enrichment Score") +
  theme_minimal() +
  ggtitle ("HUVEC vascular to HUVEC only")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "black")) +   coord_flip()

# cellchat
co <- readRDS("~/co.rds")
  
# set identities for 3 cell types 
DimPlot(co, reduction = "umap", label = TRUE)
FeaturePlot(co, features = c("CD19", "CD7", "THY1", "CDH5"))

Idents(co)<-"RNA_snn_res.0.1"
new_ids <- c("MSC", "MSC", "ALL", "ALL", "ALL", "EC", "MSC", "ALL")
names(new_ids) <- levels(co)
co <- RenameIdents(co, new_ids)
DimPlot(co, reduction = "umap", label.size = 5, repel = TRUE)

# subset seurat to create cell type identity 
co_MSC <-subset(x = co, idents = "MSC")
co_ALL <-subset(x = co, idents = "ALL")
co_EC <-subset(x = co, idents = "EC")

co_MSC$celltype <- NA
co_MSC$celltype[grepl("^CO", co$orig.ident)] <- "MSC"

co_ALL$celltype <- NA
co_ALL$celltype[grepl("^CO", co$orig.ident)] <- "ALL"

co_EC$celltype <- NA
co_EC$celltype[grepl("^CO", co$orig.ident)] <- "EC"

# merge the object
co_new<- merge(co_MSC, y = c(co_ALL, co_EC), add.cell.ids = c("MSC", "ALL", "EC"), project = "co_new")
RunUMAP(co_new, dims = 1:20)
DimPlot(co_new, reduction = "umap", group.by = 'celltype', label = FALSE)

# CellChat create
data.input <- co_new[["RNA"]]@data 
labels <- Idents(co_new)
meta <- data.frame(labels = labels, row.names = names(labels)) 
cellChat <- createCellChat(object = co_new, group.by = "celltype", assay = "RNA")

# START CELLCHAT 

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
CellChatDB.use <- subsetDB(CellChatDB)
cellChat@DB <- CellChatDB.use

cellChat <- subsetData(cellChat) 
future::plan("multisession", workers = 4) 
cellChat <- identifyOverExpressedGenes(cellChat)
cellChat <- identifyOverExpressedInteractions(cellChat)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
ptm = Sys.time()
cellChat <- computeCommunProb(cellChat, type = "triMean")
cellChat <- filterCommunication(cellChat, min.cells = 10)
cellChat <- computeCommunProbPathway(cellChat)
cellChat <- aggregateNet(cellChat)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
ptm = Sys.time()
groupSize <- as.numeric(table(cellChat@idents))
par(mfrow = c(1,2), xpd=TRUE)

# visualizations
netVisual_bubble(cellChat, sources.use = 3, targets.use = c(1), remove.isolate = FALSE,  n.colors = 3,  sort.by.source = TRUE) + 
  coord_flip() + ggtitle()
dev.off()
netVisual_chord_gene(cellChat, sources.use = 3, targets.use = 1, lab.cex = 0.5, legend.pos.y = 30, 
                     color.use =  c("#be2641", "#7da47d", "#598fe9"))

## suppl. figure 3 ----
cmi <- readRDS("~/cmi.rds")
cmi_int <- readRDS("~/cmi_int.rds")

DimPlot(cmi, group.by = "condition", cols =c("#BE2641", "#edc48e", "#e18a67")) + ggtitle("Before integration")  

DimPlot(cmi_int, group.by = "condition", cols =c("#BE2641", "#edc48e", "#e18a67")) + ggtitle("After integration")

# emt related genes
co_markers$diffexpressed <- "ns"
co_markers$diffexpressed[co_markers$avg_log2FC > 0.6 & co_markers$new_p_adj < 0.05] <- "up"

top_genes<- co_up_markers %>% filter(pct.1 > 0.2)
top_genes<- top_genes %>% filter(avg_log2FC > 2.5)

extra_genes <- c("CD44", "CD38", "ZEB2", "IGFBP4")
extra_genes <- subset(co_up_markers, geneSymbol %in% extra_genes)

top_genes <-rbind(top_genes, extra_genes)

ggplot(co_up_markers, aes(x = pct.1, y = avg_log2FC)) +
  geom_point(alpha = 15, shape = 21, color = "black", stroke = 0.2, size = 1, fill = "gray90") +
  geom_point(data = top_genes, aes(x = pct.1, y = avg_log2FC), shape = 21, color = "black", fill = "#BE2641", size = 3) +
  geom_text_repel(data = top_genes, aes(label = geneSymbol), hjust = 1, vjust = 0.5, size = 4, segment.color = NA, max.overlaps = Inf) +
  theme_minimal() +
  ylab("Avg_log2FC") +
  xlab("Pct.1") + 
  labs(title = "EMT- or adhesion-related genes upregulated in coculture") +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.4))

# emt regulatory pathways
gsea_co <- readRDS("~/gsea_co.rds")

pathways_of_interest <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", 
                          "HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_HEDGEHOG_SIGNALING", 
                          "HALLMARK_NOTCH_SIGNALING", "HALLMARK_TGF_BETA_SIGNALING")
filtered_results1 <- subset(gsea_co, pathway %in% pathways_of_interest)

ggplot(filtered_results1, aes(x = pathway, y = NES, color = -log10(pval))) +
  geom_segment(aes(xend = pathway, yend = 0), linewidth = 1) +
  geom_point(size = 4, shape = 16) +
  scale_color_gradient(low = "#D9B0B6", high = "#BE2641") +
  labs(title = "EMT-regulatory pathways", x = "Pathway", y = "Normalized Enrichment Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  coord_flip()

# upregulated pathways in co vs mono - subtype based
cmi_int_exvivo$subtype <- NA
cmi_int_exvivo$subtype[grepl("^CO1|^CO2|^MONO1|^MONO2|^CO3|^CO4|^MONO3|^MONO4", cmi_int_exvivo$orig.ident)] <- "B-ALL"
cmi_int_exvivo$subtype[grepl("^CO5|^CO6|^MONO5|^MONO6", cmi_int_exvivo$orig.ident)] <- "T-ALL"

cmi_int_exvivo$condition <- NA
cmi_int_exvivo$condition[grepl("^CO", cmi_int_exvivo$orig.ident)] <- "co"
cmi_int_exvivo$condition[grepl("^MONO", cmi_int_exvivo$orig.ident)] <- "mono"
cmi_int_exvivo$condition[grepl("^invivo", cmi_int_exvivo$orig.ident)] <- "invivo"

# subset based on subtype
Idents(cmi_int_exvivo) <- "subtype"
cmi_int_exvivoB <- subset(cmi_int_exvivo, idents = c("B-ALL"))
cmi_int_exvivoT <- subset(cmi_int_exvivo, idents = c("T-ALL"))

hallmark <- fgsea::gmtPathways("~/h.all.v2023.1.Hs.symbols.gmt")
objects <- c("cmi_int_exvivoB", "cmi_int_exvivoT")

for (object_name in objects) {
  obj <- get(object_name)
 Idents(obj) <- "condition"
  cm_markers <- FindMarkers(obj, ident.1 = "co", ident.2 = "mono", subset.ident = 'condition')
  cm_markers$geneSymbol <- rownames(cm_markers)
  cm_markers <- cm_markers %>% arrange(desc(avg_log2FC))
  fold_changes <- cm_markers$avg_log2FC
  names(fold_changes) <- cm_markers$geneSymbol
  
    # Perform GSEA analysis
  gsea_results <- fgsea(pathways = hallmark, stats = fold_changes, eps = 0.0, minSize = 15, maxSize = 500)
  topPathwaysUp <- gsea_results[order(-gsea_results$NES), ][1:10, ]
  
  # Plot the top 10 upregulated pathways
  plot <- ggplot(topPathwaysUp, aes(x = reorder(pathway, NES), y = NES, fill = -log10(padj))) +
    geom_bar(stat = "identity", width = 0.5) +
    scale_fill_gradient(low = rgb(0.75, 0.15, 0.25, 0.2), high = rgb(0.75, 0.15, 0.25, 0.7)) +
    labs(title = paste("Top 10 upregulated pathways in", object_name, "vs mono in T-ALL"),
         x = "Hallmarks", y = "Normalized Enrichment Score") +
    theme_minimal() + theme(plot.title = element_text(hjust = 1, size = 10)) + coord_flip()}

# umap - cycling, non-cycling 
# define high-cycling vs low-cycling cells
Cycling <- WhichCells(cmi_int, expression = S.Score > 0 | G2M.Score > 0)
nonCycling <- WhichCells(cmi_int, expression = S.Score <= 0 & G2M.Score <= 0)
# add cycling info to metadata column called 'Cycling_prop'
cmi_int$Cycling_prop <- ifelse(colnames(cmi_int) %in% Cycling, "Cycling", "Non-cycling")
DimPlot(cmi_int, reduction = "umap", group.by = "Cycling_prop", cols = c("#edc48e", "#BE2641")) + ggtitle("Cell cycle state")

# selected pathways G1
Idents(cmi_int_exvivo) <- "condition"
cmi_int_exvivo_co <- subset(x = cmi_int_exvivo, idents = "co")

# find markers upregulated in non cycling (G1) cells 
Idents(cmi_int_exvivo_co) <- "Phase"
cmi_np <-FindMarkers(cmi_int_exvivo_co, ident.1 = c("G1"), ident.2 = c("G2M", "S"), subset.ident = 'Phase')
cmi_np$geneSymbol <- rownames(cmi_np)
cmi_np <- cmi_np %>% arrange(desc(avg_log2FC))
fold_changes <- cmi_np$avg_log2FC
names(fold_changes) <- cmi_np$geneSymbol
gsea_cmi_np <- fgsea(pathways = hallmark, stats = fold_changes, eps = 0.0, minSize=15, maxSize=500)

plotEnrichment(hallmark[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]], fold_changes) + 
  labs(title="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", subtitle = "Upregulation in non cycling vs cycling cells")

# choose more interesting downregulated pathways, apart from proliferation
top_upregulated_pathways <- gsea_cmi_np[order(-gsea_cmi_np$NES), ][1:10, ]
top_downregulated_pathways <- c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MTORC1_SIGNALING",
                                "HALLMARK_GLYCOLYSIS", "HALLMARK_APOPTOSIS", "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
                                "HALLMARK_DNA_REPAIR", "HALLMARK_UV_RESPONSE_UP")
top_downregulated_pathways <- gsea_cmi_np[gsea_cmi_np$pathway %in% top_downregulated_pathways, ]
top_pathways <- rbind(top_upregulated_pathways, top_downregulated_pathways)

# Make the dotplot
middle_point <- mean(range(top_pathways$NES))

ggplot(top_pathways, aes(x = NES, y = reorder(pathway, NES), size = -log10(padj), color = NES)) +
  geom_point(alpha = 1) +
  scale_size_continuous(range = c(2, 6), breaks = seq(0, max(-log10(top_pathways$padj)), length.out = 4)) +
  scale_color_gradientn(colors = c("#EA6976", "#F49A6F", "#FDDD9F"),
                        values = scales::rescale(c(min(top_pathways$NES), middle_point, max(top_pathways$NES))),
                        breaks = c(min(top_pathways$NES), middle_point, max(top_pathways$NES)),
                        labels = c("Downregulated", "Not Significant", "Upregulated")) + theme_minimal() + 
  labs(title = "Selected Pathways _ G1 cells", x = "Normalized Enrichment Score (NES)", y = "Pathway")

## suppl. figure 4 ----
cmi_int$subtype2 <- NA
cmi_int$subtype2[grepl("^CO1|^CO2|^MONO1|^MONO2|^invivo1|^invivo2", cmi_int$orig.ident)] <- "B-ALL"
cmi_int$subtype2[grepl("^CO3|^CO4|^MONO3|^MONO4|^invivo3|^invivo4", cmi_int$orig.ident)] <- "B-ALL"
cmi_int$subtype2[grepl("^CO5|^CO6|^MONO5|^MONO6|^invivo5|^invivo6", cmi_int$orig.ident)] <- "T-ALL"

Idents(cmi_int) <- "subtype2"
cmi_int_B <- subset(x = cmi_int, idents = "B-ALL")
cmi_int_T <- subset(x = cmi_int, idents = "T-ALL")

# plot proliferation status per subtype 
objects <- c("cmi_int_B", "cmi_int_T")

for (i in seq_along(objects)) {
  object_name <- objects[i]
  obj <- get(object_name)
  Idents(obj) <- "condition"
  id_composition <- as.data.frame(table(Idents(obj), obj$Cycling_prop))
  id_composition$Percentage <- ave(id_composition$Freq, id_composition$Var1, FUN = function(x) x / sum(x) * 100)
  
  # Set the plot title based on the iteration index
  plot_title <- if (i == 1) {
    "Cell cycle state_B-ALL"} 
  else {
    "Cell cycle state_T-ALL"}
  
  # Create the ggplot object and assign it to a variable
  plot <- ggplot(id_composition, aes(x = Var1, y = Percentage, fill = Var2)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = c("#edc48e", "#BE2641")) +
    labs(x = "Condition", y = "Percentage", title = plot_title, fill = "") + 
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 0, vjust = 2, hjust = 0.5), plot.title = element_text(hjust = 0.5))
  print(plot)}
