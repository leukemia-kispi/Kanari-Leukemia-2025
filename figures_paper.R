# Code for figures in manuscript - Magdalini Kanari
# RDS objects created in supporting_cells_paper.R & leukemia_cells_paper.R

# load packages 
library(Seurat)
library(cowplot)
library(patchwork)
library(ggplot2)
library(dplyr)
library(paletteer)

## HUVECs ----
sup_cells_HUVEC <- readRDS("~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Analysis/Supporting_Cells/Combined_MSCHUVEC/Subset_from_MSCHUVEC/HUVEC_condition_comparison_FINAL_USED/sup_cells_HUVEC.rds")

# create columns for samples
sup_cells_HUVEC$condition <- NA
sup_cells_HUVEC$condition[grepl("^HUVEC", sup_cells_HUVEC$orig.ident)] <- "HUVEC only"
sup_cells_HUVEC$condition[grepl("^MSC_HUVEC", sup_cells_HUVEC$orig.ident)] <- "HUVEC vascular"

sup_cells_HUVEC$condition2 <- NA
sup_cells_HUVEC$condition2[grepl("^CO1", sup_cells_HUVEC$orig.ident)] <- "HUVEC leukemic - CO1"
sup_cells_HUVEC$condition2[grepl("^CO2", sup_cells_HUVEC$orig.ident)] <- "HUVEC leukemic - CO2"
sup_cells_HUVEC$condition2[grepl("^CO3", sup_cells_HUVEC$orig.ident)] <- "HUVEC leukemic - CO3"
sup_cells_HUVEC$condition2[grepl("^CO4", sup_cells_HUVEC$orig.ident)] <- "HUVEC leukemic - CO4"
sup_cells_HUVEC$condition2[grepl("^CO5", sup_cells_HUVEC$orig.ident)] <- "HUVEC leukemic - CO5"
sup_cells_HUVEC$condition2[grepl("^CO6", sup_cells_HUVEC$orig.ident)] <- "HUVEC leukemic - CO6"

# umap projections based on the sample
palette <- paletteer_c("grDevices::ag_GrnYl", 6)
DimPlot(object = subset(sup_cells_HUVEC,  orig.ident %in% c("HUVEC", "MSC_HUVEC")), 
        group.by = "condition", pt.size = 0.5,  cols = c("#009380", "#4FB978")) +ggtitle("Control HUVEC") +
  theme(plot.title = element_text(size = 14, face = "bold", color = "black"))
DimPlot(object = subset(sup_cells_HUVEC,  orig.ident %in% c("CO1", "CO2", "CO3", "CO4", "CO5", "CO6")), 
        shuffle = TRUE, group.by = "condition2", pt.size = 0.5) + ggtitle("Leukemic HUVEC") + 
  theme(plot.title = element_text(size = 14, face = "bold", color = "black")) +
  scale_color_manual(values = as.character(palette))

# set new identities for umap and dotplot 
Idents(sup_cells_HUVEC) <- "RNA_snn_res.0.2"
new.cluster.ids_HUVEC <- c("vascular", "mesenchymal-like", "vascular", "pre-vascular")
names(new.cluster.ids_HUVEC) <- levels(sup_cells_HUVEC)
sup_cells_HUVEC <- RenameIdents(sup_cells_HUVEC, new.cluster.ids_HUVEC)
h <-DimPlot(sup_cells_HUVEC, reduction = "umap", label.size = 5, pt.size = 0.5, cols = c("#009380",  "#196170", "#4FB978"))+
  labs(title = "Identified HUVEC populations") +  theme(plot.title = element_text(size = 16, face = "bold", color = "black"))
LabelClusters(h, id = "ident",  fontface = "bold", color = "gray5", size = 5, box = TRUE, fill = "white")

dp <- DotPlot(sup_cells_HUVEC, features = c("CDH5", "PECAM1", "VWF", "ICAM1", "LYVE1", "EPHB2","NPR1", "NR2F2", "KRT7",
                                            "RGS5", "ANKRD1", "TYMS", "ODC1", "CD34", "MARCKSL1", "CHI3L1", "DEPP1", "SOX18"), 
              dot.scale = 4, cols = c("gray90", "#009380"),col.min = 0, col.max = 2) + 
  theme(axis.text.x = element_text(angle = 20, hjust=1, size =12)) +
  labs (title = "Identification genes for HUVECs")
dp + coord_flip() + theme(plot.title = element_text(hjust=0.3), axis.text.y = element_text(size = 8))

# huvec composition based on condition - stacked plot
sup_cells_HUVEC$condition <- NA
sup_cells_HUVEC$condition[grepl("^CO", sup_cells_HUVEC$orig.ident)] <- "HUVEC leukemic"
sup_cells_HUVEC$condition[grepl("^MSC_HUVEC", sup_cells_HUVEC$orig.ident)] <- "HUVEC vascular"
sup_cells_HUVEC$condition[grepl("^HUVEC", sup_cells_HUVEC$orig.ident)] <- "HUVEC only"

id_composition <- as.data.frame(table(sup_cells_HUVEC$condition, Idents(sup_cells_HUVEC)))
id_composition$Percentage <- ave(id_composition$Freq, id_composition$Var1, FUN = function(x) x/sum(x) * 100)
id_composition$Var2 <- factor(id_composition$Var2, levels = c("vascular", "mesenchymal-like", "pre-vascular"))

ggplot(id_composition, aes(x = factor(Var1, levels = c("HUVEC only", "HUVEC vascular", "HUVEC leukemic")), 
                           y = Percentage, fill = Var2)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("#009380",  "#196170", "#4FB978")) +
  labs(x = "", y = "Percentage", title = "Identified HUVEC populations per sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 2, hjust = 0.5, size = 9, face = "bold", color = "black"),
        axis.text.y = element_text(size = 9, face = "bold", color = "black"),
        axis.title.y = element_text(size = 12, face = "bold", color = "black"),
        plot.title = element_text(hjust = 0, size = 16, face = "bold", color = "black"))

# load marker files (see find markers in supporting_cells_paper.R)
sup_cells_HUVEC_3 <- readRDS ("~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Analysis/Supporting_Cells/Combined_MSCHUVEC/Subset_from_MSCHUVEC/HUVEC_condition_comparison_FINAL_USED/sup_cells_HUVEC_3.rds")
Reactome <- fgsea::gmtPathways("~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Reference_gmt_files/c2.cp.reactome.v2023.1.Hs.symbols.gmt")

sup_cells_HUVEC_3$geneSymbol <- rownames(sup_cells_HUVEC_3)
sup_cells_HUVEC_3 <- sup_cells_HUVEC_3 %>% arrange(desc(avg_log2FC))
fold_changes <- sup_cells_HUVEC_3$avg_log2FC
names(fold_changes) <- sup_cells_HUVEC_3$geneSymbol
gsea3 <- fgsea(pathways = Reactome, stats = fold_changes, eps = 0.0, minSize=15, maxSize=500)

path_3 <- gsea3[order(-gsea3$NES), ][1:5, ]
path_3$source <- "HUVEC leukemic to HUVEC vascular"

path_3$pathway <- gsub("_", " ", path_3$pathway)
path_3$pathway <- str_wrap(path_3$pathway, width = 40)

ggplot(path_3, aes(x = pathway, y = NES, fill = -log10(padj))) +
  geom_segment(aes(xend = pathway, yend = 0), color = "#009380", linewidth = 14, alpha = 0.3) +
  geom_text(aes(label = pathway, y = 0.1), size = 3.5, color = "black", angle = 0, hjust = 0) +
  geom_point(size = 12, shape = 21, color = "#009380") +  # Use shape 21 for filled circles
  scale_fill_gradient(low = "whitesmoke", high = "#009380") +  # Set the gradient from white to royalblue
  labs(x = "Pathway", y = "Normalized Enrichment Score") +
  theme_minimal() +
  ggtitle ("HUVEC leukemic to HUVEC vascular")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1, color = "black", face = "bold", size = 12),
        axis.title.x = element_text(size = 12, face = "bold", color = "black"),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "black")) +   coord_flip()

## MSCs ----
sup_cells_MSC_CDH5negCD19neg <- readRDS("~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Supporting_Cells/Combined_MSCHUVEC/Subset_from_MSCHUVEC/CDH5negCD19neg_MSCs_FINAL_USED/sup_cells_MSC_CDH5negCD19neg.rds")

# create columns for samples
sup_cells_MSC_CDH5negCD19neg$condition <- NA
sup_cells_MSC_CDH5negCD19neg$condition[grepl("^MSC", sup_cells_MSC_CDH5negCD19neg$orig.ident)] <- "MSC only"
sup_cells_MSC_CDH5negCD19neg$condition[grepl("^MSC_HUVEC", sup_cells_MSC_CDH5negCD19neg$orig.ident)] <- "MSC vascular"
sup_cells_MSC_CDH5negCD19neg$condition[grepl("^CO", sup_cells_MSC_CDH5negCD19neg$orig.ident)] <- "MSC leukemic"

sup_cells_MSC_CDH5negCD19neg$condition2 <- NA
sup_cells_MSC_CDH5negCD19neg$condition2[grepl("^CO1", sup_cells_MSC_CDH5negCD19neg$orig.ident)] <- "MSC leukemic - CO1"
sup_cells_MSC_CDH5negCD19neg$condition2[grepl("^CO2", sup_cells_MSC_CDH5negCD19neg$orig.ident)] <- "MSC leukemic - CO2"
sup_cells_MSC_CDH5negCD19neg$condition2[grepl("^CO3", sup_cells_MSC_CDH5negCD19neg$orig.ident)] <- "MSC leukemic - CO3"
sup_cells_MSC_CDH5negCD19neg$condition2[grepl("^CO4", sup_cells_MSC_CDH5negCD19neg$orig.ident)] <- "MSC leukemic - CO4"
sup_cells_MSC_CDH5negCD19neg$condition2[grepl("^CO5", sup_cells_MSC_CDH5negCD19neg$orig.ident)] <- "MSC leukemic - CO5"
sup_cells_MSC_CDH5negCD19neg$condition2[grepl("^CO6", sup_cells_MSC_CDH5negCD19neg$orig.ident)] <- "MSC leukemic - CO6"

# umap projections based on the sample
palette <- paletteer_c("grDevices::YlGnBu", 6)
DimPlot(object = subset(sup_cells_MSC_CDH5negCD19neg,  orig.ident %in% c("MSC","MSC_HUVEC")), 
        group.by = "condition", shuffle = TRUE, cols = c("#0098AF", "#01357B")) + ggtitle("Control MSC") + 
  theme(plot.title = element_text(size = 14, face = "bold", color = "black"))
DimPlot(object = subset(sup_cells_MSC_CDH5negCD19neg,  orig.ident %in% c("CO1", "CO2", "CO3", "CO4", "CO5", "CO6")), 
        group.by = "condition2", shuffle = TRUE) + ggtitle("Leukemic MSC") + 
  scale_color_manual(values = as.character(palette)) +
  theme(plot.title = element_text(size = 14, face = "bold", color = "black"))

# set new identities
Idents(sup_cells_MSC_CDH5negCD19neg)<-"RNA_snn_res.0.5"
new_ids <- c("pre-fibro", "pre-osteo", "pre-fibro-chondro", "osteo-CAR", "pre-chondro", "smooth-muscle", "perivascular", "osteo-CAR", "pre-fibro-chondro",
             "pre-adipo", "smooth-muscle", "cycling", "hybrid", "perivascular")
names(new_ids) <- levels(sup_cells_MSC_CDH5negCD19neg)
sup_cells_MSC_CDH5negCD19neg <- RenameIdents(sup_cells_MSC_CDH5negCD19neg, new_ids)

palette <- paletteer_c("grDevices::BluYl", 10)
p <-DimPlot(sup_cells_MSC_CDH5negCD19neg, reduction = "umap", pt.size = 0.2, alpha = 0.4) + 
  ggtitle("Identified MSC populations") +  scale_color_manual(values = as.character(palette)) + 
  theme(plot.title = element_text(size = 17, face = "bold", hjust = 2))
LabelClusters(p, id = "ident",  fontface = "bold", color = "gray5", size = 4, box = TRUE, fill = "white")

# feature plots for marker expression
FeaturePlot(sup_cells_MSC_CDH5negCD19neg, features = "PDGFRA", cols = c("whitesmoke", "#01357B"), min.cutoff = 0, max.cutoff = 2) + 
  ggtitle("PDGFRA expression") +  theme(plot.title = element_text(size = 14))
  FeaturePlot(sup_cells_MSC_CDH5negCD19neg, features = "CXCL12", cols = c("whitesmoke", "#01357B"), min.cutoff = 0, max.cutoff = 2) + 
  ggtitle("CXCL12 expression") +  theme(plot.title = element_text(size = 14))

# circle plot
Idents(sup_cells_MSC_CDH5negCD19neg)<-"RNA_snn_res.0.5"
new_ids <- c("fibro-lineage", "osteo-lineage", "fibro-chondro-lineage", "osteo-lineage", "chondro-lineage", "smooth-muscle", "perivascular", "osteo-lineage", "fibro-chondro-lineage",
             "adipo-lineage", "smooth-muscle", "cycling", "hybrid", "perivascular")
names(new_ids) <- levels(sup_cells_MSC_CDH5negCD19neg)
sup_cells_MSC_CDH5negCD19neg <- RenameIdents(sup_cells_MSC_CDH5negCD19neg, new_ids)

id_composition <- as.data.frame(table(sup_cells_MSC_CDH5negCD19neg$condition, Idents(sup_cells_MSC_CDH5negCD19neg)))
id_composition$Percentage <- ave(id_composition$Freq, id_composition$Var1, FUN = function(x) x/sum(x) * 100)
id_composition$Var2 <- factor(id_composition$Var2, levels = c("fibro-lineage", "osteo-lineage", "fibro-chondro-lineage", "chondro-lineage", "smooth-muscle", "perivascular",
                                                              "adipo-lineage", "cycling", "hybrid"))
id_composition_filtered <- subset(id_composition, Var1 == "MSC only")

ggplot(id_composition_filtered, aes(x = 1, y = Percentage, fill = Var2)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  scale_fill_manual(values = palette) +
  labs(x = NULL, y = NULL, title = "MSC leukemic") +
  theme_void() +
  coord_polar(theta = "y") +
  guides(fill = guide_legend(reverse = TRUE)) +  # Reverse legend order if needed
  theme(axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5, vjust = -2, size = 24, face = "bold"))

# gsea
sup_cells_MSC_CDH5negCD19neg_markers_all <- readRDS("~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Analysis/Supporting_Cells/Combined_MSCHUVEC/Subset_from_MSCHUVEC/CDH5negCD19neg_MSCs_FINAL_USED/sup_cells_MSC_CDH5negCD19neg_markers_all.rds")
sup_cells_MSC_CDH5negCD19neg_markers_all <- sup_cells_MSC_CDH5negCD19neg_markers_all %>% arrange(desc(avg_log2FC))
fold_changes <- sup_cells_MSC_CDH5negCD19neg_markers_all$avg_log2FC
sup_cells_MSC_CDH5negCD19neg_markers_all$geneSymbol <-rownames(sup_cells_MSC_CDH5negCD19neg_markers_all)
names(fold_changes) <- sup_cells_MSC_CDH5negCD19neg_markers_all$geneSymbol

GO <- fgsea::gmtPathways("~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Reference_gmt_files/c5.all.v2023.1.Hs.symbols.gmt")
gsea_all <- fgsea(pathways = GO, stats = fold_changes, eps = 0.0, minSize=15, maxSize=500)

### Extract the pathways 
top_pathways1 <- gsea_all[order(-gsea_all$NES), ][1:10, ]
top_pathways1$pathway <- gsub("_", " ", top_pathways1$pathway)
# top_pathways1$pathway <- str_wrap(top_pathways1$pathway, width = 55)

ggplot(top_pathways1, aes(x = reorder(pathway, NES), y = NES, fill = -log10(padj))) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.25, binwidth = 0.2) +
  labs(title = "Leukemic & vascular MSC to MSC only") +
  theme_minimal() + theme(plot.title = element_text(size = 12, hjust = 0)) +
  scale_fill_gradient(low = "whitesmoke", high = "#01357B") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0, color = "black", face = "bold"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 1.5, color = "black", face = "bold", size = 14),
        axis.title.x = element_text(color = "black", face = "bold")) + coord_flip()

## ALLs ----
cmi <- readRDS("~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Analysis/Leukemia_Doublets_excluded/BASIC_PIPELINE_merged_co_mono_invivo/7.co_mono_invivo(cmi)/cmi_FINAL_USED/cmi.rds")
cmi_int <- readRDS("~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Analysis/Leukemia_Doublets_excluded/BASIC_PIPELINE_merged_co_mono_invivo/7.co_mono_invivo(cmi)/cmi_FINAL_USED/cmi_int.rds")

# create columns for samples 
cmi$condition <- NA
cmi$condition[grepl("^CO", cmi_int$orig.ident)] <- "co"
cmi$condition[grepl("^MONO", cmi_int$orig.ident)] <- "mono"
cmi$condition[grepl("^invivo", cmi_int$orig.ident)] <- "invivo"

cmi_int$condition <- NA
cmi_int$condition[grepl("^CO", cmi$orig.ident)] <- "co"
cmi_int$condition[grepl("^MONO", cmi$orig.ident)] <- "mono"
cmi_int$condition[grepl("^invivo", cmi$orig.ident)] <- "invivo"

cmi_int$subtype <- NA
cmi_int$subtype[grepl("^CO1|^CO2|^MONO1|^MONO2|^invivo1|^invivo2", cmi_int$orig.ident)] <- "TCF3-HLF"
cmi_int$subtype[grepl("^CO3|^CO4|^MONO3|^MONO4|^invivo3|^invivo4", cmi_int$orig.ident)] <- "TCF3-PBX1"
cmi_int$subtype[grepl("^CO5|^CO6|^MONO5|^MONO6|^invivo5|^invivo6", cmi_int$orig.ident)] <- "T-ALL"

DimPlot(cmi, group.by = "condition", cols =c("#BE2641", "#cc6e77", "#dab5ad"), pt.size = 0.12) + ggtitle("Before integration") + theme(element_text(size = 16, face = "bold"))
DimPlot(cmi_int, group.by = "condition",cols =c("#BE2641", "#cc6e77", "#dab5ad"), pt.size = 0.12) + ggtitle("Condition") + theme(element_text(size = 16, face = "bold")) 
DimPlot(cmi_int, group.by = "subtype", cols =c("#BE2641", "#cc6e77", "#dab5ad"), pt.size = 0.12) + ggtitle ("Subtype") + theme(element_text(size = 16, face = "bold"))

# dot plots for representative upregulated markers 
co_markers <- readRDS("~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Analysis/Leukemia_Doublets_excluded/BASIC_PIPELINE_merged_co_mono_invivo/7.co_mono_invivo(cmi)/cmi_FINAL_USED/co_markers.rds")

DefaultAssay(cmi_int) <- "RNA"
invivo_markers_B <- c("VPREB1", "TCF4", "MEF2C", "PAX5", "EBF1", "MZB1", "CD52", "MGMT", "AUTS2")
invivo_markers_T <- c("SIT1", "SIVA1", "PTEN", "MZB1", "PDE3B", "CDKN2A", "MGMT", "AUTS2", "FYN")
invivo_markers_B_data <- subset(co_markers, geneSymbol %in% invivo_markers_B)
invivo_markers_T_data <- subset(co_markers, geneSymbol %in% invivo_markers_T)

# divide into B and T based on pdx
numbers <- gsub("\\D", "", cmi_int$orig.ident)
cmi_int$pdx <- paste0("pdx", numbers)

Idents(cmi_int) <- "pdx"
cmi_int_B <- subset(cmi_int, idents = c("pdx1", "pdx2", "pdx3", "pdx4"))
cmi_int_T <- subset(cmi_int, idents = c("pdx5", "pdx6"))

Idents(cmi_int_B) <- "condition"
DotPlot(cmi_int_BT, features = invivo_markers_data$geneSymbol, dot.min = 0, dot.scale = 5, scale.by = "size",
        cols = c("whitesmoke", "#BE2641")) +
  theme(axis.text.x = element_text(hjust = 0.5, size = 9, face = "bold"),
        text = element_text(size = 9),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.3),
        axis.text.y = element_text(size = 8))+
  ggtitle("Leukemia associated genes - B-ALL") +  coord_flip()

Idents(cmi_int_T) <- "condition"
DotPlot(cmi_int_T, features = invivo_markers_data$geneSymbol, dot.min = 0, dot.scale = 5, scale.by = "size",
        cols = c("whitesmoke", "#BE2641")) +
  theme(axis.text.x = element_text(hjust = 0.5, size = 9, face = "bold"),
        text = element_text(size = 9),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.3),
        axis.text.y = element_text(size = 8))+
  ggtitle("Leukemia associated genes - T-ALL") +  coord_flip()

# GSEA dot plot
gsea_co <- readRDS("~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Analysis/Leukemia_Doublets_excluded/BASIC_PIPELINE_merged_co_mono_invivo/7.co_mono_invivo(cmi)/cmi_FINAL_USED/gsea_co.rds")

top_upregulated_pathways <- gsea_co[order(-gsea_co$NES), ][1:10 ]
top_upregulated_pathways$pathway <- gsub("^HALLMARK_", "", top_upregulated_pathways$pathway)

ggplot(top_upregulated_pathways, aes(x = NES, y = reorder(pathway, NES), fill = -log10(padj))) +
  geom_point(shape = 21, color = "black", stroke = 0.3, size = 3) + 
  scale_size_continuous(range = c(2, 4), breaks = seq(0, max(-log10(top_upregulated_pathways$padj)), length.out = 4)) +
  scale_fill_gradientn(colors = c("whitesmoke", "#BE2641")) +
  labs(title = "Upregulated hallmark pathways _ Co vs mono",
       x = "NES", 
       y = "", ) + theme_minimal()+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 1),
        text = element_text(size = 10, color = "black"),
        axis.text =  element_text(size = 10, color = "black"),
        axis.title.x =  element_text(size = 10, color = "black", face = "bold"))

# emt scoring based on subtype - vln plot
EMT_score_cmi_int_co <- readRDS("~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Analysis/Leukemia_Doublets_excluded/BASIC_PIPELINE_merged_co_mono_invivo/7.co_mono_invivo(cmi)/cmi_FINAL_USED/EMT_score_cmi_int_co.rds")

EMT_score_cmi_int_co$subtype2 <- NA
EMT_score_cmi_int_co$subtype2[grepl("^CO1|^CO2", EMT_score_cmi_int_co$orig.ident)] <- "B-ALL"
EMT_score_cmi_int_co$subtype2[grepl("^CO3|^CO4", EMT_score_cmi_int_co$orig.ident)] <- "B-ALL"
EMT_score_cmi_int_co$subtype2[grepl("^CO5|^CO6", EMT_score_cmi_int_co$orig.ident)] <- "T-ALL"

VlnPlot(EMT_score_cmi_int_co, group.by = "subtype2", features = "EMT_Features_cmi1", pt.size = 0) +
  #geom_jitter(size = 0.4, position = position_jitter(width = 0.4, height = 0), alpha = 0.2, color = "black") +  
  scale_fill_manual(values = c("B-ALL" = "#BE2641", "T-ALL" = "#dab5ad")) +  theme_minimal() +
  ggtitle("EMT Signature - Coculture") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text =  element_text(size = 10, color = "black", face = "bold")) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 1)
  
# cell cycle - calculate for B and T , show graph with co and in vivo samples
Idents(cmi_int_B) <- "condition"
cmi_int_B <- subset(x = cmi_int_B, idents = c("co", "invivo"))

id_composition <- as.data.frame(table(Idents(cmi_int_B), cmi_int_B$Cycling_prop))
id_composition$Percentage <- ave(id_composition$Freq, id_composition$Var1, FUN = function(x) x/sum(x) * 100)

ggplot(id_composition, aes(x = Var1, y = Percentage, fill = Var2)) +
  geom_bar(stat = "identity", color = "black", alpha = 1) +
  scale_fill_manual(values = c("#dab5ad", "#BE2641")) +
  labs(x = "Condition", y = "Percentage [%]", title = "Cell cycle B-ALL", fill = "") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10, color = "black", face = "bold"),
        axis.title = element_text(size = 10, color = "black", face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

# Repeat all above for cell cycle with cmi_int_T sample
Idents(cmi_int_T) <- "condition"
cmi_int_T <- subset(x = cmi_int_T, idents = c("co", "invivo"))

id_composition <- as.data.frame(table(Idents(cmi_int_T), cmi_int_T$Cycling_prop))
id_composition$Percentage <- ave(id_composition$Freq, id_composition$Var1, FUN = function(x) x/sum(x) * 100)

ggplot(id_composition, aes(x = Var1, y = Percentage, fill = Var2)) +
  geom_bar(stat = "identity", color = "black", alpha = 1) +
  scale_fill_manual(values = c("#dab5ad", "#BE2641")) +
  labs(x = "Condition", y = "Percentage [%]", title = "Cell cycle B-ALL", fill = "") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10, color = "black", face = "bold"),
        axis.title = element_text(size = 10, color = "black", face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

# enrichment plot for EMT comparing non cycling vs cycling cells
Idents(cmi_int_exvivo_co) <- "Phase"
cmi_np <-FindMarkers(cmi_int_exvivo_co, ident.1 = c("G1"), ident.2 = c("G2M", "S"), subset.ident = 'Phase')
cmi_np$geneSymbol <- rownames(cmi_np)
cmi_np <- cmi_np %>% arrange(desc(avg_log2FC))
fold_changes <- cmi_np$avg_log2FC
names(fold_changes) <- cmi_np$geneSymbol

# run the GSEA
hallmark <- fgsea::gmtPathways("~/kispi/Research/13_Sequencing_analysis/MagdaliniKanari/Magda_scRNA_ECTICA/Reference_gmt_files/h.all.v2023.1.Hs.symbols.gmt")
gsea_cmi_np <- fgsea(pathways = hallmark, stats = fold_changes, eps = 0.0, minSize=15, maxSize=500)

plotEnrichment(hallmark[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]], fold_changes) + 
  labs(title="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
       subtitle = "Upregulation in non cycling vs cycling cells") +
  geom_line(color = "#BE2641", size = 1.2) +        
  geom_area(fill = "whitesmoke", alpha = 0.4) +        
  theme_minimal() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 1),
        plot.subtitle = element_text(size = 10, hjust = -0.32),
        text = element_text(color = "black"))

