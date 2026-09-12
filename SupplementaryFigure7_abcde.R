library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(paletteer)
library(viridis)
library(ggridges)

setwd("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/figure/")

crispr_obj <- readRDS(file = "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/seurat_result/sgRNA/crispr_obj_seurat_sgRNA.rds")
K562<- Read10X("/public/home/wangxinlong/Project/crispr1/data/20250102ABFC20221434-76_K562/sc20241231K562_addsequence/Ab-1_L2_S2401B2401_DBEC_MolsPerCell_MEX")
K562_seurat<-CreateSeuratObject(counts =K562$`Gene Expression`)


##### SupplementaryFigure7 a
VlnPlot(crispr_obj, features = "nCount_sgRNA", group.by = "orig.ident", pt.size = 0) +
  geom_boxplot(width=.2,col="black",fill="white") + 
  scale_fill_manual(values ='#65B8AD' ) + 
  NoLegend()
ggsave("allCell_nCount_sgRNA.pdf", plot = last_plot(), width = 4, height = 5, units = "in", dpi = 300)

##### SupplementaryFigure7 c
VlnPlot(K562_seurat,features=c("nCount_RNA"),ncol=1, pt.size = 0)+
  geom_boxplot(width=.2,col="black",fill="white")+ 
  scale_fill_manual(values ='#65B8AD' ) + 
  NoLegend()
median(K562_seurat$nCount_RNA)
ggsave("Nothresholdapplied_allCell_12468.nCount_RNA.pdf", plot = last_plot(), width = 4, height = 5, units = "in", dpi = 300)





library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(paletteer)
library(viridis)
library(ggridges)



### K562 seurat base analysis
setwd("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/seurat_result/seurat/")

crispr_obj <- readRDS(file = "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/sgRNA_identity/filter_after/crispr_obj.rds")

crispr_obj[["percent.mt"]] <- PercentageFeatureSet(crispr_obj, pattern = "^MT-")

VlnPlot(crispr_obj, features = "nFeature_RNA", pt.size = 0.1)
ggsave("violin_nFeature_RNA_qc_pre.pdf", plot = last_plot(), width = 5, height = 5, units = "in", dpi = 300)
VlnPlot(crispr_obj, features = "nCount_RNA", pt.size = 0.1)
ggsave("violin_nCount_RNA_qc_pre.pdf", plot = last_plot(), width = 5, height = 5, units = "in", dpi = 300)
VlnPlot(crispr_obj, features = "percent.mt", pt.size = 0.1)
ggsave("violin_percent_mt_qc_pre.pdf", plot = last_plot(), width = 5, height = 5, units = "in", dpi = 300)


FeatureScatter(crispr_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
ggsave("pct_counts_mt_qc_pre.pdf", plot = last_plot(), width = 5, height = 4, units = "in", dpi = 300)

FeatureScatter(crispr_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("n_genes_by_counts_qc_pre.pdf", plot = last_plot(), width = 5, height = 4, units = "in", dpi = 300)

crispr_obj <- subset(crispr_obj, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 30)



##### Assign the sgRNA identity to each cell and store it in meta.data, then remove cells without an sgRNA identity
sgTypeCounts_df_filter2 <- read.csv("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/sgRNA_identity/filter_after/cell_sgRNA_identity.csv", row.names = 1,  check.names = FALSE, stringsAsFactors = FALSE)
sgRNA_values <- as.character(sgTypeCounts_df_filter2["sgRNA_identity", ])
barcodes <- colnames(sgTypeCounts_df_filter2)
sgRNA_identity_vector <- rep(NA, ncol(crispr_obj))
matching_indices <- match(barcodes, colnames(crispr_obj))
valid_indices <- matching_indices[!is.na(matching_indices)]
valid_sgRNA_values <- sgRNA_values[!is.na(matching_indices)]
sgRNA_identity_vector[valid_indices] <- valid_sgRNA_values
crispr_obj@meta.data$sgRNA_identity <- sgRNA_identity_vector

crispr_obj
na_cells <- is.na(crispr_obj@meta.data$sgRNA_identity)
crispr_obj <- subset(crispr_obj, cells = colnames(crispr_obj)[!na_cells])
crispr_obj

table(crispr_obj@meta.data[["sgRNA_identity"]])

VlnPlot(crispr_obj, features = "nFeature_RNA", pt.size = 0.1)
ggsave("violin_nFeature_RNA_qc.pdf", plot = last_plot(), width = 5, height = 5, units = "in", dpi = 300)
VlnPlot(crispr_obj, features = "nCount_RNA", pt.size = 0.1)
ggsave("violin_nCount_RNA_qc.pdf", plot = last_plot(), width = 5, height = 5, units = "in", dpi = 300)
VlnPlot(crispr_obj, features = "percent.mt", pt.size = 0.1)
ggsave("violin_percent_mt_qc.pdf", plot = last_plot(), width = 5, height = 5, units = "in", dpi = 300)

FeatureScatter(crispr_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
ggsave("pct_counts_mt_qc.pdf", plot = last_plot(), width = 5, height = 4, units = "in", dpi = 300)

FeatureScatter(crispr_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("n_genes_by_counts_qc.pdf", plot = last_plot(), width = 5, height = 4, units = "in", dpi = 300)


crispr_obj <- NormalizeData(crispr_obj, 
                            normalization.method = "LogNormalize",
                            scale.factor = 1e4)
crispr_obj <- FindVariableFeatures(crispr_obj)

top10 <- head(VariableFeatures(crispr_obj), 10)
plot1 <- VariableFeaturePlot(crispr_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
combined_plot <- plot1 + plot2
ggsave("highly_variable_genes.pdf", plot = last_plot(), width = 12, height = 5, units = "in", dpi = 300)

all.genes <- rownames(crispr_obj)
crispr_obj <- ScaleData(crispr_obj, features = all.genes)

crispr_obj <- RunPCA(crispr_obj)
ElbowPlot(crispr_obj)
ggsave("ElbowPlot.pdf", plot = last_plot(), width = 4, height = 3, units = "in", dpi = 300)

crispr_obj <- FindNeighbors(crispr_obj, dims = 1:15)
crispr_obj <- RunUMAP(crispr_obj, dims = 1:15)


DimPlot(crispr_obj, reduction = "umap") + theme(legend.position = "bottom")
ggsave("umap.pdf", plot = last_plot(), width = 5, height = 5, units = "in", dpi = 300)

crispr_obj$sgRNA_identity <- factor(crispr_obj$sgRNA_identity, 
                                    levels = c("NT1", "NT2", "AAVS", 
                                               "ITGB1-sg1", "ITGB1-sg2", "ITGB1-sg3",
                                               "RPS6-sg1", "RPS6-sg2", "RPS6-sg3",
                                               "GATA1-sg1", "GATA1-sg2", "GATA1-sg3"))
DimPlot(crispr_obj, reduction = "umap", group.by = "sgRNA_identity") + theme(legend.position = "bottom")
ggsave("umap_groupby_sgRNA.pdf", plot = last_plot(), width = 6, height = 7, units = "in", dpi = 300)
DimPlot(crispr_obj, reduction = "umap", split.by = "sgRNA_identity", ncol = 3) +
  theme(legend.position = "bottom")
ggsave("umap_splitby_sgRNA.pdf", plot = last_plot(), width = 10, height = 12, units = "in", dpi = 300)

saveRDS(crispr_obj, file = "crispr_obj_seurat.rds")



##### sgRNA
setwd("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/seurat_result/sgRNA/")

crispr_obj <- readRDS(file = "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/seurat_result/seurat/crispr_obj_seurat.rds")

table(crispr_obj@meta.data[["sgRNA_identity"]])


# Set the assay to sgRNA and perform normalization and scaling
DefaultAssay(crispr_obj) <- "sgRNA"
crispr_obj <- NormalizeData(crispr_obj, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(crispr_obj)
crispr_obj <- ScaleData(crispr_obj, features = all.genes)

crispr_obj <- RunPCA(crispr_obj, features = all.genes, reduction.name = "sgRNA_PCA", reduction.key = "sgRNAPCA_")
print(crispr_obj[["sgRNA_PCA"]], dims = 1:5, nfeatures = 10)  
ElbowPlot(crispr_obj, reduction = "sgRNA_PCA") 

crispr_obj <- FindNeighbors(crispr_obj, dims = 1:11, reduction = "sgRNA_PCA")  
crispr_obj <- FindClusters(crispr_obj, resolution = 0.5)  
table(crispr_obj@meta.data[["seurat_clusters"]])  

crispr_obj <- RunUMAP(crispr_obj, dims = 1:11, min.dist = 0.5, reduction = "sgRNA_PCA", reduction.name = "sgRNA_UMAP", reduction.key = "sgRNAUMAP_")

saveRDS(crispr_obj, file = "crispr_obj_seurat_sgRNA.rds")



setwd("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/figure/")

crispr_obj <- readRDS(file = "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/seurat_result/sgRNA/crispr_obj_seurat_sgRNA.rds")
crispr_obj_seurat <- readRDS(file = "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/seurat_result/seurat/crispr_obj_seurat.rds")
K562<- Read10X("/public/home/wangxinlong/Project/crispr1/data/20250102ABFC20221434-76_K562/sc20241231K562_addsequence/Ab-1_L2_S2401B2401_DBEC_MolsPerCell_MEX")
K562_seurat<-CreateSeuratObject(counts =K562$`Gene Expression`)
pal <-c("#E3A8A2",'#75A3D1','#E9AD95', "#65B8AD")
pal <-c('#75A3D1','#E3A8A2')
pal <-c( '#65B8AD','#DC7C56')



##### SupplementaryFigure7 b
#nCount_sgRNA
# # Set the assay to sgRNA and perform normalization and scaling
DefaultAssay(crispr_obj) <- "sgRNA"
crispr_obj <- NormalizeData(crispr_obj, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(crispr_obj)
crispr_obj <- ScaleData(crispr_obj, features = all.genes)
crispr_obj <- RunPCA(crispr_obj, features = all.genes, reduction.name = "sgRNA_PCA", reduction.key = "sgRNAPCA_")
print(crispr_obj[["sgRNA_PCA"]], dims = 1:5, nfeatures = 10)  
ElbowPlot(crispr_obj, reduction = "sgRNA_PCA") 

crispr_obj <- FindNeighbors(crispr_obj, dims = 1:11, reduction = "sgRNA_PCA")  
crispr_obj <- FindClusters(crispr_obj, resolution = 0.5)  
table(crispr_obj@meta.data[["seurat_clusters"]])  

crispr_obj <- RunUMAP(crispr_obj, dims = 1:11, min.dist = 0.5,reduction = "sgRNA_PCA", reduction.name = "sgRNA_UMAP", reduction.key = "sgRNAUMAP_")

mypal <- c("#AED6F1", '#6BB3E6', "#75A3D1", '#1E75B3', '#E3A8A2', '#CF5B5B', '#CC7177', '#B23E4D', '#51132F', '#BDE0B6', '#77BB88', '#55AA6A', "#65B8AD", '#4AA195', '#326C64', '#E69D70', '#DA712F', "#B65A20", '#F78426', "#E26A08", '#B15306', '#7F3C1A')

sgRNA_colors <- c(
  "ITGB1-sg1" = "#AED6F1",  
  "ITGB1-sg2" = '#6BB3E6',  
  "ITGB1-sg3" = "#4F8AC4",  
  "RPS6-sg1"  = "#EB9593",
  "RPS6-sg2"  = "#C26171",
  "RPS6-sg3"  = '#9A3C4C',  
  "GATA1-sg1" =  '#65B8AD',  
  "GATA1-sg2" = "#4AA195",  
  "GATA1-sg3" = "#326C64",  
  "NT1"      = "#FCB99C" ,
  "NT2" = "#FCB99C" ,
  "AAVS" ="#FCB99C" 
)

DimPlot(
  crispr_obj, 
  reduction = "sgRNA_UMAP", 
  group.by = "sgRNA_identity", 
  cols = sgRNA_colors,
  label = TRUE, 
  repel = TRUE
) + 
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 0.5)
  ) 
ggsave("sgRNA_UMAP.pdf", plot = last_plot(), width = 6.5, height = 5, units = "in", dpi = 300)





##### SupplementaryFigure7
library(dplyr)
library(ggplot2)

setwd("/public/home/wangxinlong/Project/crispr1/result/compare_seq_K562/process_plot_data/figure_3/")

file_path <- "/public/home/wangxinlong/Project/crispr1/result/compare_seq_K562/process_plot_data/merge/combined_merged_all.tsv"

combined_merged <- read.delim(file_path, stringsAsFactors = FALSE)

combined_merged <- combined_merged[combined_merged$group != "GSM4367984", ]
combined_merged$group <- factor(combined_merged$group, levels = c("PerturbPro", "Chromium", "Perturb", "ECCITE"))
combined_merged <- combined_merged %>% arrange(group)

custom_colors <- c(
  "Chromium" = "darkorange",
  "PerturbPro" = "steelblue3",
  "Perturb" = "mediumpurple3",
  "ECCITE" = "seagreen3"
)



##### SupplementaryFigure7 d
# x = nRead_RNA, y = nCount_RNA
ggplot(combined_merged, aes(x = nRead_RNA, y = nCount_RNA, color = group)) +
  geom_point(alpha = 0.5, size = 1) +
  labs(
    x = "Read count",
    y = "Number of UMIs",
    color = "Group"
  ) +
  scale_color_manual(values = custom_colors) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
ggsave("nReadRNA_nCountRNA.pdf", plot = last_plot(), width = 8, height = 6, units = "in", dpi = 300)



##### SupplementaryFigure7 e
# boxplot UMIs
ggplot(combined_merged, aes(x = group, y = nCount_RNA, fill = group)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 1) +
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = round(..y.., 0)),  
    vjust = -0.5,                  
    size = 4,
    color = "black"
  ) +
  labs(
    x = NULL,
    y = "Number of UMIs"
  ) +
  scale_fill_manual(values = custom_colors) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  )
ggsave("boxplot_UMIs.pdf", plot = last_plot(), width = 6, height = 6, units = "in", dpi = 300)



# boxplot percent.mt
ggplot(combined_merged, aes(x = group, y = percent.mt, fill = group)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 1) +
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = round(..y.., 0)),  
    vjust = -0.5,                  
    size = 4,
    color = "black"
  ) +
  labs(
    x = NULL,
    y = "Percent_MT"
  ) +
  scale_fill_manual(values = custom_colors) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  )
ggsave("boxplot_Percent_MT.pdf", plot = last_plot(), width = 6, height = 6, units = "in", dpi = 300)



# boxplot percent_intron
ggplot(combined_merged, aes(x = group, y = percent_intron, fill = group)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 1) +
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = round(..y.., 0)),  
    vjust = -0.5,                  
    size = 4,
    color = "black"
  ) +
  labs(
    x = NULL,
    y = "Percent_Intron"
  ) +
  scale_fill_manual(values = custom_colors) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  )
ggsave("boxplot_Percent_Intron.pdf", plot = last_plot(), width = 6, height = 6, units = "in", dpi = 300)
