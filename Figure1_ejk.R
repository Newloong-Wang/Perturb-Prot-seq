##### Figure1 e
setwd("/public/home/wangxinlong/Project/crispr1/data/293TbulkRNA/featurecounts/")

library(stringr)
library(tidyverse)
library(paletteer)
library(corrplot)
library(tidyverse)
library(readxl)
library(dplyr)
library(pheatmap)

filename<-list.files()
var_name<-c()
raw.all<-data.frame()
for (i in 1:length(filename)){
  var_name[i]<-gsub('?_.*counts.txt','',filename[i])
  var_name[i]<-str_replace_all(var_name[i], "-", ".")
  assign(var_name[i],read.table(filename[i],sep="\t",header=TRUE,col.names=c("Geneid","Chr","Start","End","Strand","Length",var_name[i])) %>% dplyr::select(-Chr,-Start,-End,-Strand)) 
  print(i)
  a<-get(var_name[1])
  raw.all<-data.frame(a[,c("Geneid","Length")])
}
for(i in var_name){
  a<-get(i)
  raw.all<-merge(raw.all,a,by=c("Geneid","Length"))
  rm(i)
}


cor.sp<-cor(raw.all[3:length(colnames(raw.all))],method="pearson")
corrplot(cor.sp,method="color",type="full",addCoef.col = "white",sig.level = 0.01)

# TPM normalization
tpmMatrix <- raw.all %>%
  pivot_longer(c(-Geneid, -Length),
               names_to = "Group",
               values_to = "SampleCounts") %>%
  group_by(Group) %>%
  mutate(SampleTPM = (((SampleCounts/Length)*1e6)/sum(SampleCounts/Length))) %>%
  pivot_wider(id_cols = "Geneid",
              names_from = "Group",
              values_from = "SampleTPM")

print(colSums(tpmMatrix[2:10]))

zero_rows <- tpmMatrix %>%
  dplyr::select(-Geneid) %>%        
  apply(1, function(x) all(x == 0))  
sum(zero_rows)  

tpmMatrix_filtered <- tpmMatrix[!zero_rows, ]

tpm_numeric <- tpmMatrix_filtered[, -1]

cor_matrix <- cor(tpm_numeric, method = "pearson")

sample_names <- colnames(cor_matrix)

sample_group <- ifelse(grepl("fresh", sample_names), "fresh",
                       ifelse(grepl("k1", sample_names), "k1", "k4"))

annotation_col <- data.frame(Group = factor(sample_group))
rownames(annotation_col) <- sample_names

ann_colors <- list(
  Group = c(fresh = "#FFCC00", k1 = "#1F78B4", k4 = "#33A02C")
)



p <- pheatmap(cor_matrix,
              cluster_rows = TRUE,
              cluster_cols = FALSE,
              clustering_distance_rows = "correlation",
              clustering_method = "complete",
              silent = TRUE)

row_order <- p$tree_row$order

ordered_matrix <- cor_matrix[row_order, row_order]

annotation_ordered <- annotation_col[row_order, , drop = FALSE]

pdf("/public/home/wangxinlong/Project/crispr1/data/293TbulkRNA/pearson.pdf", width = 7, height = 6)
pheatmap(ordered_matrix,
         annotation_col = annotation_ordered,
         annotation_row = annotation_ordered,
         annotation_colors = ann_colors,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         clustering_distance_rows = "correlation",
         clustering_method = "complete",
         display_numbers = FALSE,
         color = colorRampPalette(c("#EFF3FF", "#3182bd"))(100),
         border_color = NA,
         angle_col = 90,
         cellwidth = 30,
         cellheight = 30)
dev.off()





##### Figure1 j
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(paletteer)
library(writexl)
setwd('/public/home/wangxinlong/Project/crispr1/result/lung/')

###seurat
lung<-Read10X("/public/home/wangxinlong/Project/crispr1/data/20241217_sclung/Ab-5_L4_Q0028W0072_DBEC_MolsPerCell_MEX/")
lung_obj<-CreateSeuratObject(counts =lung$`Gene Expression`)
lung_obj[['Protein']] = CreateAssayObject(counts = lung$`Antibody Capture`, project="Lung")

Assays(lung_obj)
rownames(lung_obj[["Protein"]])
lung_obj[["percent.mt"]] <- PercentageFeatureSet(lung_obj, pattern = "^mt")

plot1 <- FeatureScatter(lung_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(lung_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

lung_obj_fliter <- subset(lung_obj, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 20)

VlnPlot(lung_obj_fliter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

lung_obj_fliter<-NormalizeData(lung_obj_fliter,normalization.method = "LogNormalize",scale.factor = 10000)
lung_obj_fliter<-FindVariableFeatures(lung_obj_fliter,selection.method = "vst",nfeatures=700)

marker.list<-c("Vmf","Il1b","Aplnr","Gpihbp1","Cd209a","Flt3","Tbx2","Ednrb","Apln","Lamp3","Scgb1a1","Scgb3a2","Bpifb1","Sftpd","Abca3","Foxj1","Ccdc113","Sftpc","Car4","Vwf","Slc6a2","Thy1","Pdpn","Prox1","Flt4","Mmrn1","Fxyd6","Rtkn2","Ager","Hopx","Pdpn","Retnlg","Ear2","Marco","Abcg1","Dcn","Col14a1","Msln","Gucy1a1","Notch3","Enpp2","Pdgfrb","Pi16","Cygb","Acta2","Tagln","Npnt","Col13a1","Itga8","Cxcl14","Ms4a2","Plac8","Ly6c2","Ptprc","Epcam","Scgb1a1","Mki67","C1qa","C1qb","C1qc","Cd3g","Cd3d","Cd3e","Ms4a1","Bank1","Nkg7","Klrd1",'Hopx','Cd4','Cd19','Cd79a','Itgam','Cd7','Flt3','Hhip')
current_var_features <- VariableFeatures(lung_obj_fliter)
updated_var_features <- unique(c(current_var_features, marker.list))
VariableFeatures(lung_obj_fliter) <- updated_var_features
VariableFeatures(lung_obj_fliter)

top10<-head(VariableFeatures(lung_obj_fliter),10)
plot1<-VariableFeaturePlot(lung_obj_fliter)
plot2<-LabelPoints(plot=plot1,points = top10)
plot1+plot2

all.genes<-rownames(lung_obj_fliter)
lung_obj_fliter<-ScaleData(lung_obj_fliter,features=all.genes)
lung_obj_fliter<-RunPCA(lung_obj_fliter,features=VariableFeatures(object = lung_obj_fliter),npcs =100)
VizDimLoadings(lung_obj_fliter,dims=1:2,reduction = "pca")
ElbowPlot(lung_obj_fliter,ndims=80)

lung_obj_fliter<-FindNeighbors(lung_obj_fliter,dims=1:26)
lung_obj_fliter<-FindClusters(lung_obj_fliter,resolution = 1.5)
lung_obj_fliter<-RunUMAP(lung_obj_fliter,dims = 1:26)
DimPlot(lung_obj_fliter,reduction = "umap",label =TRUE)
table(lung_obj_fliter@meta.data[["RNA_snn_res.1.5"]])

FeaturePlot(lung_obj_fliter, features = c("nFeature_RNA"))



#Cell type annotation
new.cluster.ids <- c(
  "AT2",          # 0
  "gCap",         # 1
  "Col13a1",      # 2
  "Col13a1",      # 3
  "Pericyte",     # 4
  "aCap",         # 5
  "Club",         # 6
  "AT2",          # 7
  "T",            # 8
  "Col14a1",      # 9
  "Col13a1",      # 10
  "Monocyte",     # 11
  "V_SMC",        # 12
  "Ciliated",     # 13
  "AT2",          # 14
  "EC_Vein",      # 15
  "Pericyte",     # 16
  "aM",           # 17
  "NK",           # 18
  "EC_Artery",    # 19
  "Col13a1",      # 20
  "AT1",          # 21
  "Mesothelial",  # 22
  "B",            # 23
  "Innate_Lymphoid", # 24
  "Pericyte",     # 25
  "A_SMC",        # 26
  "Neutrophil",  # 27
  "EC_Lymph"      # 28
)
names(new.cluster.ids) <- levels(lung_obj_fliter)
lung_obj_fliter <- RenameIdents(lung_obj_fliter, new.cluster.ids)
lung_obj_fliter$celltype <- Idents(lung_obj_fliter)
DimPlot(lung_obj_fliter,reduction = "umap",label = TRUE, group.by = "celltype")

ordered.season<-c("AT1","AT2",'Ciliated',"Club","aCap","gCap","EC_Artery","EC_Vein",'EC_Lymph',"Pericyte",'Col13a1','Col14a1','V_SMC','A_SMC',"Mesothelial","T",'B','Neutrophil','aM','NK','Innate_Lymphoid','Monocyte')
Idents(lung_obj_fliter) <- factor(Idents(lung_obj_fliter), levels= ordered.season)

mypal <- c("#AED6F1", "#6BB3E6", "#75A3D1", "#1E75B3", 
           "#B23E4D", "#E3A8A2", "#CC7177", "#CF5B5B", "#93344F", 
           "#BDE0B6", "#77BB88", "#448855", "#65B8AD", "#3A7E75", "#2A5B54", 
           "#f8b595", "#B65A20", "#643112", "#B5420D","#F86930", "#CE4A0D", 
           '#F59061')
DimPlot(
  lung_obj_fliter,
  reduction = "umap",
  label = TRUE, 
  pt.size = 1, 
  cols = mypal, 
  label.size = 4, 
  label.box = FALSE, 
  repel = TRUE
) + 
  labs(x = "UMAP1", y = "UMAP2") +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid")
  ) +
  guides(col = guide_legend(ncol = 1))
ggsave("UMAP1.pdf", plot = last_plot(), width = 11, height = 9, units = "in", dpi = 300)



# celltype_cellbarcode
cellbarcode <- rownames(lung_obj_fliter@meta.data)
celltype <- lung_obj_fliter@meta.data$celltype
cell_info <- data.frame(cellbarcode = cellbarcode,
                        celltype = celltype,
                        stringsAsFactors = FALSE)
cell_summary <- aggregate(cellbarcode ~ celltype, data = cell_info,
                          FUN = function(x) paste(x, collapse = ","))
cell_summary$celltype <- factor(cell_summary$celltype, levels = ordered.season)
cell_summary <- cell_summary[order(cell_summary$celltype), ]
write.table(cell_summary, file = "celltype_cellbarcode.tsv", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
write_xlsx(cell_summary, path = "celltype_cellbarcode.xlsx")

DefaultAssay(lung_obj_fliter) <- "RNA"
marker.genes<-FindAllMarkers(lung_obj_fliter, group.by = "celltype", only.pos = T,min.pct = 0.3,logfc.threshold = 1)
marker.genes <- marker.genes %>% filter(p_val_adj < 0.05)



# marker gene ----
ordered.season <- c("AT1", "AT2", "Ciliated", "Club", 
                    "aCap", "gCap", "EC_Artery", "EC_Vein", "EC_Lymph", 
                    "Pericyte", "Col13a1", "Col14a1", "V_SMC", "A_SMC", 
                    "Mesothelial", "T", "B", "Neutrophil", "aM", "NK", 
                    "Innate_Lymphoid", "Monocyte")

Idents(lung_obj_fliter) <- factor(Idents(lung_obj_fliter), levels = ordered.season)
lung_obj_fliter@commands$NormalizeData.RNA 

DotPlot(lung_obj_fliter, 
        features = c(
          "Ager", "Hopx",      # AT1
          "Lamp3", "Sftpc",    # AT2
          "Foxj1", "Rsph1",    # Ciliated
          "Scgb1a1", "Scgb3a2",# Club
          "Ednrb", "Car4",     # aCap
          "Aplnr", "Gpihbp1",  # gCap
          "Gja5", "Olfm2",     # EC_Artery
          "Slc6a2", "Ephb4",   # EC_Vein
          "Prox1", "Flt4",     # EC_Lymph
          "Gucy1a1", "Notch3", # Pericyte
          "Cxcl14", "Col13a1", # Col13a1
          "Col14a1", "Pi16",   # Col14a1
          "Cnn1", "Ntrk3",     # V_SMC
          "Igf1", "Hhip",      # A_SMC
          "Msln", "Wt1",       # Mesothelial
          "Cd3e", "Cd4",       # T 
          "Cd79a", "Cd19",     # B 
          "Il1b", "Retnlg",    # Neutrophil
          "Ear2", "Abcg1",     # aM
          "Nkg7", "Klrd1",     # NK
          "Gata3", "Cd69",     # Innate Lymphoid Cells
          "Itgam", "Csf1r"     # Monocyte
        )) + 
  scale_color_gradient(low = "white", high = "#B23E4D") +
  coord_flip() + 
  theme_minimal() + 
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(), 
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5)
  )
ggsave("dotplot_2markerGenes.pdf", plot = last_plot(), width = 10, height = 12, units = "in", dpi = 300)



##### Figure1 k
# Protein ----
library(scales)
library(colorspace)
library(ggsignif)
                            
DefaultAssay(lung_obj_fliter) <- "Protein"
lung_obj_fliter <- NormalizeData(lung_obj_fliter, normalization.method = "CLR", margin = 2)                          
fill_colors <- c("#E3A8A2", "#75A3D1", "#E69D70", "#65B8AD")
border_colors <- darken(fill_colors, amount = 0.3)

lung_obj_fliter$cell_type <- ifelse(lung_obj_fliter$celltype %in% c("AT1", "AT2", "Ciliated", "Club"), "Epithelial", 
                                    ifelse(lung_obj_fliter$celltype %in% c("aCap", "gCap", "EC_Artery", "EC_Vein", "EC_Lymph"), "Endothelial", 
                                           ifelse(lung_obj_fliter$celltype %in% c("Col13a1", "Col14a1", "V_SMC", "A_SMC", "Mesothelial",'Pericyte'), "Mesenchymal", 
                                                  "Immune")))
lung_obj_fliter$celltype <- factor(lung_obj_fliter$celltype, levels = ordered.season)


# Epcam Protein
expression_data <- FetchData(lung_obj_fliter, vars = c("Epcam-940269-pAbO", "cell_type"))
comparisons <- list(
  c("Epithelial", "Endothelial"),
  c("Epithelial", "Immune"),
  c("Epithelial", "Mesenchymal")
)
ggplot(expression_data, aes(x = cell_type, y = `Epcam-940269-pAbO`, fill = cell_type)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_signif(
    comparisons = comparisons,
    test = "wilcox.test",
    test.args = list(exact = FALSE),
    map_signif_level = function(p) sprintf("p = %.2e", p),  
    step_increase = 0.1
  ) +
  scale_fill_manual(values = fill_colors) +
  theme_classic() +
  labs(y = "Expression Level", x = "Cell Type") +
  theme(
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  NoLegend()
ggsave("vln_Epcam_Pro_1_expression.pdf", plot = last_plot(), width = 6, height = 5, units = "in", dpi = 300)


# CD31 Protein
expression_data <- FetchData(lung_obj_fliter, vars = c("CD31-940068-pAbO", "cell_type"))
comparisons <- list(
  c("Endothelial", "Epithelial"),
  c("Endothelial", "Immune"),
  c("Endothelial", "Mesenchymal")
)
ggplot(expression_data, aes(x = cell_type, y = `CD31-940068-pAbO`, fill = cell_type)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_signif(
    comparisons = comparisons,
    test = "wilcox.test",
    test.args = list(exact = FALSE),
    map_signif_level = function(p) sprintf("p = %.2e", p),  
    step_increase = 0.1
  ) +
  scale_fill_manual(values = fill_colors) +
  theme_classic() +
  labs(y = "Expression Level", x = "Cell Type") +
  theme(
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)    
  ) +
  NoLegend()
ggsave("vln_CD31_Pro_1-expression.pdf", plot = last_plot(), width = 6, height = 5, units = "in", dpi = 300)


# Hopx Protein
Epi_clusters <- subset(lung_obj_fliter, subset = cell_type == "Epithelial")
expression_data <- FetchData(Epi_clusters, vars = c("Hopx-940299-pAbO", "celltype"))
comparisons <- list(
  c("AT1", "AT2"),
  c("AT1", "Ciliated"),
  c("AT1", "Club")
)
ggplot(expression_data, aes(x = celltype, y = `Hopx-940299-pAbO`, fill = celltype)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_signif(
    comparisons = comparisons,
    test = "wilcox.test",
    test.args = list(exact = FALSE),
    map_signif_level = function(p) sprintf("p = %.2e", p),  
    step_increase = 0.1
  ) +
  scale_fill_manual(values = mypal) +
  theme_classic() +
  labs(y = "Expression Level", x = "Cell Type") +
  theme(
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)    
  ) +
  NoLegend()
ggsave("vln_Hopx_Pro_expression.pdf", plot = last_plot(), width = 6, height = 5, units = "in", dpi = 300)


# Scgb1a1 Protein
Epi_clusters <- subset(lung_obj_fliter, subset = cell_type == "Epithelial")
expression_data <- FetchData(Epi_clusters, vars = c("Scgb1a1-940272-pAbO", "celltype"))
comparisons <- list(
  c("Club", "AT1"),
  c("Club", "AT2"),
  c("Club", "Ciliated")
)
ggplot(expression_data, aes(x = celltype, y = `Scgb1a1-940272-pAbO`, fill = celltype)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_signif(
    comparisons = comparisons,
    test = "wilcox.test",
    test.args = list(exact = FALSE),
    map_signif_level = function(p) sprintf("p = %.2e", p),  
    step_increase = 0.1
  ) +
  scale_fill_manual(values = mypal) +
  theme_classic() +
  labs(y = "Expression Level", x = "Cell Type") +
  theme(
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)    
  ) +
  NoLegend()
ggsave("vln_Scgb1a1_Pro_expression.pdf", plot = last_plot(), width = 6, height = 5, units = "in", dpi = 300)
