##### SupplementaryFigure11 b
### go_plot
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(reshape2)

dir.create("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_3group/GO/")
setwd("/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_3group/GO/")



### positive

# All
file_path <- "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale/GATA1_NT_go_kegg/GATA1_NT_compareCluster_BP.tsv"
go_data <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
target_descriptions <- c(
  "positive regulation of mononuclear cell proliferation",
  "regulation of mononuclear cell proliferation",
  "positive regulation of leukocyte proliferation"
)
item <- go_data %>% filter(Cluster == "positive" & Description %in% target_descriptions)
item$group <- "All"
go_positive <- item

# remove_top20
file_path <- "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_top20middlebottoom20/go_kegg_mixscale_GATA1_top20middlebottoom20/mixscale_GATA1_compareCluster_BP.tsv"
go_data <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
target_descriptions <- c(
  "positive regulation of mononuclear cell proliferation",
  "positive regulation of leukocyte proliferation",
  "regulation of mononuclear cell proliferation"
)
item <- go_data %>% filter(Cluster == "positive" & Description %in% target_descriptions)
item$group <- "remove_top20"
go_positive <- rbind(go_positive, item)

# remove_top30
file_path <- "/public/home/wangxinlong/Project/crispr1/result/20250102ABFC20221434-76_K562/mixscale_GATA1_top30middlebottoom30/go_kegg_mixscale_GATA1_top30middlebottoom30/mixscale_GATA1_compareCluster_BP.tsv"
go_data <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
target_descriptions <- c(
  "positive regulation of mononuclear cell proliferation",
  "positive regulation of leukocyte proliferation",
  "regulation of mononuclear cell proliferation"
)
item <- go_data %>% filter(Cluster == "positive" & Description %in% target_descriptions)
item$group <- "remove_top30"
go_positive <- rbind(go_positive, item)

write.table(go_positive,
            file = "go_positive_3group.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


# UpregulatedGOtermsinGATA1sg_3group
go_positive$group <- factor(go_positive$group, levels = c("All", "remove_top20", "remove_top30"))
go_positive <- go_positive[order(go_positive$Description, go_positive$group), ]

go_positive$y_label <- factor(seq_len(nrow(go_positive)), levels = rev(seq_len(nrow(go_positive))))

group_colors <- c(
  "All" = "#FAD4D4",         
  "remove_top20" = "#F08080", 
  "remove_top30" = "#DC143C"
)

go_positive$label_combined <- paste0(
  go_positive$Description,
  str_dup(" ", 5),
  "p.adjust=", signif(go_positive$p.adjust, 2)
)

ggplot(go_positive, aes(x = Count, y = y_label, fill = group)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = label_combined),
            x = 0.04, hjust = 0, vjust = 0.5, size = 4, color = "black") +
  geom_text(aes(label = geneID),
            x = 0.04, hjust = 0, vjust = 4, size = 3.5,
            color = "#EF3B2C") +
  scale_fill_manual(values = group_colors) +
  scale_y_discrete(labels = NULL) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.x = element_text(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom"
  ) +
  labs(title = "Upregulated GO terms in GATA1sg", x = "Count", y = "") +
  coord_cartesian(clip = "off", expand = TRUE)
ggsave(filename = "UpregulatedGOtermsinGATA1sg_3group.pdf", plot = last_plot(), width = 10, height = 8.3, dpi = 300)
