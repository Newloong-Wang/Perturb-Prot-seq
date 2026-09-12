##### SupplementaryFigure2 e 
library(tidyverse)

df <- read.csv("/public/home/wangxinlong/Project/crispr1/data/293TbulkRNA/筛选差异结果.csv")

df$group <- df$group %>%
  str_replace_all("groupA", "fresh") %>%
  str_replace_all("groupB", "k1") %>%
  str_replace_all("groupC", "k4")

df <- df %>%
  separate(group, into = c("group1", "group2"), sep = "_")

groups <- c("fresh", "k1", "k4")

mat <- matrix(0, nrow = 3, ncol = 3, dimnames = list(groups, groups))
for (i in 1:nrow(df)) {
  mat[df$group1[i], df$group2[i]] <- df$num[i]
}

mat_long <- as.data.frame(as.table(mat))
colnames(mat_long) <- c("group1", "group2", "DEG")
mat_long$group1 <- factor(mat_long$group1, levels = rev(groups))  
mat_long$group2 <- factor(mat_long$group2, levels = groups)       

mat_long <- mat_long %>%
  mutate(
    upper = as.integer(group2) > (4 - as.integer(group1)),
    diagonal = group1 == group2
  )

ggplot(mat_long, aes(x = group2, y = group1)) +
  geom_tile(aes(fill = upper | diagonal), color = "black") +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "grey80"), guide = "none") +
  geom_text(data = filter(mat_long, !upper), aes(label = DEG), color = "black", size = 6) +
  geom_text(data = filter(mat_long, diagonal), aes(label = DEG), color = "black", size = 6) +
  coord_fixed() +
  scale_x_discrete(position = "top") +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0, size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 18, color = "black", face = "plain"),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  labs(title = "Differentially expressed genes", x = NULL, y = NULL)
ggsave("/public/home/wangxinlong/Project/crispr1/data/293TbulkRNA/DEG_number.pdf", plot = last_plot(), width = 6, height = 6, units = "in", dpi = 300)





##### SupplementaryFigure2 f
library(ggplot2)
library(dplyr)
library(tidyr)

average <- read.csv("/public/home/wangxinlong/Project/crispr1/result/20241118bamanotation_bulkRNA/statistic/average_statistic_summary.csv",  
                    check.names = FALSE, stringsAsFactors = FALSE)

average_long <- average %>%
  tidyr::pivot_longer(
    cols = -Group,
    names_to = "Region",
    values_to = "Reads"
  ) %>%
  group_by(Group) %>%
  mutate(Percent = Reads / sum(Reads) * 100)

average_long$Region <- factor(average_long$Region, 
                              levels = c("Exonic", "Intronic", "Intergenic", "Ambiguous", "Unmapped"))

custom_colors <- c(
  "#b2b0af",  # Exonic
  "#dda851",  # Intronic
  "#accbe0",  # Intergenic
  "#7dcbae",  # Ambiguous
  "#ebb9a5"   # Unmapped
)

ggplot(average_long, aes(x = Group, y = Percent, fill = Region)) +
  geom_bar(stat = "identity", width = 0.7) +  
  geom_text(aes(label = sprintf("%.1f", Percent)),  
            position = position_stack(vjust = 0.5), 
            size = 4, color = "black") +
  scale_fill_manual(values = custom_colors) +  
  labs(y = "Percent of reads", x = "Sample", fill = "Region") +
  scale_y_continuous(breaks = seq(0, 100, 25), labels = c("0", "25", "50", "75", "100")) + 
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(hjust = 0.5, vjust = 0.5),  
    legend.position = "right",  
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12)
  )

ggsave("/public/home/wangxinlong/Project/crispr1/result/20241118bamanotation_bulkRNA/picture/generegion.pdf", plot = last_plot(), width = 8, height = 8, units = "in", dpi = 300)
