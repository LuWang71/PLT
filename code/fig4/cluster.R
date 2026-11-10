library(tidyr)
library(pheatmap)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(ggplot2)
group$Timepoint0 <- ifelse(group$Group == "W/ Infection", 
                           paste0("W_", group$Timepoint), 
                           ifelse(group$Group == "W/O Infection", 
                                  paste0("WO_", group$Timepoint), 
                                  group$Timepoint))

bacteria = OTUs[OTUs$Kingdom == 'Bacteria',]
fungi = OTUs[OTUs$Kingdom == 'Eukaryota',]
#archaea = OTUs[OTUs$Kingdom == 'Archaea',]
bac  = bacteria
bac = bac[,-c(1,3)]
bac = aggregate(. ~ Genus, data = bac, FUN = sum)
rownames(bac) = bac$Genus
bac = bac[,-1]

fun  = fungi
fun = fun[,-c(1,3)]
fun = aggregate(. ~ Genus, data = fun, FUN = sum)
rownames(fun) = fun$Genus
fun = fun[,-1]

bac_names = c("Anaerobutyricum", "Schaalia", "Segatella", "Mediterraneibacter",
        "Phocaeicola", "Intestinibacter", "Alistipes", "Escherichia",
        "Enterobacter", "Morganella", "Dysgonomonas", "Bifidobacterium",
        "Raoultella", "Bacteroides", "Lacticaseibacillus", "Lactobacillus",
        "Sarcina", "Enterocloster", "Weissella", "Akkermansia", "Hungatella",
        "Pseudomonas", "Enterococcus", "Roseburia", "Faecalibacterium",
        "Ligilactobacillus", "Thomasclavelia", "Citrobacter", "Megamonas",
        "Rothia", "Blautia", "Abiotrophia", "Parabacteroides", "Staphylococcus",
        "Mycobacterium", "Streptococcus", "Kluyvera", "Klebsiella", "Yersinia",
        "Haemophilus", "Veillonella", "Stenotrophomonas")

fun_names = c("Fusarium", "Malassezia", "Geosmithia", "Uncinocarpus",
              "Meyerozyma", "Kwoniela", "Saccharomyces", "Puccinia", "Botrytis",
              "Pyricularia", "Ramularia", "Trichosporon", "Mycotypha", "Letharia",
              "Colletotrichum", "Neurospora", "Candida", "Rhizophagus",
              "Phycomyces", "Suillus", "Talaromyces", "Tetrapisispora",
              "Alternaria", "Xylona", "Pichia", "Penicillium",
              "unclassified_Saccharomycetales", "Lodderomyces", "Aspergillus",
              "Nakaseomyces", "Drepanopeziza", "Dissoconium", "Diutina",
              "Kluyveromyces", "Clavispora", "Meira")
names = c(bac_names,fun_names)
fun0 = fun[rownames(fun) %in% fun_names,]
bac0 = bac[rownames(bac) %in% bac_names,]

fun0 = t(fun0)
bac0 = t(bac0)

cor_matrix <- cor(fun0, bac0, method = "spearman")
brks <- c(seq(-0.6, -0.01, by = 0.001), seq(0, 0.6, by = 0.001))

pheatmap(cor_matrix,
         color = colorRampPalette(c("#2171B5", "white", "#D94801"))(100),
         main = "Fungi-Bacteria Interaction Heatmap", # 热图标题
         # color = c(colorRampPalette(colors = c("#2171B5", "white"))(length(brks)/2),
         #           colorRampPalette(colors = c("white", "#D94801"))(length(brks)/2)),
         # legend_breaks = seq(-0.6, 0.6, 0.1),
         # breaks = brks,
         cluster_rows = TRUE, 
         cluster_cols = TRUE)
#-----------------------------------------------------------------------

pheatmap_result = pheatmap(cor_matrix,
                           color = colorRampPalette(c("#2171B5", "white", "#D94801"))(100), 
                           main = "Fungi-Bacteria Interaction Heatmap", 
                           # color = c(colorRampPalette(colors = c("#2171B5", "white"))(length(brks)/2),
                           #           colorRampPalette(colors = c("white", "#D94801"))(length(brks)/2)),
                           # legend_breaks = seq(-0.6, 0.6, 0.1),
                           # breaks = brks,
                           cluster_rows = TRUE, 
                           cluster_cols = TRUE)

row_clusters <- cutree(pheatmap_result$tree_row, k = 4) 
col_clusters <- cutree(pheatmap_result$tree_col, k = 4)


print(row_clusters)
print(col_clusters)



fungi_1 = c('Ramularia','Trichosporon','Clavispora','Diutina',
            'Malassezia','Geosmithia',
            'Pichia','Neurospora','Pyricularia')
fungi_2 = c('Botrytis','Rhizophagus','Aspergillus','Fusarium',
            'Suillus','Colletotrichum','Penicillium','Letharia',
            'Talaromyces','Puccinia','Tetrapisispora')
fungi_3 = c('Dissoconium','Saccharomyces','Xylona','Meyerozyma',
            'Nakaseomyces','Kluyveromyces','Lodderomyces')
fungi_4 = c('Candida','Uncinocarpus','Drepanopeziza','Alternaria',
            'Mycotypha','Phycomyces','Meira','unclassified_Saccharomycetales')

bac_1 = c('Sarcina','Akkermansia','Segatella','Alistipes','Enterocloster',
          'Hungatella','Faecalibacterium','Roseburia','Blautia','Mediterraneibacter',
          'Anaerobutyricum','Intestinibacter','Dysgonomonas','Weissella',
          'Ligilactobacillus')
bac_2 = c('Thomasclavelia','Bacteroides','Phocaeicola','Parabacteroides')
bac_3 = c('Pseudomonas','Kluyvera','Stenotrophomonas','Schaalia',
          'Bifidobacterium','Veillonella','Morganella','Megamonas',
          'Staphylococcus','Lactobacillus','Abiotrophia','Haemophilus',
          'Lacticaseibacillus')
bac_4 = c('Mycobacterium','Klebsiella','Yersinia','Escherichia',
          'Raoultella','Citrobacter','Enterobacter','Streptococcus',
          'Enterococcus','Rothia')
fun_names = c(fungi_1,fungi_2,fungi_3,fungi_4)
bac_names = c(bac_1,bac_2,bac_3,bac_4)

all = rbind(fun0,bac0)

all = na.omit(all)
names = c(fun_names,bac_names)
mapping <- list(
  fungi_1 = fungi_1,
  fungi_2 = fungi_2,
  fungi_3 = fungi_3,
  fungi_4 = fungi_4,
  bacteria_1 = bac_1,
  bacteria_2 = bac_2,
  bacteria_3 = bac_3,
  bacteria_4 = bac_4
)


library(mclust)
library(aricode)
library(clusteval)

set.seed(123)
n_iterations <- 50


fungi_genera <- colnames(fun0)
bacteria_genera <- colnames(bac0)  


fungi_groups <- list(
  fungi_1 = fungi_1,
  fungi_2 = fungi_2,
  fungi_3 = fungi_3,
  fungi_4 = fungi_4
)

bacteria_groups <- list(
  bacteria_1 = bac_1,
  bacteria_2 = bac_2,
  bacteria_3 = bac_3,
  bacteria_4 = bac_4
)

original_fungi_clusters <- c()
for(i in 1:length(fungi_groups)) {
  for(genus in fungi_groups[[i]]) {
    if(genus %in% fungi_genera) {
      original_fungi_clusters[genus] <- i
    }
  }
}

original_bacteria_clusters <- c()
for(i in 1:length(bacteria_groups)) {
  for(genus in bacteria_groups[[i]]) {
    if(genus %in% bacteria_genera) {
      original_bacteria_clusters[genus] <- i
    }
  }
}

original_row_clusters <- original_fungi_clusters
original_col_clusters <- original_bacteria_clusters


pb <- txtProgressBar(min = 0, max = n_iterations, style = 3)

bootstrap_row_clusters <- list()
bootstrap_col_clusters <- list()
bootstrap_cor_matrices <- list()

consistency_results <- data.frame(
  Iteration = 1:n_iterations,
  Row_ARI = NA,   
  Col_ARI = NA, 
  Row_Jaccard = NA  ,
  Row_RI = NA  ,
  Row_NMI = NA,
  Col_NMI = NA,
  Col_Jaccard = NA ,
  Col_RI = NA  
)

for(i in 1:n_iterations) {

  sampled_samples <- sample(sample_names, size = n_samples, replace = TRUE)
  if(any(duplicated(sampled_samples))) {
    unique_samples <- make.unique(sampled_samples)
  } else {
    unique_samples <- sampled_samples
  }
  
  bootstrap_fungi <- matrix(0, nrow = length(fungi_genera), ncol = n_samples)
  rownames(bootstrap_fungi) <- fungi_genera
  colnames(bootstrap_fungi) <- unique_samples
  
  bootstrap_bacteria <- matrix(0, nrow = length(bacteria_genera), ncol = n_samples)
  rownames(bootstrap_bacteria) <- bacteria_genera
  colnames(bootstrap_bacteria) <- unique_samples
  
  for(j in 1:n_samples) {
    original_sample <- sampled_samples[j]  
    new_sample <- unique_samples[j]      
    
    if(original_sample %in% colnames(all)) {
      bootstrap_fungi[, j] <- as.numeric(all[fungi_genera, original_sample])
      bootstrap_bacteria[, j] <- as.numeric(all[bacteria_genera, original_sample])
    }
  }
  
  bootstrap_cor <- cor(t(bootstrap_fungi), t(bootstrap_bacteria), 
                       method = "spearman", 
                       use = "pairwise.complete.obs")
  heatmap_result00 = pheatmap(bootstrap_cor)
  bootstrap_cor_matrices[[i]] <- bootstrap_cor
  
  # bootstrap_row_dist <- dist(bootstrap_cor, method = "euclidean")
  # bootstrap_row_hc <- hclust(bootstrap_row_dist, method = "complete")
  bootstrap_row_clusters[[i]] <- cutree(heatmap_result00$tree_row, k = 4) 
  
  # bootstrap_col_dist <- dist(t(bootstrap_cor), method = "euclidean")
  # bootstrap_col_hc <- hclust(bootstrap_col_dist, method = "complete")
  bootstrap_col_clusters[[i]] <- cutree(heatmap_result00$tree_col, k = 4)

  common_rows <- intersect(names(original_row_clusters), 
                           names(bootstrap_row_clusters[[i]]))
  if(length(common_rows) > 1) {
    consistency_results$Row_ARI[i] <- adjustedRandIndex(
      original_row_clusters[common_rows],
      bootstrap_row_clusters[[i]][common_rows]
    )
    consistency_results$Row_NMI[i] <- NMI(
      original_row_clusters[common_rows],
      bootstrap_row_clusters[[i]][common_rows]
    )
    consistency_results$Row_RI[i] <- RI(
      original_row_clusters[common_rows],
      bootstrap_row_clusters[[i]][common_rows]
    )
    consistency_results$Row_Jaccard[i] <- external_validation( original_row_clusters[common_rows], 
                                                               bootstrap_row_clusters[[i]][common_rows], method = "jaccard_index") 
  }

  common_cols <- intersect(names(original_col_clusters), 
                           names(bootstrap_col_clusters[[i]]))
  if(length(common_cols) > 1) {
    consistency_results$Col_ARI[i] <- adjustedRandIndex(
      original_col_clusters[common_cols],
      bootstrap_col_clusters[[i]][common_cols]
    )
    consistency_results$Col_NMI[i] <- NMI(
      original_col_clusters[common_cols],
      bootstrap_col_clusters[[i]][common_cols]
    )
    consistency_results$Col_RI[i] <- RI(
      original_col_clusters[common_cols],
      bootstrap_col_clusters[[i]][common_cols]
    )
    consistency_results$Col_Jaccard[i] <- external_validation( original_col_clusters[common_cols], 
                         bootstrap_col_clusters[[i]][common_cols], method = "jaccard_index") 
  }
  
  setTxtProgressBar(pb, i)
}
close(pb)


library(ggplot2)
library(tidyr)

plot_data <- consistency_results %>%
  dplyr::select(Row_ARI, Row_NMI, Row_Jaccard,Row_RI, Col_ARI, Col_NMI, Col_Jaccard, Col_RI) %>%
  pivot_longer(cols = everything(), 
               names_to = "Metric", 
               values_to = "Value") %>%
  mutate(
    Type = ifelse(grepl("Row", Metric), "Fungi", "Bacteria"),
    Metric = gsub("Row_|Col_", "", Metric)
  )


plot_data$Metric = factor(plot_data$Metric,levels = c('RI','NMI','ARI','Jaccard'))
p <- ggplot(plot_data, aes(x = Metric, y = Value, fill = Type)) +

  geom_boxplot(alpha = 0.8, 
               outlier.shape = 21, 
               outlier.size = 2,
               outlier.fill = "gray50",
               outlier.alpha = 0.5,
               outlier.stroke = 0.3,
               width = 0.6,
               position = position_dodge(width = 0.8),
               lwd = 0.8,
               fatten = 2) +
  

  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, 
                                              dodge.width = 0.8,
                                              jitter.height = 0),
              alpha = 0.2, 
              size = 1.2,
              aes(color = Type),
              show.legend = FALSE) +
  

  # geom_text(data = median_data,
  #           aes(x = Metric, y = median_value + 0.03, 
  #               label = sprintf("%.2f", median_value),
  #               group = Type),
  #           position = position_dodge(width = 0.8),
  #           size = 3.5,
  #           fontface = "bold") +

scale_fill_manual(values = c("Fungi" = "#4477AA", "Bacteria" = "#EE6677")) +
  scale_color_manual(values = c("Fungi" = "#4477AA", "Bacteria" = "#EE6677")) +
  

  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"),
    expand = c(0.02, 0.02) 
  ) +
  
  # geom_hline(yintercept = c(0.5, 0.7, 0.9), 
  #            linetype = "dashed", 
  #            color = "gray70", 
  #            alpha = 0.5,
  #            linewidth = 0.3) +
  

  labs(
    title = "Comparison of Clustering Consistency Metrics",
    subtitle = paste0("Based on ", n_iterations, " Bootstrap iterations"),
    x = "Metric",
    y = "Consistency Value",
    fill = "Group",
    caption = "Higher values indicate better consistency"
  ) +
  

  theme_minimal(base_size = 12) +
  theme(

    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", 
                              margin = margin(b = 5)),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40",
                                 margin = margin(b = 15)),
    plot.caption = element_text(hjust = 0.5, size = 9, color = "gray60",
                                margin = margin(t = 10)),
    
  
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 13, face = "bold", margin = margin(r = 10)),
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 11),
    axis.ticks = element_line(color = "gray70"),
    axis.ticks.length = unit(0.2, "cm"),
    

    # panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    # panel.background = element_rect(fill = "white", color = NA),
    # plot.background = element_rect(fill = "white", color = NA),
    

    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.size = unit(1, "cm"),
    legend.background = element_rect(fill = "white", color = "gray90", 
                                     linewidth = 0.3),
    legend.margin = margin(t = 10, b = 5, l = 10, r = 10),
    

    panel.border = element_rect(fill = NA, color = "gray70", linewidth = 0.5),
    

    strip.background = element_rect(fill = "gray95", color = "gray70"),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  

  # stat_compare_means(aes(group = Type), 
  #                    label = "p.signif", 
  #                    method = "wilcox.test",
  #                    hide.ns = TRUE) +
  

  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           alpha = 0.03, fill = "gray80") +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf,
           alpha = 0.03, fill = "gray80")


print(p)

print(p)



