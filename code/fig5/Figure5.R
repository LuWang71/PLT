# 初始化一个空数据框用于存储所有结果
all_results <- data.frame()
merged_df =  read.delim('./nm/keggid.csv',sep=',')
groups = read.delim('./data/metadata2.txt',sep='\t')
groups$Timepoint0 <- ifelse(groups$Group == "w/",
                            paste0("w_", groups$Timepoint),
                            ifelse(groups$Group == "w/o",
                                   paste0("wo_", groups$Timepoint),
                                   groups$Timepoint))
groups = groups[groups$Timepoint == 'T0',]
w = intersect(groups[groups$Group == 'w/',]$SampleID,colnames(metab))
wo = intersect(groups[groups$Group == 'w/o',]$SampleID,colnames(metab))
#groups = groups[,c(1,7)]
metab = read_excel("F:/projects/肝移植/imeta/二审/稿子/ori_filtered_KNN.xlsx")
metab = metab[,-c(2:10)]
times = c('w','wo')
for (CST in CSTs) {
  for (time in times) {
     # CST='bac_1'
     # time = 'w'
    # 动态提取变量并执行相关操作
    micro = OTUs[OTUs$Genus %in% get(CST), c('Species', get(time))]
    metab0 = metab[, c('Metabolite Name', get(time))]
    
    rownames(micro) = NULL
    micro = column_to_rownames(micro, var = 'Species')
    metab0 = column_to_rownames(metab0, var = 'Metabolite Name')
    
    metab_t <- as.data.frame(t(metab0))
    micro_t <- as.data.frame(t(micro))
    
    library(Hmisc)
    cor_results <- rcorr(as.matrix(metab_t), as.matrix(micro_t), type = "spearman")
    
    # 提取相关系数和p值
    cor_matrix <- cor_results$r[colnames(metab_t), colnames(micro_t)]
    pval_matrix <- cor_results$P[colnames(metab_t), colnames(micro_t)]
    
    adj_pval_matrix <- matrix(p.adjust(pval_matrix, method = "BH"), 
                              nrow = nrow(pval_matrix),
                              dimnames = dimnames(pval_matrix))
    
    # 筛选显著关联（FDR < 0.05且|r| > 0.3）
    significant_pairs <- which(adj_pval_matrix < 0.05 & abs(cor_matrix) > 0.3, arr.ind = TRUE)
    
    # 如果没有显著性结果，输出提示
    if (length(significant_pairs) == 0) {
      cat(paste("CST：", CST, '在', time, "时无显著性\n"))
      next  # 跳过此次循环，进入下一个CST和时间点组合
    }
    
    significant_metab <- rownames(cor_matrix)[significant_pairs[, 1]]
    significant_metab <- unique(significant_metab)  # 去重
    
    library(dplyr)
    significant_results <- 
      # 创建所有可能组合的数据框
      expand_grid(
        Metabolite = rownames(cor_matrix),
        Microbe = colnames(cor_matrix)
      ) %>%
      # 添加相关系数和校正后p值
      mutate(
        Correlation = map2_dbl(Metabolite, Microbe, ~cor_matrix[.x, .y]),
        Pvalue = map2_dbl(Metabolite, Microbe, ~pval_matrix[.x, .y]),
        AdjPvalue = map2_dbl(Metabolite, Microbe, ~adj_pval_matrix[.x, .y])
      ) %>%
      # 筛选显著结果
      filter(AdjPvalue < 0.05 & abs(Correlation) > 0.3) %>%
      # 按相关性排序
      arrange(desc(abs(Correlation)))
    ids = merged_df[merged_df$Metabolite %in% significant_results$Metabolite,]
    
    # 进行 KEGG 富集分析
    kegg <- enricher(gene = ids$KEGGid,
                     TERM2GENE = com_ano[c('pathway_id', 'KEGG.ID')], 
                     TERM2NAME = com_ano[c('pathway_id', 'pathway_name')],
                     pvalueCutoff = 0.05, 
                     pAdjustMethod = 'BH', 
                     qvalueCutoff = 0.2, 
                     maxGSSize = 500)
    
    result_df <- as.data.frame(kegg)
    
    if (nrow(result_df) == 0) {
      cat(paste("在", time, "下", CST, "没有up通路\n"))
    } else {
      result_df$time = time
      result_df$CST = CST
      result_df$change = 'up'
    }
    
    # 如果有结果，将其添加到总结果数据框
    if (nrow(result_df) > 0) {
      all_results <- rbind(all_results, result_df)
    }
  }
}

# 最终输出
write.table(all_results, './imeta/二审/combined_results_all.csv', sep = ',', row.names = FALSE)
#----------------------------------------
all_results = all_results[,-15]
data = all_results[,c(2,3,9,12,13,14)]


data$p.adjust <- as.numeric(data$p.adjust)
data$GeneRatio_num <- as.numeric(gsub("/.*", "", data$GeneRatio))  # 提取GeneRatio中的数字

filtered_df <- data %>%
  group_by(CST, Description) %>%
  summarise(n_time = n_distinct(time), .groups = "drop") %>%
  filter(n_time == 1)  # 只保留只出现在一个 time 中的通路

# 第二步：筛选原始数据中符合条件的行
df_filtered <- data %>%
  inner_join(filtered_df, by = c("CST", "Description"))
df_filtered$time = factor(df_filtered$time,levels = c('wo','w'))
for (i in unique(df_filtered$CST)){
 # i= unique(df_filtered$CST)[1]
  sub = df_filtered[df_filtered$CST == i,]
  ggplot(sub, 
         aes(x = time, 
             y = Description, 
             color = p.adjust, 
             size = GeneRatio_num)) +
    geom_point() +
    scale_color_continuous(
      low = "#EF767B", 
      high = "#BFBBBA") +  # 设置颜色渐变
    scale_size_continuous(
      range = c(2, 10)) +  # 设置气泡大小范围
    #facet_wrap(~CST, ncol = 4) +  # 每行4个分面
    theme_minimal() +
    labs(title = "Bubble Plot of KEGG Pathways", 
         x = "CST", 
         y = "Description", 
         color = "p.adjust", 
         size = "Gene Ratio") +
    theme(
      panel.grid = element_blank(),               # 移除背景网格线
      panel.border = element_rect(                # 添加每个分面的边框
        colour = "black", 
        fill = NA, 
        linewidth = 0.5),
      strip.background = element_rect(            # 分面标签背景
        fill = "#f0f0f0", 
        colour = "black"),
      strip.text = element_text(face = "bold"),   # 分面标签字体加粗
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  ggsave(paste0('./imeta/二审/稿子/kegg/',i,'.KEGG.pdf'),width = 6,height = 6.3)
}
