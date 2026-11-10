
times = c('w','wo')
for (CST in CSTs) {
  for (time in times) {
    micro = OTUs[OTUs$Genus %in% get(CST), c('Species', get(time))]
    metab0 = metab[, c('Metabolite Name', get(time))]
    
    rownames(micro) = NULL
    micro = column_to_rownames(micro, var = 'Species')
    metab0 = column_to_rownames(metab0, var = 'Metabolite Name')
    
    metab_t <- as.data.frame(t(metab0))
    micro_t <- as.data.frame(t(micro))
    
    library(Hmisc)
    cor_results <- rcorr(as.matrix(metab_t), as.matrix(micro_t), type = "spearman")

    cor_matrix <- cor_results$r[colnames(metab_t), colnames(micro_t)]
    pval_matrix <- cor_results$P[colnames(metab_t), colnames(micro_t)]
    
    adj_pval_matrix <- matrix(p.adjust(pval_matrix, method = "BH"), 
                              nrow = nrow(pval_matrix),
                              dimnames = dimnames(pval_matrix))
    significant_pairs <- which(adj_pval_matrix < 0.05 & abs(cor_matrix) > 0.3, arr.ind = TRUE)

    if (length(significant_pairs) == 0) {
   
      next 
    }
    
    significant_metab <- rownames(cor_matrix)[significant_pairs[, 1]]
    significant_metab <- unique(significant_metab) 
    
    library(dplyr)
    significant_results <- 

      expand_grid(
        Metabolite = rownames(cor_matrix),
        Microbe = colnames(cor_matrix)
      ) %>%

      mutate(
        Correlation = map2_dbl(Metabolite, Microbe, ~cor_matrix[.x, .y]),
        Pvalue = map2_dbl(Metabolite, Microbe, ~pval_matrix[.x, .y]),
        AdjPvalue = map2_dbl(Metabolite, Microbe, ~adj_pval_matrix[.x, .y])
      ) %>%

      filter(AdjPvalue < 0.05 & abs(Correlation) > 0.3) %>%

      arrange(desc(abs(Correlation)))
    ids = merged_df[merged_df$Metabolite %in% significant_results$Metabolite,]
    

    kegg <- enricher(gene = ids$KEGGid,
                     TERM2GENE = com_ano[c('pathway_id', 'KEGG.ID')], 
                     TERM2NAME = com_ano[c('pathway_id', 'pathway_name')],
                     pvalueCutoff = 0.05, 
                     pAdjustMethod = 'BH', 
                     qvalueCutoff = 0.2, 
                     maxGSSize = 500)
    
    result_df <- as.data.frame(kegg)
    
    if (nrow(result_df) == 0) {
    } else {
      result_df$time = time
      result_df$CST = CST
      result_df$change = 'up'
    }
    

    if (nrow(result_df) > 0) {
      all_results <- rbind(all_results, result_df)
    }
  }
}


#----------------------------------------
all_results = all_results[,-15]
data = all_results[,c(2,3,9,12,13,14)]


data$p.adjust <- as.numeric(data$p.adjust)
data$GeneRatio_num <- as.numeric(gsub("/.*", "", data$GeneRatio)) 

filtered_df <- data %>%
  group_by(CST, Description) %>%
  summarise(n_time = n_distinct(time), .groups = "drop") %>%
  filter(n_time == 1)


df_filtered <- data %>%
  inner_join(filtered_df, by = c("CST", "Description"))
df_filtered$time = factor(df_filtered$time,levels = c('wo','w'))
for (i in unique(df_filtered$CST)){
  sub = df_filtered[df_filtered$CST == i,]
  ggplot(sub, 
         aes(x = time, 
             y = Description, 
             color = p.adjust, 
             size = GeneRatio_num)) +
    geom_point() +
    scale_color_continuous(
      low = "#EF767B", 
      high = "#BFBBBA") +  
    scale_size_continuous(
      range = c(2, 10)) +

    theme_minimal() +
    labs(title = "Bubble Plot of KEGG Pathways", 
         x = "CST", 
         y = "Description", 
         color = "p.adjust", 
         size = "Gene Ratio") +
    theme(
      panel.grid = element_blank(),            
      panel.border = element_rect(           
        colour = "black", 
        fill = NA, 
        linewidth = 0.5),
      strip.background = element_rect(           
        fill = "#f0f0f0", 
        colour = "black"),
      strip.text = element_text(face = "bold"),  
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}
