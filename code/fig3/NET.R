library(microeco)
library(fdrtool)
library(meconetcomp)
library(magrittr)
library(igraph)
library(ggplot2)

fungi_1 = c('Ramularia','Trichosporon','Clavispora','Diutina',
            'Malassezia','Geosmithia',
            'Pichia','Neurpspora','Pyricularia')
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

names = c(fun_names,bac_names)
df = df[df$Genus %in% names,]
rownames(df) = df$Species
df = df[,-c(1,2,3)]
rownames(df) = gsub('-','_',rownames(df))
patterns <- c("\\[", "\\]", "\\(", "\\)", "\\.", "\\'", "\\/", "\\:", "\\=", "\\+")
rownames(df) <- Reduce(function(x, pattern) gsub(pattern, "", x), patterns, rownames(df) )

info = rbind(bac,fun)
info$fdr = p.adjust(info$kw_pvalue)
info = info[info$key %in% rownames(df),]
info = info[info$fdr < 0.01,]
df = df[rownames(df)%in%info$key,]


tax = read.delim('./data/tax.txt')
rownames(tax) = tax$Species
rownames(tax) = gsub('-','_',rownames(tax))
patterns <- c("\\[", "\\]", "\\(", "\\)", "\\.", "\\'", "\\/", "\\:", "\\=", "\\+")
rownames(tax) <- Reduce(function(x, pattern) gsub(pattern, "", x), patterns, rownames(tax) )
tax = tax[rownames(tax) %in% rownames(df),]
tax = tidy_taxonomy(tax)



mt <- microtable$new(sample_table = group, otu_table = df, tax_table = tax)

mt


networks <- list()

# W_T0


tmp <-clone(mt)
tmp$sample_table %<>%subset(Group == "W_T0")
tmp$tidy_dataset()
#tmp$filter_taxa(freq = 0.5)
tmp <- trans_network$new(dataset = tmp, cor_method ="spearman",use_WGCNA_pearson_spearman = T,nThreads = 10)
tmp$cal_network(COR_p_thres =0.01, COR_cut =0.6)
networks$W_T0 <- tmp

# W_T1

tmp <-clone(mt)
tmp$sample_table %<>%subset(Group == "W_T1")
tmp$tidy_dataset()
#tmp$filter_taxa(freq = 0.5)
tmp <- trans_network$new(dataset = tmp, cor_method ="spearman",use_WGCNA_pearson_spearman = T,nThreads = 10)
tmp$cal_network(COR_p_thres =0.01, COR_cut =0.6)
networks$W_T1 <- tmp

# W_T2

tmp <-clone(mt)
tmp$sample_table %<>%subset(Group == "W_T2")
tmp$tidy_dataset()
#tmp$filter_taxa(freq = 0.5)
tmp <- trans_network$new(dataset = tmp, cor_method ="spearman",use_WGCNA_pearson_spearman = T,nThreads = 10)
tmp$cal_network(COR_p_thres =0.01, COR_cut =0.6)
networks$W_T2 <- tmp

# W_T3

tmp <-clone(mt)
tmp$sample_table %<>%subset(Group == "W_T3")
tmp$tidy_dataset()
#tmp$filter_taxa(freq = 0.5)
tmp <- trans_network$new(dataset = tmp, cor_method ="spearman",use_WGCNA_pearson_spearman = T,nThreads = 10)
tmp$cal_network(COR_p_thres =0.01, COR_cut =0.6)
networks$W_T3 <- tmp

# WO_T0

tmp <-clone(mt)
tmp$sample_table %<>%subset(Group == "WO_T0")
tmp$tidy_dataset()
#tmp$filter_taxa(freq = 0.5)
tmp <- trans_network$new(dataset = tmp, cor_method ="spearman",use_WGCNA_pearson_spearman = T,nThreads = 10)
tmp$cal_network(COR_p_thres =0.01, COR_cut =0.6)
networks$WO_T0 <- tmp

# WO_T1

tmp <-clone(mt)
tmp$sample_table %<>%subset(Group == "WO_T1")
tmp$tidy_dataset()
#tmp$filter_taxa(freq = 0.5)
tmp <- trans_network$new(dataset = tmp, cor_method ="spearman",use_WGCNA_pearson_spearman = T,nThreads = 10)
tmp$cal_network(COR_p_thres =0.01, COR_cut =0.6)
networks$WO_T1 <- tmp

# WO_T3

tmp <-clone(mt)
tmp$sample_table %<>%subset(Group == "WO_T3")
tmp$tidy_dataset()
#tmp$filter_taxa(freq = 0.5)
tmp <- trans_network$new(dataset = tmp, cor_method ="spearman",use_WGCNA_pearson_spearman = T,nThreads = 10)
tmp$cal_network(COR_p_thres =0.01, COR_cut =0.6)
networks$WO_T3 <- tmp

#--------------------------------------------------------------------------------
fungi_1 = c('Ramularia','Trichosporon','Clavispora','Diutina',
            'Malassezia','Geosmithia',
            'Pichia','Neurpspora','Pyricularia')
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
vectors <- list(
  fungi_1 = fungi_1,
  fungi_2 = fungi_2,
  fungi_3 = fungi_3,
  fungi_4 = fungi_4,
  bac_1 = bac_1,
  bac_2 = bac_2,
  bac_3 = bac_3,
  bac_4 = bac_4
)
class <- do.call(rbind, lapply(names(vectors), function(name) {
  data.frame(Element = vectors[[name]], Source = name, stringsAsFactors = FALSE)
}))
class$Element = paste0('g__',class$Element)
colnames(class) = c('Genus','Module')
#-----------------------------------------------------------------------------------------


networks %<>%cal_module(undirected_method ="cluster_fast_greedy")
tmp <- cal_network_attr(networks)

networks %<>% get_node_table(node_roles = TRUE) %>% get_edge_table

color = data.frame(module = c('bac_1','bac_2','bac_3','bac_4',
                              'fungi_1','fungi_2','fungi_3','fungi_4'))
#extended_palette <- colorRampPalette(brewer.pal(8, "Paired"))
color$color = c("#A51C36", "#7ABBDB", "#84BA42", "#682487", "#DBB428", 
                "#D4562E", "#FFFFB3",'#FCCDE5')

module_colors <- setNames(color$color, color$module)



for(i in names(networks)){
  #i='W_T0'
  nodes = networks[[i]]$res_node_table
  edges = networks[[i]]$res_edge_table
  nodes = merge(nodes,class,by = 'Genus')
  nodes = nodes[,c(2,1,3:ncol(nodes))]
  abun = as.data.frame(t(networks[[i]]$data_abund))
  abun = abun[rownames(abun)%in%nodes$name,]
  nodes$abun = rowSums(abun)
  nodes <- nodes %>% select(-module)
  nodes = merge(nodes,color,by.x = 'Module',by.y = 'module')
  nodes = nodes[,c(2,1,3:ncol(nodes))]
  colnames(nodes)[1] = 'Id'
  colnames(edges)[1:2] = c('source','target')
  write.table(edges, paste0('./fig2/1229/net/',i,'_network.edge_list.csv'), sep = ',', row.names = FALSE, quote = FALSE)
  write.table(nodes, paste0('./fig2/1229/net/',i,'_network.node_list.csv'), sep = ',', row.names = FALSE, quote = FALSE)
  #networks[[i]]$save_network(filepath =paste0(tmp, "/", i, ".gexf"))
  
}






save(networks,file = './fig2/1229/net/networks.RData')
networks[["W_T0"]]$plot_network(method ="networkD3", node_color ="Kingdom")
networks[["W_T1"]]$plot_network(method ="networkD3", node_color ="Kingdom")
networks[["W_T2"]]$plot_network(method ="networkD3", node_color ="Kingdom")
networks[["W_T3"]]$plot_network(method ="networkD3", node_color ="Kingdom")
networks[["WO_T0"]]$plot_network(method ="networkD3", node_color ="Kingdom")
networks[["WO_T1"]]$plot_network(method ="networkD3", node_color ="Kingdom")
networks[["WO_T3"]]$plot_network(method ="networkD3", node_color ="Kingdom")



#------------------------------网络子集---------------------------------------------


