library(UpSetR)

data = rbind(data0,data1)
data$group = c(rep('fungi',nrow(data0)),rep('bacteria',nrow(data1)))

tmp = data[,c(1,17:25)]
info = tmp
numeric_matrix <- as.data.frame(lapply(info[,-1], function(col) ifelse(col < 0.05, 1, 0)))
rownames(numeric_matrix) <- info$key

upset(
  numeric_matrix,
  sets = colnames(numeric_matrix),
  keep.order = TRUE,
  #nintersects=50,
  order.by = "freq",
 
)
library(tibble)
data = rownames_to_column(as.data.frame(numeric_matrix),var = 'Species')
write.table(data,'./figure/suppdata/fig2/upsetR.csv',sep=',',row.names = F)
upset(
  numeric_matrix,
  sets = colnames(numeric_matrix),  
  order.by = "freq",           

)

intersection1 <- rownames(numeric_matrix[numeric_matrix[, "p_w_wo_t0"] == 1 
                                            & numeric_matrix[, "p_w_wo_t1"] == 1
                                            & numeric_matrix[, "p_w_wo_t21"] == 1
                                            & numeric_matrix[, "p_w_wo_t3"] == 1
                                         & numeric_matrix[, "p_wt01"] == 0
                                         & numeric_matrix[, "p_wt12"] == 0
                                         & numeric_matrix[, "p_wt23"] == 0
                                         & numeric_matrix[, "p_wot01"] == 0
                                         & numeric_matrix[, "p_wot13"] == 0
                                         , ])

library(pheatmap)
map = data[data$key %in% intersection1,]
rownames(map) = map$key
map = map[,c(10:16)]

map = na.omit(map)
map = map[rowSums(map[1:7])>0,]
aaa = map
aaa$sum = rowSums(aaa)
aaa = aaa[order(aaa$sum,decreasing = T),]
pheatmap(aaa[1:4,1:7],cluster_cols = F,scale ='row',gaps_col = 4,
         cellwidth = 10,cellheight = 10,
         color = c(colorRampPalette(colors = c('#2171B5',"white"))(50),colorRampPalette(colors = c("white" ,"#D94801"))(50)))

name = c('Megamonas_funiformis',
'Mediterraneibacter_gnavus',
'Alistipes_onderdonkii',
'Faecalibacterium_prausnitzii',
'Bacteroides_ovatus',
'Klebsiella_aerogenes',
'Bifidobacterium_breve',
'Hungatella_hathewayi',
'Lacticaseibacillus_paracasei',
'Streptococcus_lactarius',
'Enterococcus_faecalis',
'Bifidobacterium_pseudocatenulatum',
'Salinibaculum_sp_SYNS191',
'Halomicroarcula_marina',
'Aspergillus_viridinutans',
'Scheffersomyces_stipitis',
'Aspergillus_nomiae',
'Neurospora_tetrasperma')
info00 = info[info$key %in%name,]
