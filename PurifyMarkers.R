##################################################################################
######  Deconvolution Markers  by Ximo Pechuan i Jorge 7/10/2021
#################################################################################

library(Seurat)
####
library(ggplot2)
library(dplyr)
library(ggdendro)
library(ComplexHeatmap)
library(ggpubr)
library(tidyverse)
#markers = read.csv("../ATLASClusterMarkers.csv")
#seu = readRDS("../Final_STAMP_September21_Atlas.Rdata")

markers = read.csv("./ATLASClusterMarkers.csv")     ##  use  "Nk:Ncr1" 
seu = readRDS("./Final_STAMP_September21_Atlas.Rdata")   

# Average expresson
seu = NormalizeData(seu, normalization.method = "LogNormalize",assay = "RNA")
ExpressioMitja = AverageExpression(seu,assays = "RNA")
##  use "Ncr1_Nk" 

## For each cluster

name_cluster=unique(markers$cluster)

## check name ###
setdiff(name_cluster,colnames(ExpressioMitja$RNA))


name_cluster[c(6,8,9,10,21)]
i=5

filtered =  markers %>%  filter(cluster == name_cluster[i])
topMarkers  = filtered %>% filter(p_val_adj<0.05 & pct.1>0.5 & pct.2<0.1 & avg_log2FC>0.5 )
topMarkers  = topMarkers %>% top_n(p_val_adj,n=-1) %>% top_n(avg_log2FC,n=30)


top10 = ExpressioMitja$RNA[unique(topMarkers$gene),] %>% as.matrix()
Stop10 = scale(t(top10))
rownames(Stop10)=name_cluster

nokpp = Stop10[name_cluster[-i],]

# run this besides name_cluster[c(6,8,9,10,21)] ##

genes = names(colSums(nokpp > 0)[colSums(nokpp > 0) == 0])
## run this for name_cluster[c(6,8,9,10,21)]  ##
#genes=unique(topMarkers$gene)

top_10 = ExpressioMitja$RNA[unique(genes),] %>% as.matrix()
Finalmatrix =  scale(t(top_10))
dim(Finalmatrix)
#library(circlize)
#mat_pal =colorRamp2(c(-2, 0, 2), c("steelblue1", "white", "tomato"))

Heatmap(t(Finalmatrix),col = mat_pal,cluster_columns = T,cluster_rows = T,
        column_names_gp = grid::gpar(fontsize = 10),
        row_names_gp = grid::gpar(fontsize = 7),
        column_title = paste(as.character(i),name_cluster[i],sep = ":"),
              
        #row_split = 7,
        #column_split = 5,
        border_gp = gpar(col = "gray", lty = 2),
        row_dend_gp = gpar(col = "gray"))

#load("gene_seek.RData")
#gene_list=vector(mode = "list",length=25)  ## length of cluster ##

gene_list[[i]]=genes

save(gene_list,file="gene_seek.RData")   ## gene list ###

names(gene_list)=name_cluster


###  save data frame , output csv file #### 

#gene_df=data.frame(cluster=NUlL,gene=NULL,weight=)



gene_df = data.frame(cluster=character(),gene = character(), weight = numeric())

for (j in 1:length(gene_list)){
        sub_list=gene_list[[j]]
        df2=markers %>% filter(cluster == names(gene_list)[j]) %>% filter(gene %in% sub_list)
        gene_df= rbind(gene_df,df2)
}


getwd()
write.csv(gene_df,"sig_purify_genes.csv")

colnames(gene_df)

use_bis=gene_df[c('gene','cluster','avg_log2FC')]
