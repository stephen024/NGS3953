
source("/Users/sunq17/Desktop/intern6-21/728a/BulkRNAseq_AnalysisFunctions.R")

#library(dplyr)
#library(ComplexHeatmap)
#library(clusterProfiler)
#library(DOSE)
library(cluster)
library(factoextra)
set.seed(999)    ## used for k means ####

es_ds=readRDS('data/NGS3953_ESET_do.Rdata')
#varLabels(es_ds)
#mat=exprs(es_ds);dim(mat)
#head(sampleNames(es_ds));length(sampleNames(es_ds))

##  combine left and right ear combo treatment to one variable Combo_all##
es_ds$Treatment=replace(es_ds$Treatment,es_ds$Treatment == "Combo_0","Combo")
es_ds$Treatment=replace(es_ds$Treatment,es_ds$Treatment == "Combo_1","Combo")

## typo in excel file ###
es_ds$Phenotype=replace(es_ds$Phenotype,es_ds$Phenotype == "infamed","inflamed")
#table(es_ds$Degree)

### 
deg_obj = DGEList(exprs(es_ds),genes=fData(es_ds),group=es_ds$Treatment)
df=data.frame(Degree=es_ds$Degree,
              Phenotype=es_ds$Phenotype)

keep.exprs = filterByExpr(deg_obj, group=es_ds$Treatment)
yf =  deg_obj[keep.exprs, keep.lib.sizes=FALSE]
yf$samples = cbind(yf$samples, df)

yf=yf[,yf$samples$Pheno!="undefined"]    ##  delete undefined pheno SAM ##

levels(yf$samples$Degree)=c("low","median","strong")
levels(yf$samples$Phenotype)=c("desert","excluded","inflamed")

yf = calcNormFactors(yf)
yf$genes$Length=yf$genes$size

ExprMat = cpm(yf,log=T) 
# rpkmMat = edgeR::rpkm(yf,log=F)
# tpmMat = tpm_from_rpkm(rpkmMat)


var_genes = apply(ExprMat, 1, var)
var_genes = var_genes[order(var_genes,decreasing = T)]
nvar = 500       
topgenes = var_genes[1:nvar]
topdat = ExprMat[names(topgenes),]

###  Notice 1 : this is important for making a ordered column split (color )##
## Treatment - Pheno (D-E-I) - Degree(L-M-H)  ##
##  will be used in Heatmap function ###
new_yf=yf$samples %>% arrange(group,Phenotype,Degree)
sam_order=rownames(new_yf)
###  Done ####

sample_cor=list(
  Treatment=c(Control="gainsboro",Combo="green",TGFB="violet",PDL1="orange"),
  Phenotype=c(desert="darkgray",excluded="dodgerblue",inflamed="deeppink"),
  Degree=c(low="light green",median="lime green",strong="dark olive green")
)

### Notice 2 : This is for column annotation , trt-pheno-degree-row names   ####
### should be 1-1 mapping,  used in heatmap. Notice 1 will control order ####

sample_col=data.frame(
  Treatment=yf$samples$group,
  Phenotype=yf$samples$Phenotype,
  Degree=yf$samples$Degree,
  row.names=rownames(yf$samples)
)
sample_col$Treatment=factor(sample_col$Treatment,levels = c("Control","Combo","TGFB","PDL1"))
sample_col$Phenotype=factor(sample_col$Phenotype,levels = c("desert","excluded","inflamed"))
sample_col$Degree=factor(sample_col$Degree,levels = c("low","median","strong"))
### Done ###

top_bar=HeatmapAnnotation(df=sample_col,col=sample_cor)

###  Determine optimal clusters #####
#group=kmeans(topdat,centers=10)

# gap_stat=clusGap(topdatt,FUN=kmeans,nstart=25,K.max=20, B=50) 

# B is too small ,rows##
## if you wanna use GAP, change B value, at least 70% of our row. Need Very long time ###

#print(gap_stat, method = "firstmax")
#fviz_gap_stat(gap_stat)

## Notice change 'nstart' value to see plot  ####
fviz_nbclust(topdat, kmeans,nstart=15,method = "silhouette")+ 
  labs(caption = "nstart=15")
fviz_nbclust(topdat, kmeans,nstart=20, method = "wss")
### Done   ####

##  Define K=10  in kmeans clustering ###
scaled_mat = t(scale(t(topdat)))

c=10
group = cluster::pam(scaled_mat, k = c)

###  skip to GO analysis ###

mouse_db="org.Mm.eg.db"

exp_matrix=topdat

#----------length is number from kmeans  #
go_result=vector(mode = "list",length=c)

#gene_list_cluster = names(id_out)[id_out == 1]

## Notice 3: this is for row annotation .You can use 'go_fun.R" ####
g_in=group$clustering

for(i in 1:c){
  #i=2
  gene_vec=names(g_in)[g_in ==i]
  gene_cluster<-bitr(gene=gene_vec,fromType="SYMBOL",toType = c("ENTREZID","ENSEMBL"),
                     OrgDb =org.Mm.eg.db )
  
  over_rep=enrichGO(gene=gene_cluster$ENTREZID,OrgDb = mouse_db,
                    ont = "ALL",pAdjustMethod = "fdr",
                    pvalueCutoff = 0.01,qvalueCutoff = 0.05,
                    readable = TRUE)
  enrich_out=data.frame(Description=over_rep@result$Description[1:5],
                        FDR_Pvalue=over_rep@result$p.adjust[1:5],
                        geneID=over_rep@result$geneID[1:5])
  go_result[[i]]=enrich_out
  
}

### save csv file ##
setwd("/Users/sunq17/Desktop/new_bulk/GO_ana")
## GO result ####
go_path=paste0(paste("GO_out","k",c,sep = "_"),".csv")
go_csv=bind_rows(go_result, .id = "cluster_name")
write.csv(go_csv,file=go_path)
  
##  heat markers for all clusters  saving ###
name_vec=names(group$clustering)
heat_df=NULL
for ( j in 1:c) {
  gene_name=name_vec[group$clustering==j]
  cluster_index=rep(j,length(gene_name))
  this_df=data.frame(cluster_index=cluster_index,markers=gene_name)
  heat_df=rbind(heat_df,this_df)
}

heat_path=paste0(paste("marker_hm_k",c,sep = "_"),".csv")
write.csv(heat_df,file=heat_path)

### below is used for heatmap row split ##
for (i in 1:c){
  
  g_in[g_in==i]=go_result[[i]]$Description[1]
}

## will change group cluster object 
#group$clustering=g_in

###  Next run Heatmap ##

heat_2km=Heatmap(scaled_mat,
                 show_column_names = TRUE,cluster_columns = FALSE,
                 show_column_dend = FALSE,
                 column_names_gp = gpar(fontsize = 3.5),
                 column_split = sample_col$Treatment,
                 #column_order = sam_order,
                 column_order = col_order_id,
                 column_title ="Top 500 Gene Expression(Bulk NGS3953)",
                 column_title_side = "top",
                 show_row_names = FALSE,show_row_dend = TRUE,
                 row_dend_side = "left",
                 row_split = g_in,
              #cluster_rows =cluster_within_group(t(scaled_mat), group$cluster),
                cluster_row_slices = TRUE,row_title_rot = 0,
              #row_title = cluster_name,
              name="s_mat",top_annotation = top_bar
)


### secretome Gene DO not run ###
heat_2km=Heatmap(scaled_mat,
                 show_column_names = TRUE,cluster_columns = FALSE,
                 show_column_dend = FALSE,
                 column_names_gp = gpar(fontsize = 3.5),
                 column_split = sample_col$Treatment,
                 column_order = sam_order,
                 # column_order = col_order_id,
                 column_title ="226 secretome Gene Expression(Bulk NGS3953)",
                 column_title_side = "top",
                 show_row_names = TRUE,show_row_dend = TRUE,
                 row_dend_side = "left",
                 row_split = g_in,
                 row_labels = rownames(scaled_mat),
                 row_names_side = "right",
                 row_names_gp = gpar(fontsize =3),
                 #cluster_rows =cluster_within_group(t(scaled_mat), group$cluster),
                 cluster_row_slices = TRUE,row_title_rot = 0,
                 #row_title = cluster_name,
                 name="s_mat",top_annotation = top_bar
)
###  Notice 4: run below before heatmap if you want to ####
### change SAMXXX to tumor ID 
##  Also change switch 'column_order' in heatmap function  ###

dim(scaled_mat)
t_id=match(colnames(scaled_mat),es_ds$Sample)
col_order_id=match(sam_order,colnames(scaled_mat))

colnames(scaled_mat)=as.numeric(es_ds$Tumor_id[t_id]) 

###   Done #####

#### read Ximo Mouse genes ###
mou_markers=read.csv(file="/Users/sunq17/Desktop/new_bulk/GO_ana/Mouse Secretome.csv",header = TRUE)
row_markers=mou_markers$Symbol
length(row_markers)
setdiff(row_markers,rownames(ExprMat))
overlap_m=intersect(row_markers,rownames(ExprMat))
mou_matrix=ExprMat[overlap_m,]
topdat=mou_matrix