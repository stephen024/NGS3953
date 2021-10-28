
source("/Users/sunq17/Desktop/new_bulk/BulkRNAseq_AnalysisFunctions.R")

#library(dplyr)
#library(ComplexHeatmap)
#library(clusterProfiler)
#library(DOSE)
library(cluster)
library(factoextra)
set.seed(999)    ## used for k means ####

es_ds=readRDS('data/NGS3953_ESET_do.Rdata')

##  combine left and right ear combo treatment to one variable Combo_all##
es_ds$Treatment=replace(es_ds$Treatment,es_ds$Treatment == "Combo_0","Combo")
es_ds$Treatment=replace(es_ds$Treatment,es_ds$Treatment == "Combo_1","Combo")
## typo in excel file ###
es_ds$Phenotype=replace(es_ds$Phenotype,es_ds$Phenotype == "infamed","inflamed")

mouse_db="org.Mm.eg.db"
name_trt=c("Combo","Control","PDL1","TGFB")
num_marker=c(400,200,300,300)

col_list=vector(mode = "list",length=4)
sam_list=vector(mode = "list",length=4)
smat_list=vector(mode = "list",length=4)
go_list=vector(mode = "list",length=4)

cluster_in_list=vector(mode = "list",length=4)
  
for (k in 1:4){
  k=3
  trt_sel=name_trt[k]
  trt_ds=es_ds[,es_ds$Treatment==trt_sel]
  
  deg_obj = DGEList(exprs(trt_ds),genes=fData(trt_ds),group = trt_ds$Treatment)
  df=data.frame(Degree=trt_ds$Degree,
                Phenotype=trt_ds$Phenotype)
  
  keep.exprs = filterByExpr(deg_obj)
  yf =  deg_obj[keep.exprs, keep.lib.sizes=FALSE]
  yf$samples = cbind(yf$samples, df)
  
  yf=yf[,yf$samples$Phenotype!="undefined"]  
  levels(yf$samples$Degree)=c("low","median","strong")
  levels(yf$samples$Phenotype)=c("desert","excluded","inflamed")
  yf = calcNormFactors(yf)
  yf$genes$Length=yf$genes$size
  
  #ExprMat = cpm(yf,log=T) 
  rpkmMat = edgeR::rpkm(yf,log=F)
  tpmMat = tpm_from_rpkm(rpkmMat)
  
  var_genes = apply(tpmMat, 1, var)
  var_genes = var_genes[order(var_genes,decreasing = T)]
  nvar =  num_marker[k]     
  topgenes = var_genes[1:nvar]
  topdat =tpmMat[names(topgenes),]                ## No need to save ##
  
  new_yf=yf$samples %>% arrange(group,Phenotype,Degree)
  sam_order=rownames(new_yf)                   ## save 
  
  sam_list[[k]]=sam_order
  
  sample_col=data.frame(                       ## save 
    Treatment=yf$samples$group,
    Phenotype=yf$samples$Phenotype,
    Degree=yf$samples$Degree,
    row.names=rownames(yf$samples)
  )
  sample_col$Phenotype=factor(sample_col$Phenotype,levels = c("desert","excluded","inflamed"))
  sample_col$Degree=factor(sample_col$Degree,levels = c("low","median","strong"))
  
  col_list[[k]]=sample_col
  exp_matrix=topdat
  
  eout=enrich_out(expmat=exp_matrix,nk=10,scale = FALSE,go_db="org.Mm.eg.db",
                  top_go =1,pcut = 0.05)    ## save 
  
  go_list[[k]]=eout$go_result
  cluster_in_list[[k]]=eout$marker_label
  
  setwd("/Users/sunq17/Desktop/new_bulk")
  dir.create(path=paste("/Users/sunq17/Desktop/new_bulk",trt_sel,sep="/"))
  setwd(paste("/Users/sunq17/Desktop/new_bulk",trt_sel,sep="/"))
  
  go_enrich_save(enrich_list = eout$go_result)
  
  setwd("/Users/sunq17/Desktop/new_bulk")
  
  scaled_mat = t(scale(t(topdat)))
  
  t_id=match(colnames(scaled_mat),trt_ds$Sample)
  col_order_id=match(sam_order,colnames(scaled_mat))
  
  colnames(scaled_mat)=as.numeric(trt_ds$Tumor_id[t_id]) 
  ### scaled_mat ### save ######
  smat_list[[k]]=scaled_mat
}




### 


#name_trt  =c("Combo","Control","PDL1","TGFB")
#num_marker=c(400,200,300,300)

col_list 
sam_list
smat_list 

sample_cor=list(
  Treatment=c(Control="gainsboro",Combo="green",TGFB="violet",PDL1="orange"),
  Phenotype=c(desert="darkgray",excluded="dodgerblue",inflamed="deeppink"),
  Degree=c(low="light green",median="lime green",strong="dark olive green")
)


top_bar=HeatmapAnnotation(df=sample_col,col=sample_cor)

g_in=add_go_heat(eout$marker_label,eout$go_result)

###  Next run Heatmap ##

heat_2km=Heatmap(scaled_mat,
                 show_column_names = TRUE,cluster_columns = FALSE,
                 show_column_dend = FALSE,
                 column_names_gp = gpar(fontsize =4),
                 column_split = sample_col$Phenotype,
                 #column_order = sam_order,
                 column_order = col_order_id,
                 column_title ="Top Gene Expression(Bulk NGS3953)",
                 column_title_side = "top",
                 
                 show_row_names = TRUE,show_row_dend = TRUE,
                 row_dend_side = "left",
                 row_split = g_in,
                 row_labels = rownames(scaled_mat),
                 row_names_side = "right",
                 row_names_gp = gpar(fontsize =3),
                 
                cluster_row_slices = TRUE,row_title_rot = 0,
               
                name="s_mat",top_annotation = top_bar
)

