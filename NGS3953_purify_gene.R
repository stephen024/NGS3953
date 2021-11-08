


source("/Users/sunq17/Desktop/intern6-21/728a/BulkRNAseq_AnalysisFunctions.R")

es_ds=readRDS('/Users/sunq17/Desktop/new_bulk/data/NGS3953_ESET_do.Rdata')

##  combine left and right ear combo treatment to one variable Combo_all##
es_ds$Treatment=replace(es_ds$Treatment,es_ds$Treatment == "Combo_0","Combo")
es_ds$Treatment=replace(es_ds$Treatment,es_ds$Treatment == "Combo_1","Combo")

## typo in excel file ###
es_ds$Phenotype=replace(es_ds$Phenotype,es_ds$Phenotype == "infamed","inflamed")

### 
deg_obj = DGEList(exprs(es_ds),genes=fData(es_ds),group=es_ds$Treatment)
df=data.frame(Degree=es_ds$Degree,Phenotype=es_ds$Phenotype)

keep.exprs = filterByExpr(deg_obj, group=es_ds$Treatment)
yf =  deg_obj[keep.exprs, keep.lib.sizes=FALSE]
yf$samples = cbind(yf$samples, df)

yf=yf[,yf$samples$Pheno!="undefined"] 


yf = calcNormFactors(yf)
yf$genes$Length=yf$genes$size

ExprMat = cpm(yf,log=T) 

intersect(rownames(ExprMat),"Trdv5")  same for Tcrg-v4

### delete below ###
var_genes = apply(ExprMat, 1, var)
var_genes = var_genes[order(var_genes,decreasing = T)]
nvar = 500       
topgenes = var_genes[1:nvar]
topdat = ExprMat[names(topgenes),]
############# not run above ###

# mark_vec=unlist(gene_list,use.names = FALSE)

###  for each cell type // cluster ####

### cdc1  ,  gdtcekk   no boxplot for small number of markers  ###
for (k in 20:length(gene_list)){
   k=5
   cell_in=names(gene_list)[k]
   mark_vec=unlist(gene_list[k],use.names = FALSE)
   Ifn_jeremy = ExprMat[(rownames(ExprMat) %in% mark_vec),]  ## can not find markers "Trdv5" and "Tcrg-V4"
   
   # G-Score
   pcSig = gsScore(Ifn_jeremy)     ## can not run for one marker matrix ## 
   yf$samples$pcSigJeremie = pcSig
   scores = yf$samples
   # head(yf$samples)
   
   setwd('/Users/sunq17/Desktop/Oct/purify_box')
   png(paste(cell_in,"box.png",sep="_"), units="px", width=597, height=441)
   print(
     ggboxplot(scores, x = "group", y = "pcSigJeremie",color = "Phenotype",add="jitter")+
       ylab("Gene Score")+
       labs(title = paste("Expression of",cell_in))
   )
   dev.off()
}


# G-Score
pcSig = gsScore(Ifn_jeremy)
yf$samples$pcSigJeremie = pcSig
scores = yf$samples
head(yf$samples)

ggboxplot(scores, x = "group", y = "pcSigJeremie",color = "Phenotype",add="jitter")+
   #ylim(15, 20)
  ylab("Gene Score")
 

####
#new_yf=yf$samples %>% arrange(group,Phenotype,Degree)
#sam_order=rownames(new_yf)