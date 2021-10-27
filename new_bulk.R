

setwd("/Users/sunq17/Desktop/new_bulk/")

source("/Users/sunq17/Desktop/intern6-21/728a/BulkRNAseq_AnalysisFunctions.R")
library(readxl)
library(tidyverse)


es_ds=readRDS('data/NGS3953_ESET.Rdata')

# sample_inventory=data.frame(x=es_ds@phenoData@data$INVENTORY_SAMPLE_NAME )

## Add phenotype information @@  ## 
es_ds@phenoData@data=es_ds@phenoData@data %>% separate(INVENTORY_SAMPLE_NAME,into=c("Tumor_id","Rest "),sep = "_Tumor_")
pheno_type=c()
pheno_level=c()

m12_con_phe=read_excel('treatment bulkRNAseq manual classification for ximo.xlsx',
                       range = "C3:E55",col_names = FALSE)
colnames(m12_con_phe)=colnames(m12_con_phe) %>% str_replace("...", "X") 

m12_con_phe=m12_con_phe %>%replace_na(list(X3="median"))
table(m12_con_phe$X3)

m12_con_phe=m12_con_phe %>% mutate (X2=replace(X2, X2 == "excluded/inflamed","undefined")) 
m12_con_phe=m12_con_phe %>% mutate (X2=replace(X2, X2 == "undetermined","undefined")) 

table(m12_con_phe$X2)

m12_total=nrow(m12_con_phe) 
if (identical(es_ds$Tumor_id[1:m12_total],as.character(m12_con_phe[[1]])))
{
  pheno_in=c(pheno_in,m12_con_phe$X2)
  pheno_level= c(pheno_level,m12_con_phe$X3)
}else{
  print("Error: Mouse12 control ")
}

### Mouse 9 R -> 0  combo ## 
m9_0phe=read_excel('treatment bulkRNAseq manual classification for ximo.xlsx',
                       range = "G3:I22",col_names = FALSE)
colnames(m9_0phe)=colnames(m9_0phe) %>% str_replace("...", "X") 

m9_0phe=m9_0phe %>%replace_na(list(X3="median"))
table(m9_0phe$X3)

#m12_con_phe=m12_con_phe %>% mutate (X2=replace(X2, X2 == "excluded/inflamed","undefined")) 
m9_0phe=m9_0phe %>% mutate (X2=replace(X2, X2 == "undetermined","undefined")) 

table(m9_0phe$X2)

m9_0total=nrow(m9_0phe) 
if (identical(es_ds$Tumor_id[(m12_total+1):(m12_total+m9_0total)],as.character(m9_0phe[[1]])))
{
  pheno_in=c(pheno_in,m9_0phe$X2)
  pheno_level= c(pheno_level,m9_0phe$X3)
}else{
  print("Error: Mouse9 Right Combo ")
}

### Mouse 9 L -> 1  combo ## 
m9_1phe=read_excel('treatment bulkRNAseq manual classification for ximo.xlsx',
                   range = "K3:M10",col_names = FALSE)
colnames(m9_1phe)=colnames(m9_1phe) %>% str_replace("...", "X") 

m9_1phe=m9_1phe %>%replace_na(list(X3="median"))
table(m9_1phe$X3)

m9_1phe=m9_1phe %>% mutate (X2=replace(X2, X2 == "undetermined","undefined")) 

table(m9_1phe$X2)

m9_1total=nrow(m9_1phe) 
beg=m12_total+m9_0total
end_ind=beg+m9_1total
if (identical(es_ds$Tumor_id[(beg+1):end_ind],as.character(m9_1phe[[1]])))
{
  pheno_in=c(pheno_in,m9_1phe$X2)
  pheno_level= c(pheno_level,m9_1phe$X3)
}else{
  print("Error: Mouse9 Left Combo ")
}

### Mouse 15 -> tgfb ## 
m15_phe=read_excel('treatment bulkRNAseq manual classification for ximo.xlsx',
                  range = "S3:U32",col_names = FALSE)
colnames(m15_phe)=colnames(m15_phe) %>% str_replace("...", "X") 

m15_phe=m15_phe %>%replace_na(list(X3="median"))
table(m15_phe$X3)
table(m15_phe$X2)

m15_total=nrow(m15_phe) 
beg=m12_total+m9_0total+m9_1total
end_ind=beg+m15_total

if (identical(es_ds$Tumor_id[(beg+1):end_ind],as.character(m15_phe[[1]])))
{
  pheno_in=c(pheno_in,m15_phe$X2)
  pheno_level= c(pheno_level,m15_phe$X3)
}else{
  print("Error: Mouse15 tgfb ")
}

### Mouse 8 L  -> pdl1 ## 
m8_phe=read_excel('treatment bulkRNAseq manual classification for ximo.xlsx',
                   range = "O3:Q31",col_names = FALSE)
colnames(m8_phe)=colnames(m8_phe) %>% str_replace("...", "X") 

m8_phe=m8_phe %>%replace_na(list(X3="median"))
table(m8_phe$X3)
table(m8_phe$X2)

m8_total=nrow(m8_phe) 
beg=m12_total+m9_0total+m9_1total+m15_total
end_ind=beg+m8_total

if (identical(es_ds$Tumor_id[(beg+1):end_ind],as.character(m8_phe[[1]])))
{
  pheno_in=c(pheno_in,m8_phe$X2)
  pheno_level= c(pheno_level,m8_phe$X3)
}else{
  print("Error: Mouse8 pdl1")
}
####### check ###

length(pheno_in)
length(pheno_level)
## plug in ## 
colnames(es_ds@phenoData)
es_ds$Phenotype=pheno_in
es_ds$Degree=pheno_level

####  separate treatment into single column ## 
treat_in=rep(c("Control","Combo_0","Combo_1","TGFB","PDL1"),
    times=c(m12_total,m9_0total,m9_1total,m15_total,m8_total))
                                              
es_ds$Treatment=treat_in
saveRDS(es_ds,"/Users/sunq17/Desktop/new_bulk/data/NGS3953_ESET_do.Rdata")

# pheno = read.csv('data/NGS3953_pheno.csv',row.names = 1)
# count_m=read.csv('data/NGS3953_counts.csv',row.names = 1)  ### it is a  list ---Sun ##
# check1= readRDS("/Users/sunq17/Desktop/intern6-21/728a/data/NGS3947_ESET.Rdata")

### start  analyze ## 

es_ds=readRDS('data/NGS3953_ESET_do.Rdata')


##  combine left and right ear combo treatment to one variable Combo_all##
es_ds$Treatment=replace(es_ds$Treatment,es_ds$Treatment == "Combo_0","Combo_all")

es_ds$Treatment=replace(es_ds$Treatment,es_ds$Treatment == "Combo_1","Combo_all")

## typo in excel file ###
es_ds$Phenotype=replace(es_ds$Phenotype,es_ds$Phenotype == "infamed","inflamed")

###  exprs(es_ds) return matrix type ##
###  fData(es_ds) return data frame type ##

deg_obj = DGEList(exprs(es_ds),genes=fData(es_ds),group=es_ds$Treatment)
df=data.frame(Degree=es_ds$Degree,Pheno=es_ds$Phenotype)


# y.all$samples = cbind(y.all$samples, pheno)

###### filter out lowly expressed genes
keep.exprs = filterByExpr(deg_obj, group=es_ds$Treatment)
yf =  deg_obj[keep.exprs, keep.lib.sizes=FALSE]
yf$samples = cbind(yf$samples, df)

yf_pheno=yf[,yf$samples$Pheno!="undefined"]  ## delete undetermined phenotype ##
# count(yf$samples$Pheno=="undefined")

## tg_ds=yf[,yf$samples$group=="TGFB"]
yf=yf_pheno
yf = calcNormFactors(yf)

pheno_degree= paste(yf$samples$Pheno,yf$samples$Degree,sep = "_")

group_pheno=paste(yf$samples$group,yf$samples$Pheno,sep = "_")

treat_group=factor(yf$samples$group)

# covar_in=factor(group_pheno)
pheno_in=factor(yf$samples$Pheno)
# covar_in=factor(pheno_degree)
design=model.matrix(~0+covar_in)
colnames(design) = levels(covar_in)

#design = model.matrix(~0+treat_group)
colnames(design) = levels(treat_group)

## design of degree ##
design = model.matrix(~0+treat_group+covar_in) 
# colnames(design)
colnames(design) [1:4]= levels(treat_group)
colnames(design)[5:8] = levels(covar_in)[2:5]

# Run voom 
vmf=voomWithQualityWeights(yf, design, plot=T)
# Run Limma
fit=lmFit(vmf, design)

# Find differentially expressed genes
contrast.matrix = makeContrasts("Tgfb_exclu_infla" = TGFB_excluded - TGFB_inflamed,
                                "PDL1" = PDL1_inflamed - PDL1_excluded,
                                #"Combo_all" = Combo_all - Control,
                                levels=design)

## across all treatment ## 
contrast.matrix = makeContrasts("inf_st_ex" = inflamed_strong - excluded_median,
                                "inf_me_ex" = inflamed_median- excluded_median,
                                "inf_lo_ex" = inflamed_low - excluded_median,
                                levels=design)
contrast.matrix = makeContrasts("inf_st_med" = inflamed_strong - inflamed_median,
                                "inf_st_low" = inflamed_strong- inflamed_low,
                                "inf_med_low" = inflamed_median -inflamed_low,
                                levels=design)

contrast.matrix = makeContrasts("inf_ex" = inflamed - excluded, levels=design)

fit2 = contrasts.fit(fit,contrast.matrix)
fit2 = eBayes(fit2)
topTable(fit2)

dt = decideTests(fit2,lfc = 0,p.value = 0.05)
summary(decideTests(fit2))

## check genes in contrast $$ ##
p_home=fit2$p.value
is.data.frame(fit2$p.value)
typeof(fit2$p.value)
dim(fit2$p.value)
vennDiagram(dt)
##########################################PCA plot
ExprMat = cpm(yf,log=T)

annotations =data.frame(yf$samples$Pheno,yf$samples$group)

colnames(annotations)=c("Phenotype","treat_group")
rownames(annotations)=rownames(yf$samples)
# input ready
species = "Mus musculus"
logfold = 0
Contrastes = colnames(decideTests(fit2))[colSums(abs(decideTests(fit2))) > 0]
#Prefix = "Treatment_effect"
#Prefix="Treat_pheno_var"
Prefix = "inflamed_group"
Prefix = "inf_exc"
dir.create("plot_2home")
setwd("/Users/sunq17/Desktop/new_bulk/plot_2home/")


# all plots
BulkPlots(ExprMat,Prefix,Contrastes,annotations,covar_in,species,logfold)



##########################################PCA plot

var_genes = apply(ExprMat, 1, var)
var_genes = var_genes[order(var_genes,decreasing = T)]
nvar = 500       
topgenes = var_genes[1:nvar]
topdat = ExprMat[names(topgenes),]


### proofread typo ###
es_ds$Phenotype=replace(es_ds$Phenotype,es_ds$Phenotype == "infamed","inflamed")
####################################### PCA
res.pca = PCA(t(topdat), scale.unit = TRUE, ncp = 10, graph = F)

# pdf("PCA_NGS3947.pdf")
fviz_pca_ind(res.pca,
             axes = c(1,2),
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = as.character(es_ds$Phenotype), # color by groups
             palette = "viridis", 
             addEllipses = TRUE,
             label = "var",
             col.var = "black", repel = TRUE,
             legend.title = "",
             ellipse.level  = 0.90,
             pointshape = 21
) 


####  PCA per treatment 
library(dplyr)
ds_part=as_tibble(es_ds@phenoData@data)
colnames(ds_part)
head(ds_part)

ds_part %>% filter(Treatment=="Control") %>% dplyr::select(Sample)


####  study ### 
fit=glmQLFit(deg_obj,design,robust = TRUE)
qlf=glmQLFTest( fit,contrast=contrast.matrix )
summary(decideTests(qlf))


