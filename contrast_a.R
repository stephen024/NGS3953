

setwd("/Users/sunq17/Desktop/new_bulk/")

source("/Users/sunq17/Desktop/intern6-21/728a/BulkRNAseq_AnalysisFunctions.R")

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

trt_phe_de=paste(yf$samples$group,factor(pheno_degree),sep = "_")

treat_group=factor(yf$samples$group)

pheno_degree=as.factor(pheno_degree)

degree_in=factor(yf$samples$Degree)
pheno_in=factor(yf$samples$Pheno)

design=model.matrix(~0+pheno_in)
colnames(design)= levels(pheno_in)
#colnames(design)[5:8] = levels(pheno_degree)[2:5]
design=cbind(design,X1= treat_group== "TGFB" & degree_in=="median" )
design=cbind(design,X2= treat_group== "TGFB" & degree_in=="low" )


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
fit=lmFit(vmf,design)

# Find differentially expressed genes
contrast.matrix = makeContrasts("inf_st_ex" = inflamed_strong - excluded_median,
                                "inf_st_med" = inflamed_strong - inflamed_median,
                                
                                levels=design)


contrast.matrix = makeContrasts("tg_" = inflamed - excluded, levels=design)

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


