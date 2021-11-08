###  ###


# source("/Users/sunq17/Desktop/intern6-21/728a/BulkRNAseq_AnalysisFunctions.R")
library(boot)

### load csv file to have a vector of genes and its cluster(cell type)###
marker_df=read.csv("sig_purify_genes.csv")
#colnames(marker_df)
#length(marker_df$gene)

es_ds=readRDS('/Users/sunq17/Desktop/new_bulk/data/NGS3953_ESET_do.Rdata')

##  combine left and right ear combo treatment to one variable Combo_all##
es_ds$Treatment=replace(es_ds$Treatment,es_ds$Treatment == "Combo_0","Combo")
es_ds$Treatment=replace(es_ds$Treatment,es_ds$Treatment == "Combo_1","Combo")

## typo in excel file ###
es_ds$Phenotype=replace(es_ds$Phenotype,es_ds$Phenotype == "infamed","inflamed")

### 
deg_obj = DGEList(exprs(es_ds),genes=fData(es_ds),group=es_ds$Treatment)
df=data.frame(Degree=es_ds$Degree,Phenotype=es_ds$Phenotype)

###  keep.exprs = filterByExpr(deg_obj, group=es_ds$Treatment) ### 
### Because you have gene list , no need to use this  ### 

gene_in= rownames(deg_obj) %in% marker_df$gene 


#sum(gene_in)
#setdiff(marker_df$gene,rownames(deg_obj))

yf =  deg_obj[gene_in, ,keep.lib.sizes=FALSE]

yf$samples = cbind(yf$samples, df)

yf=yf[,yf$samples$Pheno!="undefined",keep.lib.sizes=FALSE] 

yf = calcNormFactors(yf)

ExprMat = cpm(yf,log=T) 

n_cluster=unique(marker_df$cluster)
# pcSig_list=vector(mode = "list",length=length(n_cluster)) 

for (i in 1:length(n_cluster)) {
  #i=2
  gene_cluster=marker_df %>% filter(cluster==n_cluster[i]) %>% dplyr::select(gene)
  
  clus_ind=rownames(ExprMat) %in% gene_cluster$gene
  #sum(clus_ind)
  
  topdat = ExprMat[clus_ind,]
  
  pcSig = gsScore(topdat)
  
  name_col=paste("PC",n_cluster[i],sep = "_")
  #yf$samples$name_col = pcSig
  
  yf$samples=cbind(yf$samples,pcSig)
  
  colnames(yf$samples)[colnames(yf$samples)=='pcSig']=name_col
}



### cell type starting from 6th column ##
pc_matrix=yf$samples[,6:30]
library(heatmaply)
heatmaply::heatmaply_cor(x=cor(pc_matrix,method = "spearman"),fontsize_row = 6,fontsize_col = 6)


library(ggcorrplot)

corr=round(cor(pc_matrix,method = "spearman"),3)
p_mat=cor_pmat(pc_matrix,method = "spearman")
corr.plot <- ggcorrplot(
  corr=corr, hc.order = TRUE, type = "lower", outline.col = "white",
  p.mat = p_mat,tl.cex = 7
)
#dev.off()
corr.plot



yf$samples$Complete = paste(yf$samples$group,yf$samples$Pheno,sep = "_")

## delete phenotype desert ###
## select two cell type you are interested in ##

no_ex_yf=yf$samples %>% filter(Phenotype %in% c("excluded","inflamed"))%>%dplyr::select(c('PC_Cd8','PC_Tregs','Complete')) 

no_ex_yf$Complete=factor(no_ex_yf$Complete)

table(no_ex_yf$Complete)

prop_fun<-function(data,idx,cor_type="spearman"){
  df=data[idx,]
  c(cor(df[,1],df[,2],method =cor_type))
}

## begin bootstrap for correlation coefficient  ##
im_pheno=unique(no_ex_yf$Complete)

bci_list=vector(mode = "list",length=length(im_pheno)) 

boot_df = data.frame(cluster=character(),boot_out=numeric())

for (k in 1:length(im_pheno)) {
  im_name=im_pheno[k]
  
  sub_yf=no_ex_yf %>% filter(Complete==im_name)
  myboot<-boot::boot(sub_yf,prop_fun,R=500)
 
  bci=boot::boot.ci(myboot, index=1, type='perc')$percent[4:5]
  
  bci_list[[k]]=bci
  bo_df=data.frame(cluster=rep(im_name,500),boot_out=myboot$t)
  
  boot_df=rbind(boot_df,bo_df)
}

names(boot_df)

##  boxplot of bootstrapped coorealtion coefficient ##
ggboxplot(boot_df, x = "cluster", y = "boot_out")+
  ylab("Coef")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))



#### below is to see proportion of Treg over Cd8  ####
  
prop_df = data.frame(cluster=character(),prop=numeric())

for (k in 1:length(im_pheno)) {
  im_name=im_pheno[k]
  
  sub_yf=no_ex_yf %>% filter(Complete==im_name)
   
  treg=sub_yf$PC_Tregs/sub_yf$PC_Cd8
 
   pp_df=data.frame(cluster=rep(im_name,nrow(sub_yf)),prop=treg)
  
  prop_df=rbind(prop_df,pp_df)
}

ggboxplot(prop_df, x = "cluster", y = "prop")+
  ylab("Treg over Cd8")+ylim(-2,2)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

### linear fit btw Tregs and Cd8 ####
ggscatter(no_ex_yf, x = "PC_Cd8", y = "PC_Tregs",
          add = "reg.line",                         
          conf.int = TRUE,                          
          color = "Complete", palette = "jco",           # Color by Treatment_pheno 
          shape = "Complete"                             # Change point shape by groups 
)+ylim(-2,2)
 


##  fit Treg and Cd8 in control_inflamed dataset ##
fit_ds=no_ex_yf %>% filter(Complete=="Control_inflamed")  
fit<-lm(PC_Tregs~PC_Cd8,data=fit_ds)    
###  R square is very low , bad fit only btw treg and cd8##
#summary(fit)
#par(mfrow = c(2, 2)) 
#plot(fit) 

pred_df = data.frame(matrix(ncol=5,nrow=0))
col_name=c("PC_Cd8","PC_Tregs","Complete","predicted","Diff")
colnames(pred_df)=col_name

for(j in 1:length(im_pheno)){

  im_name=im_pheno[j]
  
  con_ex=no_ex_yf %>% filter(Complete==im_name)
  con_ex$predicted=predict(object=fit,newdata=con_ex %>% dplyr::select(c(1,2)))
  
  con_ex=con_ex%>% mutate(Diff=PC_Tregs-predicted)
  
  pred_df=rbind(pred_df,con_ex)
}

diff_box<-ggplot(data=pred_df,aes(x=Complete,y=Diff,color=Complete))+
  geom_boxplot(fill = NA)+
  labs(title = paste("Deviation from control_inflamed fitted line"))+
  theme(plot.title = element_text(hjust = 0.5),legend.key = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

diff_box

### linear model and deviation plot ###
##  each treatment_pheno has a plot ##

for(j in 1:length(im_pheno)){
  
  im_name=im_pheno[j]
  
  con_ex=no_ex_yf %>% filter(Complete==im_name)
  con_ex$predicted=predict(object=fit,newdata=con_ex%>% select(c(1,2)))
  
  con_ex=con_ex%>% mutate(Diff=PC_Tregs-predicted)
  
  png(paste(im_name,"box.png",sep="_"), units="px", width=597, height=441)
  print(
    ggplot(fit_ds, aes(x = PC_Cd8, y = PC_Tregs)) +
      geom_point(colour="tan1")+
      geom_smooth(method = "lm", se = FALSE, color = "#C0C0C0")+
      geom_segment(data=con_ex,mapping = aes(x=PC_Cd8,y=predicted,xend = PC_Cd8, yend = PC_Tregs), 
                   alpha = .3) +
      geom_point(data=con_ex%>%select(1,2),colour="dodger blue")+
      labs(title = paste(im_name,"Deviation"),subtitle ="Orange:control_inflamed")+
      theme(plot.title = element_text( size = 12, face = "bold"),
            plot.subtitle = element_text(colour = "Orange",face = "italic"),
            
            panel.background = element_rect(fill = "white", size = 1.5,colour = "grey", 
                                            linetype = "solid"),
            panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                            colour = "gainsboro"), 
            panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                            colour = "gainsboro"))
  )
  dev.off()
}


