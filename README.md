# NGS3953
BulkRNAseq-STAMP
### Description 
ex_new_bulk.R     
&nbsp;&nbsp; &nbsp;&nbsp; &nbsp;  _Provided that Excel file, copy cells in a range by hand for reading data._ <br/>
&nbsp;&nbsp; &nbsp;&nbsp; &nbsp;   _Darya_tsv_excel_eset data._<br/>

new_bulk.R        
&nbsp;&nbsp; &nbsp;&nbsp; &nbsp;   _Inspect input from Excel file, run differential expression analysis for treatment\_phenotype contrast._ <br/> 
contrast_a.R    
&nbsp;&nbsp; &nbsp;&nbsp; &nbsp;   _Differential expression analysis for phenotype\_degree contrast._ <br/>
pca_per_treat.R    
&nbsp;&nbsp; &nbsp;&nbsp; &nbsp;   _Subset of ExpressionSet object, determine optimal number of clusters, GO term_  <br/>
go_ana.R   
&nbsp;&nbsp; &nbsp;&nbsp; &nbsp;  _Heatmap annotation(treatment,pheno,degree) in order, row annotation._ <br/>
go_per_trt.R.   
&nbsp;&nbsp; &nbsp;&nbsp; &nbsp;   _Saving CSV file per treatment in a loop, heatmap of each treatment._ <br/>
correlation_cell.R.   
&nbsp;&nbsp; &nbsp;&nbsp; &nbsp;   _Given significant markers, correlation of differenct cell type._ <br/>
purify_pca.R.   
&nbsp;&nbsp; &nbsp;&nbsp; &nbsp;    _Residual plot , fit PC of Treg by PC of Cd8 in different treatment groups._<br/>
PurifyMarkers.R.   
&nbsp;&nbsp; &nbsp;&nbsp; &nbsp; _Purify sig markers from Atlas._ <br/>
NGS3953_purify_gene.R.  
&nbsp;&nbsp; &nbsp;&nbsp; &nbsp; _Avg expression value of cell type, pheno boxplot side by side among treatment._ <br/> 
