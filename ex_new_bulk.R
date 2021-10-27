####  In excel , copy cells in range by hand for reading data ###


es_ds=readRDS('data/NGS3953_ESET.Rdata')

# sample_inventory=data.frame(x=es_ds@phenoData@data$INVENTORY_SAMPLE_NAME )

# es_ds$Check=a


expand_inventory=sample_inventory %>% separate(x,"_", into =c("Tumor_id","T_name","Mouse","Treat","Time", "Tr_name")) 
table(expand_inventory$Mouse)

## verify mouse 12 _control  by # id ## 

m12_con_phe=read_excel('treatment bulkRNAseq manual classification for ximo.xlsx',
                       range = cell_cols("C:E"),col_names = FALSE,skip = 2)

m12_id=expand_inventory %>% filter( Mouse=="Mouse12R") %>% select(Tumor_id)

identical(m12_id$Tumor_id,m12_con_phe[[1]][3:nrow(m12_con_phe)])
# m12_id$Tumor_id
# m12_con_phe[[1]][3:nrow(m12_con_phe)]

##  mouse 9 R combo ## 
m9r_cob_phe=read_excel('treatment bulkRNAseq manual classification for ximo.xlsx',
                       range = cell_cols("G:I"),col_names = FALSE)

m9_id=expand_inventory %>% filter( Mouse=="Mouse9R") %>% select(Tumor_id)
m9_right=nrow(m9r_cob_phe)

identical(m9_id$Tumor_id[1:(m9_right-2)],m9r_cob_phe[[1]][3:m9_right])


##  mouse 9 L combo ## 
m9l_cob_phe=read_excel('treatment bulkRNAseq manual classification for ximo.xlsx',
                       range = cell_cols("K:M"),col_names = FALSE)

index_m9l=m9l_cob_phe[[1]][3:nrow(m9l_cob_phe)]