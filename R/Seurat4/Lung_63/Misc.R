library(Seurat)
library(dplyr)
library(tidyr)
library(magrittr)
library(readxl)
library(cowplot)
library(stringr)
library(openxlsx)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

meta.data = readRDS(file = "output/Lung_63_20220408_meta.data_v4.rds")
meta.data1 = readRDS(file = "../scRNAseq-Lung/output/20211222/meta.data_SCINA_Lung30_Azimuth_Cell_Types_2021.rds")

meta.data = meta.data[meta.data$orig.ident %in% meta.data1$orig.ident,]
table(rownames(meta.data) %in% rownames(meta.data1))
meta.data$barcode = rownames(meta.data)
meta.data1$barcode = rownames(meta.data1)

meta.data2 = inner_join(meta.data1[,c("Cell_subtype","barcode")],
                        meta.data[,c("Cell_label","barcode","Doublets")],by = "barcode")
#meta.data3 = meta.data2[,c("Doublets.x","Doublets.y")]
#colnames(meta.data3) = c("Cell_subtype","Cell_label")
#meta.data2 = rbind(meta.data2[,c("Cell_subtype","Cell_label")],
#                   meta.data3[,c("Cell_subtype","Cell_label")])
#lvl = levels(meta.data2$Cell_subtype)
#meta.data2$Cell_subtype %<>% as.character()
#meta.data2$Cell_subtype[is.na(meta.data2$Cell_subtype)] = "Unlabel"
#meta.data2$Cell_subtype %<>% factor(levels = c(lvl,"Unlabel"))
df <- table(meta.data2$Cell_subtype, meta.data2$Cell_label) %>% as.data.frame.matrix()
df_doublet =  table(meta.data2$Cell_subtype, meta.data2$Doublets) %>% as.data.frame.matrix()
df %<>% cbind(df_doublet)
df_colSums = colSums(df)
df_output = rbind(df, "Total cell number"=c(df_colSums))
df_output %<>% rbind("---------------------" = NA)

df1 <- table(meta.data2$Cell_subtype, meta.data2$Cell_label) %>% prop.table(margin = 2)%>% 
    as.data.frame.matrix()
df1_doublet =  table(meta.data2$Cell_subtype, meta.data2$Doublets) %>% prop.table(margin = 2)%>% 
    as.data.frame.matrix()
df1 %<>% cbind(df1_doublet)
df1[is.na(df1)] = 0
df1 = round(df1*100, 2)

df_output %<>% rbind(df1)  
openxlsx::write.xlsx(df_output, file =  paste0(path,"20220722_annotation_comparsion.xlsx"),
                     colNames = TRUE,row.names = T,borders = "surrounding")
#======================================

meta.data2 = inner_join(meta.data1[,c("orig.ident","Cell_subtype","barcode")],
                        meta.data[,c("Cell_label","barcode","Doublets")],by = "barcode")
meta.data_list <- split(meta.data2,f = meta.data2$orig.ident)
res_list <- lapply(meta.data_list, function(data){
    df <- table(data$Cell_subtype, data$Cell_label) %>% as.data.frame.matrix()
    df_doublet =  table(data$Cell_subtype, data$Doublets) %>% as.data.frame.matrix()
    df %<>% cbind(df_doublet)
    df_colSums = colSums(df)
    df_output = rbind(df, "Total cell number"=c(df_colSums))
    df_output %<>% rbind("---------------------" = NA)
    
    df1 <- table(data$Cell_subtype, data$Cell_label) %>% prop.table(margin = 2)%>% 
        as.data.frame.matrix()
    df1_doublet =  table(data$Cell_subtype, data$Doublets) %>% prop.table(margin = 2)%>% 
        as.data.frame.matrix()
    df1 %<>% cbind(df1_doublet)
    df1[is.na(df1)] = 0
    df1 = round(df1*100, 2)
    
    df_output %<>% rbind(df1)  
})
res_list$Total = df_output
res_list = res_list[c("Total",names(meta.data_list))]
openxlsx::write.xlsx(res_list, file =  paste0(path,"20220722_annotation_comparsion_by_sample.xlsx"),
                     colNames = TRUE,row.names = T,borders = "surrounding")
