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
# load data
object = readRDS(file = "data/Lung_SCT_63_20220606.rds")
meta.data = readRDS(file = "output/Lung_63_20220408_meta.data_v2.rds")
DefaultAssay(object) = "SCT"
#meta.data = object@meta.data
#umap = object[["umap"]]@cell.embeddings
meta.data_exp = FetchData(object, vars = c("GNLY","SFTPC","SFTPB","SAA1","SAA2","SCGB1A1","SCGB3A2","PRF1",
                                           "CD3E","CD3D","CD8A","CD8B","MARCO"))

#saveRDS(umap, "output/20210901/umap_dist.0.3_spread.1.rds")
#saveRDS(meta.data_exp, "output/20210901/meta.data_exp.rds")
object[["umap"]] = NULL
file.name = "output/20220420/3000/umap_cs100_dist.0.1_spread.1.5.rds"
umap  = readRDS(file.name)[[1]]
umap@key = "UMAP_"
colnames(umap@cell.embeddings) = c("UMAP_1","UMAP_2")
object[["umap"]] <- umap


SCT_snn_res <- grep("SCT_snn_res",colnames(meta.data),value =T)
for(res in SCT_snn_res){
    meta.data[,res] %<>% as.character() %>% as.integer() %>% as.factor()
}

meta.data %<>% cbind(object[["umap"]]@cell.embeddings)
meta.data %<>% cbind(meta.data_exp)
table(rownames(meta.data) == colnames(object))
#======== rename ident =================
#meta.data = readRDS(file = "output/Lung_63_20220408_meta.data_v3.rds")

df_annotation <- readxl::read_excel("doc/Annotations/20220602_Annotations-3000-100-01-1.5-RS-06-03-22 HY.xlsx",
                                    sheet = "Sheet1")
res = "SCT_snn_res.5"
df_annotation1 = df_annotation[!is.na(df_annotation[,"cluster"]),]

cluster_check = c()
for(i in 1:nrow(df_annotation1)){
    cl = df_annotation1$SCT_snn_res.5[i]
    clusters <- df_annotation1[i,"cluster"] %>% pull %>% strsplit(", ") %>% .[[1]] %>% as.integer
    cluster_check[i] <- cl %in% clusters
}
table(cluster_check)
df_annotation1[which(!cluster_check),]
df_annotation2 = df_annotation1[!is.na(df_annotation1[,"SCT_snn_res.5"]),]
df_annotation2 = df_annotation2[!duplicated(df_annotation2[,"SCT_snn_res.5"]),]
meta.data[,"Cell_label"] = plyr::mapvalues(meta.data[,res],
                                             from = pull(df_annotation2[,res]),
                                             to = pull(df_annotation2[,"Cell_label"])
                                             ) %>% as.character()

colnames(meta.data) %<>% gsub("UMAP_","UMAP",.)
keep = !is.na(df_annotation$modify_condition)
df_annotation4 = df_annotation[keep,c("SCT_snn_res.5","modify_condition","modify_label")]
colnames(df_annotation4) %<>% gsub("modify_","",.)

for(m in 1:nrow(df_annotation4)){
    cl = pull(df_annotation4[m,res])
    change_to = pull(df_annotation4[m,"label"])
    
    select_id = meta.data %>% dplyr::filter(!!as.name(res) %in% cl) %>% 
                          dplyr::filter(eval(parse(text = df_annotation4$condition[m])))
    meta.data[rownames(select_id),"Cell_label"] = change_to
    print(paste (res,"at",cl,df_annotation4$condition[m],"------->",change_to))
}

table((meta.data$Cell_label %in% "MC") & (meta.data$UMAP1 < -5 | meta.data$UMAP1  > 3 |  meta.data$UMAP2 > 0 | meta.data$UMAP2  < -5)
)
table(meta.data$Cell_label %in% "TASC")
#  164 170 ======
#reName = meta.data$Cell_label %in% c("164","170")
#meta.data[reName,"Cell_label"] = meta.data[reName,"Cell_label"]

df_annotation <- readxl::read_excel("doc/Annotations/20220602_Annotations-3000-100-01-1.5-RS-06-03-22 HY.xlsx",
                                    sheet = "63-sample annotations")

df_annotation = df_annotation %>% dplyr::filter(!is.na(Cell_label)) %>% dplyr::filter(!duplicated(Cell_label))

Cell_types <- c("Cell_type","Family","Superfamily","Pathology")
for(Cell_type in Cell_types){
    meta.data[,Cell_type] = plyr::mapvalues(meta.data$Cell_label,
                                            from = pull(df_annotation[,"Cell_label"]),
                                            to = pull(df_annotation[,Cell_type]))
    meta.data[,Cell_type] %<>% factor(levels = unique(pull(df_annotation[,Cell_type])))
}
meta.data$Cell_label %<>% factor(levels = df_annotation$Cell_label)
saveRDS(meta.data[,-c(46:60)], file = "output/Lung_63_20220408_meta.data_v3.rds")


meta.data1 = readRDS(file = "output/Lung_63_20220408_meta.data_v4.rds")

Cell_types <- c("Cell_label","Cell_type","Family","Superfamily","Pathology")
for(Cell_type in Cell_types){
    meta.data1[,Cell_type] = meta.data[,Cell_type]
}
    
saveRDS(meta.data1, file = "output/Lung_63_20220408_meta.data_v4.rds")


# add color
df_color = t(data.frame(
    c("AT1","#C946D4"),#yes
    c("AT2","#A794D7"),#yes
    c("B","#B4C7E7"),#yes
    c("BC","#70AD47"),#yes
    c("C-s","#FFE699"),#yes
    c("C1","#FFC000"),#yes
    c("CD8-T1","#5B9BD5"),
    c("c-DC","#0070C0"),
    c("Cr","#B8C6DA"),#yes
    c("En-a","#00B0F0"),#yes
    c("En-c1","#A7E5BC"),#yes
    c("En-ca","#BDD7EE"),#yes
    c("En-l","#7CAFDD"),#yes
    c("En-SM","#FFF2CC"),
    c("En-v","#C7BCE7"),
    c("Fb1","#FFA919"),#yes
    c("Fb2","#C00000"),#yes
    c("Fb3","#FA7FA9"),#yes
    c("Fb4","#77D900"),#yes
    c("G-Muc","#00B050"),#yes
    c("G-Ser","#ADCDEA"),#yes
    c("Gli","#8FAADC"),#yes
    c("H","#EBE621"),#yes
    c("IC","#FBE5D6"),#yes
    c("Ion","#0F23FF"),#yes
    c("M1","#FF6600"),#yes
    c("M1-2","#D29F26"),#yes
    c("M2","#BFBFBF"),#yes
    c("MC","#2FE6F9"),
    c("ME","#E169CD"),#yes
    c("Mon","#FBC8EF"),#yes
    c("NE","#FF0000"),#yes
    c("Neu","#FF7C88"),#yes
    c("NK","#FF8B6A"),#yes
    c("p-C","#EB6C11"),#yes
    c("PC","#ED7D31"),#yes
    c("pDC","#E31A1C"),
    c("Pr","#747474"),#yes
    c("S-Muc","#B38600"),
    c("S1","#FFBAD1"),#yes
    c("SM1","#C8D8E0"),#yes
    c("SM2","#0087B6"),#yes
    c("SM3","#FFFF00"),#yes
    c("T-ifn","#FD5A00"),#yes
    c("TASC","#FF3990"),#yes
    c("Tcn","#9EF971"),#yes
    c("T-NK","#FF42A4"),#yes
    c("Trm","#f1c232"),#fc8e66
    c("T-un","#b1bcc5"),
    c("Un","#DEEBF7")))#yes
df_color %<>% as.data.frame
rownames(df_color) = NULL
colnames(df_color) = c("Cell_label","Cell_label.colors")
table(duplicated(df_color$Cell_label.colors))

meta.data_u = meta.data[!duplicated(meta.data$Cell_label),]
meta.data1 = meta.data_u[!(meta.data_u$Cell_label %in% df_color$Cell_label),];print(dim(meta.data1))
meta.data2 = meta.data_u[(meta.data_u$Cell_label %in% df_color$Cell_label),];print(dim(meta.data2))

meta.data$Cell_label.colors =  plyr::mapvalues(meta.data$Cell_label,
                                                 from = df_color$Cell_label,
                                                 to = df_color$Cell_label.colors)

df_annotation <- readxl::read_excel("doc/Annotations/20210917_20UMAP res0.8 annotations.xlsx",
                                    sheet = "Sheet1")
Cell_labels <- c("Cell_label","Cell_type","UMAP_land","Family","Superfamily")

df_annotation = df_annotation[order(df_annotation$Cell_label),Cell_labels]
df_annotation = df_annotation[!duplicated(df_annotation$Cell_label),]
for(Cell_label in Cell_labels[2:5]){
    meta.data[,Cell_label] = plyr::mapvalues(meta.data$Cell_label,
                                            from = pull(df_annotation[,"Cell_label"]),
                                            to = pull(df_annotation[,Cell_label]))
}

