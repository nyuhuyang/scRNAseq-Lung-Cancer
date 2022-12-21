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


meta.data %<>% cbind(object[["umap"]]@cell.embeddings)
meta.data %<>% cbind(meta.data_exp)
table(rownames(meta.data) == colnames(object))
#======== rename ident =================
#meta.data = readRDS(file = "output/Lung_63_20220408_meta.data_v3.rds")
SCT_snn_res <- grep("SCT_snn_res",colnames(meta.data),value =T)
for(res in SCT_snn_res){
    meta.data[,res] %<>% as.character() %>% as.integer() %>% as.factor()
}

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

#===================================================================================
object = readRDS(file = "data/Lung_SCT_63_20220606.rds")
meta.data = readRDS(file = "output/Lung_63_20220408_meta.data_v4.rds")
table(rownames(object@meta.data) == rownames(meta.data))
DefaultAssay(object) = "SCT"
#meta.data = object@meta.data
#umap = object[["umap"]]@cell.embeddings
meta.data_exp = FetchData(object, vars = unique(c("ACKR1","ACP5","ACTA2","AGER","APOC1","BCHE","C1QB",
                                           "CA4","CCL18","CCL3","CCL4","CCNO","CCL3",
                                           "CCL4","CCR7","CCR3","CCR4","CD1A","CD1C","CD3D","CD3E",
                                           "CD4","CD83","CD8A","CD8B","CDH5","CLDN5",
                                           "COX4I2","CXCL10","CXCL8","CXCL9","CXCR6",
                                           "DES","EDNRB","EREG","FABP4","FCER1A","GJA5",
                                           "GNLY","HES6","IFNG","IGFBP3","IL1B",
                                           "IL1RL1","KCNK3","KLRD1",
                                           "KRT15","KRT5","LTB","LTF","MAL","MARCO",
                                           "MSR1","MUC5AC","MUC5B","NKG7",
                                           "NOTCH3","PRB3","PRR4","S100A8","S100A9",
                                           "SAA1","SAA2","SAA4","SCGB1A1","SCGB3A2",
                                           "SELE","SERPINB3","SERPINB4","SFTPB","SFTPC",
                                           "SPP1","TAGLN")))
object[["umap"]] = NULL
file.name = "output/20220420/3000/umap_cs100_dist.0.1_spread.1.5.rds"
umap  = readRDS(file.name)[[1]]
umap@key = "UMAP_"
colnames(umap@cell.embeddings) = c("UMAP_1","UMAP_2")
object[["umap"]] <- umap
meta.data.v5 <- meta.data
meta.data %<>% cbind(object[["umap"]]@cell.embeddings)
meta.data %<>% cbind(meta.data_exp)
table(rownames(meta.data) == colnames(object))

df_annotation <- readxl::read_excel("doc/Annotations/63-samples annotation optimization - part 1 and 2-YH v2 RS.xlsx",
                                    sheet = "Sheet1")

colnames(meta.data) %<>% gsub("UMAP_","UMAP",.)
df_annotation$condition[is.na(df_annotation$condition)] = "TRUE"
df_annotation$`total cells` = NULL
df_annotation$`list.of.sample.IDs` = NULL
Cell_label = unique(c(df_annotation$Cell_label,df_annotation$label))
table(Cell_label %in% meta.data$Cell_label)
Cell_label[!(Cell_label %in% meta.data$Cell_label)]
total_cells = c()
list_of_sample_IDs = c()
meta.data$Cell_label %<>% as.character()
for(m in 1:nrow(df_annotation)){
    cl = pull(df_annotation[m,"Cell_label"])
    step = pull(df_annotation[m,"steps"])
    change_to = pull(df_annotation[m,"label"])
    
    select_id = meta.data %>% dplyr::filter(Cell_label %in% cl) %>% 
        dplyr::filter(eval(parse(text = df_annotation$condition[m])))
    if(nrow(select_id) == 0) {
        total_cells[m] = nrow(select_id)
        list_of_sample_IDs[m] = "0"
        next
    }
    meta.data[rownames(select_id),"Cell_label"] = change_to
    total_cells[m] = nrow(select_id)
    list_of_sample_IDs[m] = gsub("-.*","",rownames(select_id)) %>% 
                            table %>% 
                            as.data.frame.table %>% 
                            apply(1,function(x) paste(x,collapse= " x ")) %>% 
                            paste(collapse= " cells, ") %>%
                            paste("cells")
    print(paste ("Cell_label at",cl,df_annotation$condition[m],"------->",change_to,
                 "total cells = ",total_cells[m],
                 " list of sample IDs = ", list_of_sample_IDs[m]))
}
# save annotation notes
meta.data.v5$celltype.3 = meta.data$Cell_label
meta.data.v5$celltype.3 %<>% gsub("En-SM","En-sm",.)
meta.data.v5$celltype.3 %<>% gsub("En-ca-F4","En-ca,Fb-a",.)
saveRDS(meta.data.v5,file = "output/Lung_63_20220408_meta.data_v5.rds")

df_annotation %<>% cbind(data.frame("total cells" = total_cells))
df_annotation %<>% cbind(data.frame("list of sample IDs" = list_of_sample_IDs))
openxlsx::write.xlsx(df_annotation, file =  paste0(path,"63-samples annotation optimization - part 1 and 2-YH v3.xlsx"),
                     colNames = TRUE,rowNames = F,borders = "surrounding")

df_samples <- readxl::read_excel("doc/20220406-samples metadata RS.xlsx", sheet = "RS in vivo metadata")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
nrow(df_samples)
df_samples = df_samples[order(df_samples$idx),]
meta.data$orig.ident %<>% factor(levels = df_samples$sample)
df <- table(meta.data$Cell_label, meta.data$orig.ident) %>% as.data.frame.matrix()
df_colSums = colSums(df)
df_output = df_samples[,c("sample","study group")] %>% tibble::column_to_rownames("sample") %>% t
df_output %<>% rbind(df)
df_output %<>% rbind( "Total cell number"=c(df_colSums))
df_output %<>% rbind("---------------------" = NA)

df1 <- table(meta.data$Cell_label, meta.data$orig.ident) %>% prop.table(margin = 2)%>% 
    as.data.frame.matrix()
df1[is.na(df1)] = 0
df1 = round(df1*100, 2)

df_output %<>% rbind(df1)  
openxlsx::write.xlsx(df_output, file =  paste0(path,"20220805_annotation_part 1 and 2.xlsx"),
                     colNames = TRUE,row.names = T,borders = "surrounding")

table(Cell_label.v4 %in% df_color[,1])
Cell_label.v4 = sort(as.character(unique(meta.data.v4$Cell_label)))
Cell_label.v5 = sort(as.character(unique(meta.data$Cell_label)))

table(Cell_label.v4 %in% Cell_label.v5)
table(Cell_label.v5 %in% Cell_label.v4)
Cell_label.v5[!(Cell_label.v5 %in% Cell_label.v4)]
df_annotation1 <- readxl::read_excel("doc/Annotations/annotations.xlsx",
                                    sheet = "annotations")
meta.data.v5 = readRDS(file = "output/Lung_63_20220408_meta.data_v5.rds")
df_annotation1 = df_annotation1[!(df_annotation1$celltype.3 %in% "Db"),]
celltype.3 = sort(unique(df_annotation1$celltype.3))
celltype.3 = celltype.3[-which(celltype.3 %in% "Db")]

Cell_types <- c("celltype.3.colors","celltype.2","celltype.1","Family","Superfamily","Pathology")
for(Cell_type in Cell_types){
    meta.data.v5[,Cell_type] = plyr::mapvalues(meta.data.v5$celltype.3,
                                            from = pull(df_annotation1[,"celltype.3"]),
                                            to = pull(df_annotation1[,Cell_type]))
    meta.data.v5[,Cell_type] %<>% factor(levels = unique(pull(df_annotation1[,Cell_type])))
}
meta.data.v5$celltype.3  %<>% droplevels()
meta.data.v5$celltype.3 %<>% factor(levels = df_annotation1$celltype.3)
saveRDS(meta.data.v5, file = "output/Lung_63_20220408_meta.data_v5.rds")

meta.data.v5 = readRDS( "output/Lung_63_20220408_meta.data_v5.rds")
annotations <- readxl::read_excel("doc/Annotations/annotations.xlsx")
annotations = as.data.frame(annotations)
colnames(annotations) %<>% tolower()
nrow(annotations)

df_annotations = data.frame("celltype" = c(annotations$celltype.1,
                                           annotations$celltype.2,
                                           annotations$celltype.3),
                            "celltype.colors" = c(annotations$celltype.1.colors,
                                                   annotations$celltype.2.colors,
                                                   annotations$celltype.3.colors))
df_annotations1 = df_annotations[!duplicated(df_annotations$celltype),]
df_annotations1 = df_annotations1[order(df_annotations1$celltype),]
df_annotations2 = df_annotations[!duplicated(df_annotations$celltype.colors),]
df_annotations2 = df_annotations2[order(df_annotations2$celltype),]

identical(df_annotations1,df_annotations2)
table(duplicated(df_annotations1$celltype))

df_samples <- readxl::read_excel("doc/20220406-samples metadata RS.xlsx", sheet = "invivo")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
nrow(df_samples)



table(df_samples$sample %in% meta.data.v5$orig.ident)
for(i in 1:length(df_samples$sample)){
    cells <- meta.data.v5$orig.ident %in% df_samples$sample[i]
    print(df_samples$sample[i])
    print(table(cells))
    meta.data.v5[cells,"type of sample"] = df_samples$`type of sample`[i]
    meta.data.v5[cells,"category"] = df_samples$category[i]
    meta.data.v5[cells,"group3"] = df_samples$group3[i]
    meta.data.v5[cells,"group2"] = df_samples$group2[i]
    meta.data.v5[cells,"group1"] = df_samples$group1[i]
    meta.data.v5[cells,"group0"] = df_samples$group0[i]
    meta.data.v5[cells,"lung disease-1"] = df_samples$`lung disease-1`[i]
    meta.data.v5[cells,"lung disease-2"] = df_samples$`lung disease-2`[i]
    meta.data.v5[cells,"lung disease-3"] = df_samples$`lung disease-3`[i]
    meta.data.v5[cells,"sex"] = as.character(df_samples$sex[i])
}
meta.data.v5$orig.ident %<>% factor(levels = df_samples$sample)
meta.data.v5$celltype.3.colors = plyr::mapvalues(meta.data.v5$celltype.3,
                                                 from = df_annotations1$celltype,
                                                 to = df_annotations1$celltype.colors)
meta.data.v5$celltype.2.colors = plyr::mapvalues(meta.data.v5$celltype.2,
                                                 from = df_annotations1$celltype,
                                                 to = df_annotations1$celltype.colors)
meta.data.v5$celltype.1.colors = plyr::mapvalues(meta.data.v5$celltype.1,
                                                 from = df_annotations1$celltype,
                                                 to = df_annotations1$celltype.colors)
anyNA(meta.data.v5$celltype.3.colors)
anyNA(meta.data.v5$celltype.2.colors)
anyNA(meta.data.v5$celltype.1.colors)
category <- c("P-norm",
              "D-norm",
              "T-norm",
              "D-COPD",
              "D-IPF",
              "L-norm",
              "L-COPD",
              "L-COPD-Dex",
              "L-IPF",
              "L-IPF-norm",
              "L-Ad",
              "L-Sq",
              "L-Ad-Sq")

meta.data.v5$category %<>% factor(levels = category)
meta.data.v5$category.colors <- plyr::mapvalues(meta.data.v5$category,
                                                from = category,
                                                to = c("#00B0F0",
                                                       "#92D050",
                                                       "#00B050",
                                                       "#FCC4F5",
                                                       "#FF9933",
                                                       "#4472C4",
                                                       "#AFABAB",
                                                       "#FF6699",
                                                       "#DAB48E",
                                                       "#FFD966",
                                                       "#B4C7E7",
                                                       "#AE78D6",
                                                       "#92C0F2"))
meta.data.v5$condition %<>% gsub("D-norm","Normal",.)
saveRDS(meta.data.v5, file = "output/Lung_63_20220408_meta.data_v5.rds")

