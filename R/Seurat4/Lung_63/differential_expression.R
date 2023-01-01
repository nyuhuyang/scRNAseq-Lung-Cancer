invisible(lapply(c("dplyr","magrittr","tidyr","openxlsx",#"Seurat","MAST","future",
                   "gplots"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
############### step = "resolutions" ############### 
csv_names = paste0("SCT_snn_res.",c(0.01, 0.1, 0.2, 0.5, 0.8, 0.9, 1, 2, 3,4,5))
csv_index = list.files("output/20220429",pattern = ".csv") %>% gsub("_.*","",.) %>% as.integer()
table(1:879 %in% csv_index)
csv_names = list.files("output/20220429",pattern = ".csv")
deg_list <- pbapply::pblapply(csv_names, function(csv){
        tmp <- read.csv(paste0("output/20220429/",csv),row.names = 1)
        tmp = tmp[tmp$p_val_adj < 0.05,]
        tmp$gene = rownames(tmp)
        tmp %<>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE)
        tmp$resolution = sub("^.*_SCT","SCT",csv) %>% 
                sub(".csv","",.) %>% 
                sub("SCT_snn_res","SCT-snn-res",.) %>%
                sub("_.*","",.) %>%
                sub("SCT-snn-res","SCT_snn_res",.)
        tmp
})

deg = bind_rows(deg_list)
deg %<>% filter(p_val_adj < 0.05)
deg_list = split(deg, f = deg$resolution)
write.xlsx(deg_list, file = paste0(path,"Lung_63_DEG.xlsx"),
           colNames = TRUE, borders = "surrounding")

############### step = "Adj_Dex_Cont" ############### 

csv_names = list.files("output/20220327/Cell_subtype",pattern = ".csv")
deg_list <- pbapply::pblapply(csv_names, function(csv){
        tmp <- read.csv(paste0("output/20220327/Cell_subtype/",csv),row.names = 1)
        #tmp = tmp[tmp$p_val_adj < 0.05,]
        tmp %<>% arrange(desc(avg_log2FC))
        tmp$gene = rownames(tmp)
        tmp
})

deg = bind_rows(deg_list)
deg1 = filter(deg, avg_log2FC > 0) %>% group_by(cluster) %>% top_n(50, avg_log2FC)
deg2 = filter(deg) %>% group_by(cluster) %>% top_n(100, abs(avg_log2FC))
deg_list = list(deg1,deg2)
names(deg_list) = c("postive","positve_negative")
write.xlsx(deg_list, file = "output/20220327/WC_30_T_Adj_Dex_vs_Count.xlsx",
           colNames = TRUE, borders = "surrounding")
############### step == "Cell_label" ############### 
opts = data.frame(ident = c(rep("Cell_label",66),
                            rep("Cell_type",31),
                            rep("Family",9),
                            rep("Superfamily",4)),
                  num = c(1:66,
                          1:31,
                          1:9,
                          1:4)
)
Cell_category = unique(opts$ident)
deg_list <- list()
for(i in seq_along(Cell_category)){
        csv_names = list.files("output/20220616/Cell_label",pattern = Cell_category[i],full.names = T)
        all_idx = which(opts$ident %in% Cell_category[i])
        idx <- gsub("output/20220616/Cell_label/","",csv_names) %>% gsub("-.*","",.) %>% as.integer()
        print(table(all_idx %in% idx))
        all_idx[!all_idx %in% idx]
        print(paste(Cell_category[i], "missing",all_idx[!(all_idx %in% idx)]))
        deg <- pbapply::pblapply(csv_names, function(csv){
                tmp <- read.csv(csv,row.names = 1)
                tmp$gene = rownames(tmp)
                tmp = tmp[order(tmp$avg_log2FC,decreasing = T),]
                if(tmp$cluster == TRUE) tmp$cluster = "T"
                return(tmp)
        }) %>% bind_rows
        deg = deg[deg$p_val_adj < 0.05,]
        deg_list[[i]] = deg
}
names(deg_list) =Cell_category

write.xlsx(deg_list, file = paste0(path,"Lung_63_DEG_Cell.category.xlsx"),
           colNames = TRUE, borders = "surrounding")

############### step == "celltype.3" ############### 
opts = data.frame(ident = c(rep("celltype.3",71),
                            rep("celltype.2",66),
                            rep("celltype.1",32),
                            rep("Family",10),
                            rep("Superfamily",4)),
                  num = c(1:71,
                          1:66,
                          1:32,
                          1:10,
                          1:4)
)
Cell_category = unique(opts$ident)
deg_list <- list()
for(i in seq_along(Cell_category)){
        csv_names = list.files("output/20221230/celltype.3",pattern = Cell_category[i],full.names = T)
        all_idx = which(opts$ident %in% Cell_category[i])
        idx <- gsub("output/20221230/celltype.3/","",csv_names) %>% gsub("-.*","",.) %>% as.integer()
        print(table(all_idx %in% idx))
        all_idx[!all_idx %in% idx]
        print(paste(Cell_category[i], "missing",all_idx[!(all_idx %in% idx)]))
        deg <- pbapply::pblapply(csv_names, function(csv){
                tmp <- read.csv(csv,row.names = 1)
                tmp$gene = rownames(tmp)
                tmp = tmp[order(tmp$avg_log2FC,decreasing = T),]
                if(tmp$cluster == TRUE) tmp$cluster = "T"
                return(tmp)
        }) %>% bind_rows
        deg = deg[deg$p_val_adj < 0.05,]
        deg_list[[i]] = deg
}
names(deg_list) =Cell_category

write.xlsx(deg_list, file = paste0(path,"Lung_63_DEG_Cell.category.xlsx"),
           colNames = TRUE, borders = "surrounding")

############### step == "IPF" ############### 
df_annotations <- readxl::read_excel("doc/20220816 IPF comparison.xlsx", sheet = "annotations")
df_comparision <- readxl::read_excel("doc/20220816 IPF comparison.xlsx", sheet = "comparison")
df_annotations %<>% tidyr::pivot_longer(everything()) %>%
        .[!duplicated(.[,"value"]),]
colnames(df_annotations)[1] = "type"
df_annotations$type %<>% factor(levels = rev(c("celltype.3","celltype.2","celltype.1","Family","Superfamily")))
df_annotations = df_annotations[order(df_annotations$type),]
opts = df_annotations[rep(seq_len(nrow(df_annotations)),
                          time = nrow(df_comparision)), ]
opts_comparision = df_comparision[rep(seq_len(nrow(df_comparision)),
                                      each = nrow(df_annotations)), ]
opts %<>% cbind(opts_comparision)

comparisions = unique(opts$sheetName)
deg_list <- list()
for(i in seq_along(comparisions)){
        csv_names = list.files("output/20220817/IPF",pattern = comparisions[i],full.names = T)
        all_idx = which(opts$sheetName %in% comparisions[i])
        idx <- gsub("output/20220817/IPF/","",csv_names) %>% gsub("-.*","",.) %>% as.integer()
        print(table(all_idx %in% idx))
        #print(paste(comparisions[i], "missing",all_idx[!(all_idx %in% idx)]))
        deg <- pbapply::pblapply(csv_names, function(csv){
                tmp <- read.csv(csv,row.names = 1)
                #tmp$gene = rownames(tmp)
                tmp = tmp[order(tmp$avg_log2FC,decreasing = T),]
                if(tmp$celltype == TRUE) tmp$celltype = "T"
                return(tmp)
        }) %>% bind_rows
        deg = deg[deg$p_val_adj < 0.05,]
        #deg = deg[deg$avg_log2FC > 0, ]
        deg_list[[i]] = deg
}
names(deg_list) =1:10

write.xlsx(deg_list, file = paste0(path,"Lung_63_DEGs_IPF.xlsx"),
           colNames = TRUE, borders = "surrounding")

# volcano plots
for(i in 1:2){
        csv_names = list.files("output/20220817/IPF",pattern = comparisions[i],full.names = T)
        all_idx = which(opts$sheetName %in% comparisions[i])
        idx <- gsub("output/20220817/IPF/","",csv_names) %>% gsub("-.*","",.) %>% as.integer()
        print(table(all_idx %in% idx))
        save.path <- paste0(path, comparisions[i],"/")
        if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
        options(ggrepel.max.overlaps = Inf)
        
        #print(paste(comparisions[i], "missing",all_idx[!(all_idx %in% idx)]))
        pbapply::pblapply(csv_names, function(csv){
                tmp <- read.csv(csv,row.names = 1)
                if(any(tmp$celltype == TRUE)) tmp$celltype = "T"
                cluster = stringr::str_split(tmp$cluster[1],patter = " vs.")[[1]]
                cluster = paste(rev(cluster),collapse = " <----    ----> ")
                subtitle = paste(c(cluster, tmp$celltype[1]),collapse = " \n in ")
                plot1 <- VolcanoPlots(tmp, cut_off ="p_val_adj", cut_off_value = 0.05,
                             cut_off_logFC = 0.25,cut_off_ptc = 10, top = 10,
                             sort.by = "p_val_adj",cols = c("#4575B4","#74ADD1","#E0F3F8","#FEE090","#F46D43"),
                             cols.inv = FALSE, alpha=0.9, pt.size=3, font.size=18,lab.size = 4,
                             subtitle = subtitle,inplab1 = "color text",
                             legend.show = TRUE,legend.size = 18, legend.position = "bottom",force = 2)
                jpeg(paste0(save.path,tmp$celltype[1],".jpg"), units="in", width=10, height=7,res=600)
                print(plot1)
                dev.off()
                
                plot2 <- VolcanoPlots(tmp, cut_off ="p_val_adj", cut_off_value = 0.05,
                                      cut_off_logFC = 0.25,cut_off_ptc = 10, top = 10,
                                      sort.by = "p_val_adj",cols = c("#4575B4","#74ADD1","#E0F3F8","#FEE090","#F46D43"),
                                      cols.inv = FALSE, alpha=0.9, pt.size=3, font.size=18,lab.size = 4,
                                      subtitle = subtitle,inplab1 = "No labels",
                                      legend.show = TRUE,legend.size = 18, legend.position = "bottom")
                jpeg(paste0(save.path,tmp$celltype[1],"_nolab.jpg"), units="in", width=10, height=7,res=600)
                print(plot2)
                dev.off()
                
        }) 
}

############### step == "pairwise" ############### 
xlx_file <- list.files("output/20221221/A.All samples combined",pattern = ".xlsx",full.names = TRUE)
degs_combine <- pbapply::pblapply(xlx_file,function(x) {
    tmp <- readxl::read_excel(x)
    tmp$Cell_category <- basename(x) %>% sub("2022-12-19-","",.) %>% sub("_combined.*","",.)
    colnames(tmp) <- c("celltype","gene","scores","avg_log2FC","p_val","p_val_adj","pts","pts_rest","Cell_category")
    tmp <- tmp[,c("p_val","avg_log2FC","pts","pts_rest","p_val_adj","scores","gene","celltype","Cell_category")]
    tmp <- tmp[order(tmp$avg_log2FC,decreasing = TRUE),]
    tmp
}) %>% bind_rows()

degs_combine %<>% filter(p_val_adj < 0.05)

openxlsx::write.xlsx(degs_combine, file =  paste0(path,"A.All samples combined.xlsx"),
                     colNames = TRUE,rowNames = FALSE,borders = "surrounding")
########
csv_files <- list.files("output/20221228",pattern = ".csv",recursive = TRUE,full.names = TRUE)

degs <- pbapply::pblapply(csv_files,function(x) {
    tmp <- read.csv(x,row.names = 1)
    tmp$celltype %<>% as.character()
    tmp$ident1 %<>% as.character()
    tmp$ident2 %<>% as.character()
    tmp$group <- sub("output/20221228/","",x) %>% sub("/.*","",.) %>% as.character()
    tmp <- tmp[order(tmp$avg_log2FC,decreasing = TRUE),]
    tmp
}) %>% bind_rows()
degs$ident2 %<>% gsub("FALSE","F",.)
degs$celltype %<>% gsub("TRUE","T",.)

degs %<>% filter(p_val_adj < 0.05)

deg_list <- split(degs,f = degs$group)
openxlsx::write.xlsx(deg_list, file =  paste0(path,"BtoK_pairwise_DEGs.xlsx"),
                     colNames = TRUE,rowNames = FALSE,borders = "surrounding")
