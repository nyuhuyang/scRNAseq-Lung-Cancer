invisible(lapply(c("dplyr","magrittr","tidyr","openxlsx",#"Seurat","MAST","future",
                   "gplots"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#==============
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

#===========================================================
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
        csv_names = list.files("output/20220817/celltype.3",pattern = Cell_category[i],full.names = T)
        all_idx = which(opts$ident %in% Cell_category[i])
        idx <- gsub("output/20220817/celltype.3/","",csv_names) %>% gsub("-.*","",.) %>% as.integer()
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
