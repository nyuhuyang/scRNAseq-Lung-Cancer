library(dplyr)
library(tidyr)
library(magrittr)
library(readxl)
library(stringr)
library(openxlsx)
library(tidyverse)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive <- T)


# load data
meta.data <- readRDS( "output/Lung_63_20220408_meta.data_v5.rds")

df_samples <- readxl::read_excel("doc/20220406-samples metadata RS.xlsx", sheet <- "invivo")
df_samples <- as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
nrow(df_samples)

table(df_samples$sample %in% meta.data$orig.ident)
for(i in 1:length(df_samples$sample)){
    cells <- meta.data$orig.ident %in% df_samples$sample[i]
    print(df_samples$sample[i])
    print(table(cells))
    meta.data[cells,"group3"] <- df_samples$group3[i]
    meta.data[cells,"group2"] <- df_samples$group2[i]
    meta.data[cells,"group1"] <- df_samples$group1[i]
    meta.data[cells,"group0"] <- df_samples$group0[i]
}
meta.data$Doublets %<>% gsub("Doublet-High Confidence|Doublet-Low Confidence","Db",.)
#============ prepare annotation_df and group_df ============================
df_template <- readxl::read_excel("doc/63-sample Cell distribution per sample per group-YH.xlsx",
                                sheet <- "A")
Symbol <- unique(df_template$Symbol)
Groups <- colnames(df_template)[6:ncol(df_template)] %>% 
            strsplit(split <- " vs ") %>%
            unlist() %>%
            gsub("\\.\\.\\..*","",.) %>%
            unique %>% sort

#annotation_df
annotations <- readxl::read_excel("doc/Annotations/annotations.xlsx")
#group_list
df_samples <- readxl::read_excel("doc/20220406-samples metadata RS.xlsx", sheet <- "invivo")
colnames(df_samples) %<>% tolower()
group_df <- df_samples[,c("sample",paste0("group",3:0),"sex")] %>%
    pivot_longer(cols <- starts_with("group"))
colnames(group_df) <- c("sample","sex","group","symbol")
group_df$sex %<>% plyr::mapvalues(from <- c("M","F"),
                                  to <- c("male","female"))
group_df$`symbol-sex` <- paste0(group_df$symbol,"-",group_df$sex)
group_df <- as.data.frame(group_df[!is.na(group_df$symbol),
                                  c("sample","symbol","symbol-sex")])
group_df %<>% pivot_longer(cols <- starts_with("symbol"))

table(Groups %in% group_df$value)
Groups[!(Groups %in% group_df$value)]
group_df %<>% filter(value %in% Groups)

group_df$value %<>% factor(levels <- Groups)
group_df <- group_df[order(group_df$value),c("sample","value")]
group_list <- split(group_df ,f <- group_df$value) %>% 
    sapply(function(x) x$sample)

df_template <- select(df_template,-c("Order","Cell category","Family","superfamily")) %>%
                column_to_rownames("Symbol")

#============ A - Distribution per ALL cells ============================
annotations_df_A <- annotations[,c("celltype.3","celltype.2","celltype.1",
                                   "Family","Superfamily")] %>%
    pivot_longer(cols <- everything())
annotations_df_A$name %<>% factor(levels <- c("Superfamily","Family",
                                              "celltype.1","celltype.2","celltype.3"))
annotations_df_A <- annotations_df_A[order(annotations_df_A$name),]
annotations_df_A <- annotations_df_A[!duplicated(annotations_df_A$value),]
annotations_df_A %<>% as.data.frame()
Symbol[!(Symbol %in% annotations_df_A$value)]
annotations_df_A$value[!(annotations_df_A$value %in% Symbol)]
annotations_df_A %<>% rbind(data.frame("name" = "Doublets","value" = "Db"))

sample_size_list <- table(meta.data$orig.ident) %>% as.data.frame %>% as.list()
sample_size <- sample_size_list$Freq
names(sample_size) <- sample_size_list$Var1

Class <- sapply(meta.data,class)
Factor_class <- Class[Class == "factor"]
for(cl in grep("orig.ident|SCT_snn_",names(Factor_class),value =T, invert = T)){
    meta.data[,cl] %<>% as.character()
}

for(cl in colnames(annotations_df_A)){
    annotations_df_A[,cl] %<>% as.character()
}

# cell number
df_number_A <- select(df_template,-grep(" vs ",colnames(df_template)))


table(colnames(df_number_A) %in% Groups)
table(rownames(df_number_A) %in% Symbol)

for(m in 1:nrow(df_number_A)){
    cell_type_df <- annotations_df_A %>% filter(value %in% rownames(df_number_A)[m]) 
    for(n in 1:ncol(df_number_A)){
        group <- colnames(df_number_A)[n]
        num <-  sum((meta.data[,cell_type_df$name] %in% cell_type_df$value) &
                        (meta.data$orig.ident %in% group_list[[group]]))
        df_number_A[m,n] = num/sum(sample_size[group_list[[group]]])
    }
    Progress(m,nrow(df_number_A))
}
df_number_A <- df_number_A*100

# 2-tailed Mann-Whitney test
df_wix_A <- select(df_template,grep(" vs ",colnames(df_template)))[,1:15]
colnames(df_wix_A) %<>% gsub("\\.\\.\\..*","",.)


test_allcells <- function(df_res, df_data = meta.data, annotations_df, group_list, test = wilcox.test, paried = FALSE){
    for(m in 1:nrow(df_res)){
        cell_type_df <- annotations_df %>% filter(value %in% rownames(df_res)[m]) 
        df_data %>% filter(get(cell_type_df$name) %in% cell_type_df$value) -> sub_meta.data
        
        for(n in 1:ncol(df_res)){
            groups <- strsplit(colnames(df_res)[n],split = " vs ")[[1]]
            group1 <- group_list[[groups[1]]]
            group2 <- group_list[[groups[2]]]
            
            cell_numbers_list <- lapply(list(group1,group2), function(group){
                table(sub_meta.data$orig.ident) %>% as.data.frame.table() %>%
                    filter(Var1 %in% group) %>% .[,"Freq"]
            })
            if(any(sapply(cell_numbers_list,length) == 0) | all(sapply(cell_numbers_list, sd)  == 0)) {
                next
            } else {
                df_res[m,n] <- test(cell_numbers_list[[1]],cell_numbers_list[[2]],
                                             alternative = c("two.sided"),
                                             paried = paried, exact = TRUE) %>% .$"p.value"
            }
        }
        Progress(m,nrow(df_res))
    }
    return(df_res)
}
df_wix_A %<>% test_allcells(df_data = meta.data, annotations_df = annotations_df_A, group_list =group_list,
                            test = wilcox.test, paried = FALSE)

# 2-tailed t test
df_ttest_A <- select(df_template,grep(" vs ",colnames(df_template)))[,16:27]
colnames(df_ttest_A) %<>% gsub("\\.\\.\\..*","",.)

df_ttest_A %<>% test_allcells(df_data = meta.data, annotations_df = annotations_df_A, group_list =group_list,
                            test = t.test, paried = FALSE)

# 2-tailed paried Mann-Whitney test
df_paired_wix_A <- select(df_template,grep(" vs ",colnames(df_template)))[,28:33]
colnames(df_paired_wix_A) %<>% gsub("\\.\\.\\..*","",.)

df_paired_wix_A %<>% test_allcells(df_data = meta.data, annotations_df = annotations_df_A, group_list =group_list,
                              test = wilcox.test, paried = TRUE)


A_res = bind_cols(list(df_number_A,df_wix_A,df_ttest_A,df_paired_wix_A))

write.xlsx(A_res, file <- paste0(path,"A - Distribution per ALL cells.xlsx"),
           colNames = TRUE, rowNames = TRUE, borders = "surrounding")

#============ B - Distribution per ALL cells ============================
annotations_df_B1 <- annotations[,c("Superfamily")] %>%
    pivot_longer(cols <- everything()) %>%
    distinct(value, .keep_all = TRUE) %>%
    mutate(Superfamily = value, Family = value) %>%
    filter(value != "Un")

annotations_df_B2 <- annotations[,c("Family","Superfamily")] %>%
    pivot_longer(cols <- starts_with("Family")) %>%
    distinct(value, .keep_all = TRUE) %>%
    mutate(Family = value)

annotations_df_B3 <- annotations[,c("celltype.3","celltype.2","celltype.1",
                                    "Family","Superfamily")] %>%
    pivot_longer(cols <- starts_with("celltype")) %>%
    arrange(name) %>%
    distinct(value, .keep_all = TRUE) 

annotations_df_B <- bind_rows(list(annotations_df_B1[,colnames(annotations_df_B3)],
                                     annotations_df_B2[,colnames(annotations_df_B3)], 
                                     annotations_df_B3))
annotations_df_B %<>% as.data.frame()
annotations_df_B = annotations_df_B[!duplicated(annotations_df_B$value),]

Class <- sapply(meta.data,class)
Factor_class <- Class[Class == "factor"]
for(cl in grep("orig.ident|SCT_snn_",names(Factor_class),value =T, invert = T)){
    meta.data[,cl] %<>% as.character()
}

for(cl in colnames(annotations_df_B)){
    annotations_df_B[,cl] %<>% as.character()
}

# cell number
df_number_B <- select(df_template,-grep(" vs ",colnames(df_template)))
table(colnames(df_number_B) %in% Groups)
table(rownames(df_number_B) %in% Symbol)
for(m in 1:nrow(df_number_B)){
    cell_type_df <- annotations_df_B %>% filter(value %in% rownames(df_number_B)[m]) 
    meta.data %>% filter(Superfamily %in% cell_type_df$Superfamily) -> Super_meta.data
    
    sample_size_list <- table(Super_meta.data$orig.ident) %>% as.data.frame %>% as.list()
    sample_size <- sample_size_list$Freq
    names(sample_size) <- sample_size_list$Var1
    
    for(n in 1:ncol(df_number_B)){
        group <- colnames(df_number_B)[n]
        num <-  sum((Super_meta.data[,cell_type_df$name] %in% cell_type_df$value) &
                        (Super_meta.data$orig.ident %in% group_list[[group]]))
        df_number_B[m,n] = num/sum(sample_size[group_list[[group]]])
    }
    Progress(m,nrow(df_number_B))
}
df_number_B[c("Un","Db"),] = NA
df_number_B <- df_number_B*100

# 2-tailed Mann-Whitney test
df_wix_B <- select(df_template,grep(" vs ",colnames(df_template)))[,1:15]
colnames(df_wix_B) %<>% gsub("\\.\\.\\..*","",.)

test_superfamily <- function(df_res, df_data = meta.data, annotations_df, group_list, test = wilcox.test, paried = FALSE){
    for(m in 1:nrow(df_res)){
        cell_type_df <- annotations_df %>% filter(value %in% rownames(df_res)[m])
        meta.data %>% filter(Superfamily %in% cell_type_df$Superfamily) -> Super_meta.data
        
        sample_size_list <- table(Super_meta.data$orig.ident) %>% as.data.frame %>% as.list()
        sample_size <- sample_size_list$Freq
        names(sample_size) <- sample_size_list$Var1
        df_Freq <- table(Super_meta.data$orig.ident,Super_meta.data[,cell_type_df$name]) %>% as.data.frame.matrix()
        
        for(n in 1:ncol(df_res)){
            groups <- strsplit(colnames(df_res)[n],split = " vs ")[[1]]
            group1 <- group_list[[groups[1]]]
            group2 <- group_list[[groups[2]]]
            
            cell_numbers_list <- lapply(list(group1,group2), function(group){
                df_Freq[group,cell_type_df$value] /sum(sample_size[group])
            })
            if(any(sapply(cell_numbers_list,length) == 0) | all(sapply(cell_numbers_list, sd)  == 0)) {
                next
            } else {
                df_res[m,n] <- test(cell_numbers_list[[1]],cell_numbers_list[[2]],
                                             alternative = c("two.sided"),
                                             paried = paried, exact = FALSE) %>% .$"p.value"
            }
        }
        Progress(m,nrow(df_res))
    }
    return(df_res)
}
df_wix_B %<>% test_superfamily(df_data = meta.data, annotations_df = annotations_df_B, group_list =group_list,
                            test = wilcox.test, paried = FALSE)

# 2-tailed t test

df_ttest_B <- select(df_template,grep(" vs ",colnames(df_template)))[,16:27]
colnames(df_ttest_B) %<>% gsub("\\.\\.\\..*","",.)
df_ttest_B %<>% test_superfamily(df_data = meta.data, annotations_df = annotations_df_B, group_list =group_list,
                            test = t.test, paried = FALSE)
# 2-tailed paried Mann-Whitney test
df_paired_wix_B <- select(df_template,grep(" vs ",colnames(df_template)))[,28:33]
colnames(df_paired_wix_B) %<>% gsub("\\.\\.\\..*","",.)

df_paired_wix_B %<>% test_superfamily(df_data = meta.data, annotations_df = annotations_df_B, group_list =group_list,
                               test = wilcox.test, paried = TRUE)

B_res = bind_cols(list(df_number_B,df_wix_B,df_ttest_B,df_paired_wix_B))
write.xlsx(B_res, file <- paste0(path,"B - Distribution in SUPERFAMILY.xlsx"),
           colNames = TRUE, rowNames = TRUE, borders = "surrounding")



#============ C - Distribution in each FAMILY ============================
annotations_df_C <- annotations_df_B

# cell number
df_number_C <- select(df_template,-grep(" vs ",colnames(df_template)))
table(colnames(df_number_C) %in% Groups)
table(rownames(df_number_C) %in% Symbol)
for(m in 1:nrow(df_number_C)){
    cell_type_df <- annotations_df_C %>% filter(value %in% rownames(df_number_C)[m]) 
    if(!any(meta.data$Family %in% cell_type_df$Family)) next
    meta.data %>% filter(Family %in% cell_type_df$Family) -> Family_meta.data
    
    sample_size_list <- table(Family_meta.data$orig.ident) %>% as.data.frame %>% as.list()
    sample_size <- sample_size_list$Freq
    names(sample_size) <- sample_size_list$Var1
    
    for(n in 1:ncol(df_number_C)){
        group <- colnames(df_number_C)[n]
        num <-  sum((Family_meta.data[,cell_type_df$name] %in% cell_type_df$value) &
                        (Family_meta.data$orig.ident %in% group_list[[group]]))
        if(sum(sample_size[group_list[[group]]]) == 0) {
            df_number_C[m,n] = 0
            next
        } else {
            df_number_C[m,n] = num/sum(sample_size[group_list[[group]]])
            }
    }
    Progress(m,nrow(df_number_C))
}
df_number_C[c("Un","Db"),] = NA
df_number_C <- df_number_C*100

# 2-tailed Mann-Whitney test
df_wix_C <- select(df_template,grep(" vs ",colnames(df_template)))[,1:15]
colnames(df_wix_C) %<>% gsub("\\.\\.\\..*","",.)

test_family <- function(df_res, df_data = meta.data, annotations_df, group_list, test = wilcox.test, paried = FALSE,verbose =FALSE){
    for(m in 1:nrow(df_res)){
        cell_type_df <- annotations_df %>% filter(value %in% rownames(df_res)[m])
        if(!any(meta.data$Family %in% cell_type_df$Family)) next
        
        meta.data %>% filter(Family %in% cell_type_df$Family) -> Family_meta.data
        
        sample_size_list <- table(Family_meta.data$orig.ident) %>% as.data.frame %>% as.list()
        sample_size <- sample_size_list$Freq
        names(sample_size) <- sample_size_list$Var1
        df_Freq <- table(Family_meta.data$orig.ident,Family_meta.data[,cell_type_df$name]) %>% as.data.frame.matrix()
        
        for(n in 1:ncol(df_res)){
            if(verbose) print(paste("m =",m, "; n=",n))
            groups <- strsplit(colnames(df_res)[n],split = " vs ")[[1]]
            group1 <- group_list[[groups[1]]]
            group2 <- group_list[[groups[2]]]
            
            cell_numbers_list <- lapply(list(group1,group2), function(group){
                if(sum(sample_size[group]) == 0) return(sample_size[group])
                return(df_Freq[group,cell_type_df$value] /sum(sample_size[group]))
            })
            if(any(sapply(cell_numbers_list,length) == 0) | all(sapply(cell_numbers_list, sd)  == 0)) {
                next
            } else {
                df_res[m,n] <- test(cell_numbers_list[[1]],cell_numbers_list[[2]],
                                    alternative = c("two.sided"),
                                    paried = paried, exact = FALSE) %>% .$"p.value"
            }
        }
        Progress(m,nrow(df_res))
    }
    return(df_res)
}
df_wix_C %<>% test_family(df_data = meta.data, annotations_df = annotations_df_C, group_list =group_list,
                               test = wilcox.test, paried = FALSE)

# 2-tailed t test

df_ttest_C <- select(df_template,grep(" vs ",colnames(df_template)))[,16:27]
colnames(df_ttest_C) %<>% gsub("\\.\\.\\..*","",.)
df_ttest_C %<>% test_family(df_data = meta.data, annotations_df = annotations_df_C, group_list =group_list,
                                 test = t.test, paried = FALSE)
# 2-tailed paried Mann-Whitney test
df_paired_wix_C <- select(df_template,grep(" vs ",colnames(df_template)))[,28:33]
colnames(df_paired_wix_C) %<>% gsub("\\.\\.\\..*","",.)

df_paired_wix_C %<>% test_family(df_data = meta.data, annotations_df = annotations_df_C, group_list =group_list,
                                      test = wilcox.test, paried = TRUE)

C_res = bind_cols(list(df_number_C,df_wix_C,df_ttest_C,df_paired_wix_C))
write.xlsx(C_res, file <- paste0(path,"C - Distribution in each FAMILY.xlsx"),
           colNames = TRUE, rowNames = TRUE, borders = "surrounding")



