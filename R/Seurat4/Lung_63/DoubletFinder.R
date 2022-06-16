########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
library(cowplot)
library(magrittr)
library(DoubletFinder)
require(fields)
require(parallel)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
########################################################################
#
#  2. DoubletFinder 
# 
# ######################################################################

# samples

object = readRDS(file = "data/Lung_63_20220606.rds")
meta.data = readRDS(file = "output/Lung_63_20220408_meta.data_v3.rds")
table(rownames(object@meta.data) == rownames(meta.data))
object@meta.data = meta.data
DefaultAssay(object) = "RNA"
object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20, nfeatures = 3000,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
object %<>% ScaleData()
object %<>% RunPCA()
object[["RNA"]]@scale.data = matrix(0,0,0) # prepare command list
object[["SCT"]] = NULL
object$Cell_label %<>% factor()
object_list <- SplitObject(object,split.by = "orig.ident")
#remove(object);GC()

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
Sys.setenv("OMP_NUM_THREADS" = 16)
npcs <- 100
sweep.res_list <- list()
for (i in 1:length(object_list)) {
    sweep.res_list[[i]] <- paramSweep_v4(object_list[[i]], PCs = 1:npcs, sct = F)
    Progress(i,length(object_list))
}
saveRDS(sweep.res_list,file = "output/Lung_63_20220606_sweep.res_list.rds")
sweep.res_list = readRDS("output/Lung_63_20220606_sweep.res_list.rds")
sweep_list <- lapply(sweep.res_list, function(x) summarizeSweep(x, GT = FALSE))
bcmvn_list <- lapply(sweep_list,find.pK)
# find histgram local maximam
find.localMaxima <- function(x) {
    # Use -Inf instead if x is numeric (non-integer)
    y <- diff(c(-.Machine$integer.max, x)) > 0L
    rle(y)$lengths
    y <- cumsum(rle(y)$lengths)
    y <- y[seq.int(1L, length(y), 2L)]
    if (x[[1]] == x[[2]]) {
        y <- y[-1]
    }
    which(x == max(x[y]))
}

(maximal_pk <- sapply(bcmvn_list,function(x) {
    as.numeric(as.character(x[find.localMaxima(x$BCmetric),"pK"]))
    }))
maximal_pk

# http://rstudio-pubs-static.s3.amazonaws.com/329613_f53e84d1a18840d5a1df55efb90739d9.html
qplot_2axis <- function(data,x = "pK", y1 = "MeanBC", y2 = "BCmetric"){
    if(class(data[,x]) == "factor") data[,x] <- as.numeric(as.character(data[,x]))
    data_y1 <- data[,y1]
    data_y2 <- data[,y2]
    a <- range(data_y1)
    b <- range(data_y2)
    scale_factor <- diff(a)/diff(b)
    data_y2 <- ((data_y2 - b[1]) * scale_factor) + a[1]
    trans <- ~ ((. - a[1]) / scale_factor) + b[1]
    
    g <- ggplot(data = data, aes_string(x = x, y = y1))+
        geom_line()+geom_point()+
        geom_point(aes(y = data_y2),colour = "blue")+
        geom_line(aes(y = data_y2),colour = "blue")+
        scale_y_continuous(name = y1,
                           sec.axis = sec_axis(trans=trans, name=y2))+
        theme(axis.text.y.right = element_text(color = "blue"))
    
    g
    
}
#qplot_2axis(data = bcmvn_list[[2]])

Multiplet_Rate <- function(object, numBatches = 1, num10xRuns = 1){
    
    numCellsRecovered = 1.0 * ncol(object)
    m = 4.597701e-06
    r = 0.5714286
    
    numCellsLoaded = numCellsRecovered / r
    multipletRate = m * numCellsLoaded / num10xRuns
    
    singletRate = 1.0 - multipletRate;
    numSinglet = singletRate * numCellsRecovered
    numMultiplet = numCellsRecovered - numSinglet
    numIdentMultiplet = numMultiplet * (numBatches - 1) / numBatches
    numNonIdentMultiplet = numMultiplet - numIdentMultiplet
    numCells = numSinglet + numNonIdentMultiplet
    
    return(numNonIdentMultiplet/numCells)
}
Multiplet_Rate(object_list[[1]])
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
for(i in 1:length(object_list)){
    print(paste("processing",unique(object_list[[i]]$orig.ident)))
    homotypic.prop <- modelHomotypic(object_list[[i]]@meta.data$Cell_label)
    nExp_poi <- round(Multiplet_Rate(object_list[[i]])*length(colnames(object_list[[i]])))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    object_list[[i]]$Cell_label %<>% factor(levels = levels(object$Cell_label))
    ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
    object_list[[i]] <- doubletFinder_v4(seu = object_list[[i]], PCs = 1:npcs,
                                         pN = 0.25, pK = maximal_pk[i],
                                         nExp = nExp_poi, reuse.pANN = FALSE,sct = F)#,annotations = object_list[[i]]$Cell_label)
    object_list[[i]] <- doubletFinder_v3(object_list[[i]], PCs = 1:npcs,
                                         pN = 0.25, pK = maximal_pk[i],
                                         nExp = nExp_poi.adj,
                                         reuse.pANN = grep("pANN",colnames(object_list[[i]]@meta.data),value = T),
                                         sct = FALSE)
    colName = colnames(object_list[[i]]@meta.data)
    colName[grep("DF.classifications",colName)] = c("Low_confident_doublets",
                                                    "High_confident_doublets")
    colnames(object_list[[i]]@meta.data) = colName
    Progress(i,length(object_list))
}

for(i in 1:length(object_list)){
    object_list[[i]]@meta.data$row.names = rownames(object_list[[i]]@meta.data)
}
meta.data_list <- lapply(object_list, function(x) {
    temp <- x@meta.data
    temp$row.names = rownames(temp)
    return(temp)
    })
meta.data = bind_rows(meta.data_list)
rownames(meta.data) = meta.data$row.names

object = readRDS(file = "data/Lung_63_20220606.rds")
meta.data = meta.data[rownames(object@meta.data),]
colnames(object@meta.data) %<>% gsub("Doublets","Doublets.old",.)
meta.data$doublets = gsub("Doublet","Doublet-Low Confidence",meta.data$Low_confident_doublets)
meta.data[meta.data$High_confident_doublets %in% "Doublet","doublets"] = "Doublet-High Confidence"
meta.data = cbind(object@meta.data,meta.data$doublets)
colnames(meta.data)[ncol(meta.data)] = "Doublets"
table(meta.data$Doublets)


object@meta.data = meta.data
saveRDS(meta.data, file = "output/Lung_63_20220408_meta.data_v4.rds")

#  cell number / cell percentage
meta.data = readRDS(file = "output/Lung_63_20220408_meta.data_v4.rds")
meta.data %<>% subset(Doublets == "Singlet")
table(meta.data$Doublets)

df_samples <- readxl::read_excel("doc/20220406-samples metadata RS.xlsx", sheet = "RS in vivo metadata")
df_samples = as.data.frame(df_samples)
df_samples = df_samples[,c("Sample","Subject ID","Age", "Sex","Ancestry","Ever smoked",
                           "Smoking status","Region","Condition",
                           "Study group")]
df_samples$`Study group` %<>% factor(levels = c("P-norm","D-norm","T-norm","L-norm","L-COPD","L-COPD-Dex","D-COPD",
                                                "L-IPF","D-IPF","L-Ad","L-Sq","L-Ad-Sq"))
df_table1 <- table(meta.data$orig.ident, meta.data$Cell_label) %>% as.data.frame.matrix()
df_table2 <- table(meta.data$orig.ident, meta.data$Superfamily) %>% as.data.frame.matrix()
table(df_samples$Sample == rownames(df_table1))
table(df_samples$Sample == rownames(df_table2))

df_sample1 = cbind(df_samples, df_table2[df_samples$Sample,], df_table1[df_samples$Sample,])

prop_table1 <- table(meta.data$orig.ident, meta.data$Cell_label) %>% 
    prop.table(1) %>% as.data.frame.matrix()
prop_table2 <- table(meta.data$orig.ident, meta.data$Superfamily) %>% 
    prop.table(1) %>% as.data.frame.matrix()
table(df_samples$Sample == rownames(prop_table1))
table(df_samples$Sample == rownames(prop_table2))

df_sample2 = cbind(df_samples, prop_table2[df_samples$Sample,], prop_table1[df_samples$Sample,])

df_sample1 = df_sample1[order(df_sample1$`Study group`),]
df_sample2 = df_sample2[order(df_sample2$`Study group`),]

df_list <- list("cell number" = df_sample1, "cell percentage" = df_sample2)
openxlsx::write.xlsx(df_list, file =  paste0(path,"20220615_63samples_cell number per cell type per sample.xlsx"),
                     colNames = TRUE,rowNames = F,borders = "surrounding")


#  Doubltes
meta.data = readRDS(file = "output/Lung_63_20220408_meta.data_v4.rds")

n_table3 <- table(meta.data$Cell_label, meta.data$Doublets) %>% 
    as.data.frame.matrix()
prop_table3 <- table(meta.data$Cell_label, meta.data$Doublets) %>% 
    prop.table(1) %>% as.data.frame.matrix()
n_table3 %<>% cbind(prop_table3)
openxlsx::write.xlsx(n_table3, file =  paste0(path,"20220615_63samples_Doublets.xlsx"),
                     colNames = TRUE,rowNames = T,borders = "surrounding")
