########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
# conda activate r4.1.1
#devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","sctransform","harmony","magrittr"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)


########################################################################
#
#  1 Seurat Alignment
#
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# read sample summary list
df_samples <- readxl::read_excel("output/20220406/20220406_Lung_2_dexamethasone.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
df_samples %<>% filter(sample %in% c("WC_30_T_Adj_Cont","WC_30_T_Adj_Dex"))
nrow(df_samples)
#======1.2 load  Seurat =========================
object = readRDS(file = "data/Lung_2_dexamethasone_20220406.rds")

table(df_samples$sample %in% object$orig.ident)
meta.data = object@meta.data
table(df_samples$sample %in% meta.data$orig.ident)
for(i in 1:length(df_samples$sample)){
    cells <- meta.data$orig.ident %in% df_samples$sample[i]
    print(df_samples$sample[i])
    print(table(cells))
    meta.data[cells,"region"] = df_samples$region[i]
    meta.data[cells,"subject"] = df_samples$subject[i]
    meta.data[cells,"subject id"] = df_samples$`subject id`[i]
    meta.data[cells,"age"] = df_samples$age[i]
    meta.data[cells,"sex"] = df_samples$sex[i]
    meta.data[cells,"ancestry"] = df_samples$ancestry[i]
    meta.data[cells,"ever smoked"] = df_samples$`ever smoked`[i]
    meta.data[cells,"smoking status"] = df_samples$`smoking status`[i]
    meta.data[cells,"condition"] = df_samples$condition[i]
    meta.data[cells,"study group"] = df_samples$`study group`[i]
    meta.data[cells,"project"] = df_samples$project[i]
    meta.data[cells,"date"] = as.character(df_samples$date[i])
    meta.data[cells,"Mean.Reads.per.Cell"] = df_samples$mean.reads.per.cell[i]
    meta.data[cells,"Number.of.Reads"] = df_samples$number.of.reads[i]
    meta.data[cells,"Sequencing.Saturation"] = df_samples$sequencing.saturation[i]
}
meta.data$orig.ident %<>% factor(levels = df_samples$sample)
table(rownames(object@meta.data) == rownames(meta.data))
table(colnames(object) == rownames(meta.data))
object@meta.data = meta.data
# Determine the ‘dimensionality’ of the dataset  =========
npcs = 100

DefaultAssay(object) <- "RNA"
object %<>% NormalizeData()
object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20, nfeatures = 2000,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
object %<>% ScaleData(verbose = FALSE)
object %<>% RunPCA(npcs = 100, verbose = FALSE)
jpeg(paste0(path,"ElbowPlot_RNA.jpeg"), units="in", width=10, height=7,res=600)
print(ElbowPlot(object,ndims = 100))
dev.off()
object <- JackStraw(object, num.replicate = 20,dims = npcs)
object <- ScoreJackStraw(object, dims = 1:npcs)

for(i in 0:9){
    a = i*10+1; b = (i+1)*10
    jpeg(paste0(path,"JackStrawPlot_",a,"_",b,".jpeg"), units="in", width=10, height=7,res=600)
    print(JackStrawPlot(object, dims = a:b))
    dev.off()
    Progress(i, 9)
}
p.values = object[["pca"]]@jackstraw@overall.p.values
print(npcs <- max(which(p.values[,"Score"] <=0.05)))
npcs = 61

#======1.6 Performing SCTransform and integration =========================

format(object.size(object),unit = "GB")
options(future.globals.maxSize= object.size(object)*10)
object %<>% SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = TRUE)
DefaultAssay(object) <- "SCT"
object %<>% ScaleData(verbose = FALSE)
object %<>% RunPCA(verbose = T,npcs = npcs)

jpeg(paste0(path,"S1_ElbowPlot_SCT.jpeg"), units="in", width=10, height=7,res=600)
ElbowPlot(object, ndims = npcs)
dev.off()

saveRDS(object, file = "data/Lung_2_dexamethasone_20220406.rds")

#======1.8 UMAP from harmony and pca =========================
DefaultAssay(object) = "SCT"

npcs = 61
jpeg(paste0(path,"S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(object %<>% RunHarmony.1(group.by = "orig.ident", dims.use = 1:npcs,
                                     theta = 2, plot_convergence = TRUE,
                                     nclust = 50, max.iter.cluster = 100))
dev.off()

object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)
#system.time(object %<>% RunTSNE(reduction = "harmony", dims = 1:npcs))

object[["harmony.umap"]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                                 key = "harmonyUMAP_", assay = DefaultAssay(object))
#object[["harmony.tsne"]] <- CreateDimReducObject(embeddings = object@reductions[["tsne"]]@cell.embeddings,
#                                                 key = "harmonytSNE_", assay = DefaultAssay(object))

object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
#system.time(object %<>% RunTSNE(reduction = "pca", dims = 1:npcs))
object %<>% FindNeighbors(reduction = "umap",dims = 1:2)
resolutions = c( 0.01, 0.1, 0.2, 0.5,0.8)
for(i in 1:length(resolutions)){
    object %<>% FindClusters(resolution = resolutions[i], algorithm = 1)
    Progress(i,length(resolutions))
}

saveRDS(object, file = "data/Lung_2_dexamethasone_20220406.rds")

# ======1.8.5 Unimodal UMAP Projection =========================
object = readRDS("data/Lung_2_dexamethasone_20220406.rds")
DefaultAssay(object) = "SCT"
object[["raw.umap"]] <- CreateDimReducObject(embeddings = object@reductions[["umap"]]@cell.embeddings,
                                             key = "rawUMAP_", assay = DefaultAssay(object))


object.reference = readRDS("data/Lung_SCT_30_20210831.rds")
meta.data = readRDS("output/20211209/meta.data_azimuth_subtype.rds")
table(rownames(object.reference@meta.data) == rownames(meta.data))
object.reference@meta.data = meta.data
object.reference[["orig.umap"]] <- CreateDimReducObject(embeddings = object.reference@reductions[["umap"]]@cell.embeddings,
                                                 key = "origUMAP_", assay = DefaultAssay(object.reference))
object.reference %<>% RunUMAP(reduction = "pca", dims = 1:100,return.model = TRUE)

embedding = object.reference[["orig.umap"]]@cell.embeddings
colnames(embedding) = c("UMAP_1","UMAP_2")
object.reference[["umap"]]@cell.embeddings = embedding
colnames(embedding) =NULL
object.reference[["umap"]]@misc$model$embedding = embedding


UMAPPlot.1(object.reference,group.by = "Cell_subtype",do.print = T, raster=FALSE)

#length(setdiff(colnames(object),colnames(object_reference)))
#newCells <- setdiff(colnames(object),colnames(object.reference))
#object.query = object[,newCells]


anchors <- FindTransferAnchors(reference = object.reference, query = object,
                               dims = 1:npcs, reference.reduction = "pca")
object <- MapQuery(anchorset = anchors,
                     reference = object.reference, query = object,
                     refdata = list(Cell_subtype = "Cell_subtype"),
                     reference.reduction = "pca", reduction.model = "umap")
colnames(object@meta.data) %<>% gsub("predicted.","",.)


meta.data = object.reference@meta.data
meta.data = meta.data[!duplicated(meta.data$Cell_subtype),]
object$Cell_subtype.colors = plyr::mapvalues(object$Cell_subtype,
                                             from = meta.data$Cell_subtype,
                                             to = meta.data$Cell_subtype.colors)
object$Cell_type = plyr::mapvalues(object$Cell_subtype,
                                   from = meta.data$Cell_subtype,
                                   to = meta.data$Cell_type)
object$Family = plyr::mapvalues(object$Cell_subtype,
                                from = meta.data$Cell_subtype,
                                to = meta.data$Family)
object$Superfamily = plyr::mapvalues(object$Cell_subtype,
                                     from = meta.data$Cell_subtype,
                                     to = meta.data$Superfamily)
saveRDS(object, file = "data/Lung_2_dexamethasone_20220406.rds")
