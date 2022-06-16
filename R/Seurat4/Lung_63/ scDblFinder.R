library(scDblFinder)
library(SingleCellExperiment)
library(magrittr)
library(Seurat)
library(BiocSingular)
library(scran)
# ====== load single cell =============
object = readRDS(file = "data/Lung_63_20220606.rds")
meta.data = readRDS(file = "output/Lung_63_20220408_meta.data_v4.rds")
table(rownames(object@meta.data) == rownames(meta.data))
object@meta.data = meta.data
object[["RNA"]] = NULL
sce <- as.SingleCellExperiment(object)
#rm(object);
GC()
set.seed(100)

# Setting up the parameters for consistency with denoisePCA();
# this can be changed depending on your feature selection scheme.
dec <- modelGeneVarByPoisson(sce)
top.genes <- getTopHVGs(dec, prop=0.1)
length(top.genes)
dbl.dens <- computeDoubletDensity(sce, subset.row=top.genes, 
                                  d=ncol(reducedDim(sce)))
summary(dbl.dens)
sce$DoubletScore <- dbl.dens
dbl.calls <- doubletThresholding(data.frame(score=dbl.dens),
                                 method="griffiths", returnType="call")
summary(dbl.calls)
meta.data$DoubletScore = dbl.dens
meta.data$DoubletCall = dbl.calls
colnames(meta.data) %<>% sub("DoubletScore","DoubletScore.scDblFinder",.)
colnames(meta.data) %<>% sub("DoubletCall","DoubletCall.scDblFinder",.)

saveRDS(meta.data, file = "output/Lung_63_20220408_meta.data_v4.rds")
