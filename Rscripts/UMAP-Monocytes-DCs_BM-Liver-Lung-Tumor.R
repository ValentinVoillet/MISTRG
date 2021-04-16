####---- Bone Marrow, Liver, Lung and Tumor - Monocytes/DCs - MISTRG
# install.packages(here("Rpackage", "Seurat_2.3.3.tar.gz"), repos = NULL, type = "source", lib = here("Rpackage"))
library(here)
library(Seurat, lib.loc = here("Rpackage")) # Seurat - v2.3.3
library(tidyverse)

###--- Analysis of DCs
seuratObj <- readRDS(here("output", "seuratObj_BM-Liver-Lung_Tumor.rds"))

#- Set-up a seurat object
seuratObj <- SubsetData(object = seuratObj, ident.use = c(1, 4, 9, 10, 11, 12), do.clean = TRUE)
seuratObj <- CreateSeuratObject(raw.data = seuratObj@raw.data, min.cells = 1, min.genes = 0) 

#- Add meta.data
mito.genes <- grep(pattern = "^MT-", x = rownames(x = seuratObj@data), value = TRUE)
percent.mito <- Matrix::colSums(seuratObj@raw.data[mito.genes, ]) / Matrix::colSums(seuratObj@raw.data)
seuratObj@meta.data$percent.mito <- percent.mito
seuratObj@meta.data$Tissue <- plyr::mapvalues(x = sapply(str_split(rownames(seuratObj@meta.data), "[.]"), function(x) x[2]), from = 1:4, to = c("Bone Marrow", "Liver", "Lung", "Tumor"))
seuratObj@meta.data %>% View("meta.data")

#- Normalization
seuratObj <- NormalizeData(object = seuratObj, normalization.method = "LogNormalize", scale.factor = 10000)

#- Highly variable features
seuratObj <- FindVariableGenes(object = seuratObj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = .2, x.high.cutoff = 5, y.cutoff = .5)

#- Scaling
seuratObj <- ScaleData(object = seuratObj, vars.to.regress = c("nUMI", "percent.mito"))

#- Clustering
seuratObj <- RunPCA(object = seuratObj, pc.genes = seuratObj@var.genes, pcs.compute = 100)
ElbowPlot(object = seuratObj, ndims = 50)
seuratObj <- FindClusters(object = seuratObj, reduction.type = "pca", dims.use = 1:20, resolution = .75, k.param = 25, print.output = TRUE, save.SNN = TRUE)

#- UMAP
seuratObj <- RunUMAP(object = seuratObj, reduction.use = "pca", dims.use = 1:20, min_dist = .75)
DimPlot(object = seuratObj, reduction = "umap", group.by = "Tissue", pt.size = .5)
DimPlot(object = seuratObj, reduction = "umap", group.by = "res.0.75", pt.size = .5)

#- Biomarkers
markers.DCs <- FindAllMarkers(object = seuratObj,  only.pos = TRUE, min.pct = .10, logfc.threshold = .25, return.thresh = .01, test.use = "MAST", latent.vars = "nUMI")

#- Output
saveRDS(seuratObj, here("output", "seuratObj_Monocytes-DCs_BM-Liver-Lung-Tumor.rds"))

#- Seurat v4
seuratObj.update <- UpdateSeuratObject(object = seuratObj)
data.table::data.table(cell.barcode = rownames(seuratObj.update@meta.data),
                       seuratObj.update@meta.data, 
                       Embeddings(seuratObj.update, "umap")) %>%
  write_csv(here("output", "meta-data_Monocytes-DCs_BM-Liver-Lung-Tumor.csv"))

