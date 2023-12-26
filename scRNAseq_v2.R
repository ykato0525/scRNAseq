library(tidyverse)
library(Seurat)
library(patchwork)

setwd("~/project/CGAS_BRCA/Fig2_scRNAseq/Wu_etal_2021_BRCA_scRNASeq/script")
data_dir = "~/project/CGAS_BRCA/Fig2_scRNAseq/Wu_etal_2021_BRCA_scRNASeq/Wu_etal_2021_BRCA_scRNASeq/"
data <- Read10X(data.dir=data_dir,
                gene.column = 1,
                cell.column = 1,
                unique.features = TRUE,
                strip.suffix = FALSE)
s_obj <- CreateSeuratObject(counts = data, project = "breast_cancer_atlas", min.cells = 3, min.features = 200)
s_obj


# メタデータからバーコードのリストを作成する
metadata = read.csv("../Wu_etal_2021_BRCA_scRNASeq/metadata.csv", header = 1, row.names = 1)
celltype_unique = unique(metadata$celltype_major)
print(celltype_unique)
metadata.cancer = metadata[metadata$celltype_major=="Cancer Epithelial",]
label.cancer = row.names(metadata.cancer)
s_obj.cancer = subset(s_obj, cells=label.cancer)



s_obj.cancer[["percent.mt"]] <- PercentageFeatureSet(s_obj, pattern = "^MT-")
VlnPlot(s_obj.cancer, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

s_sub <- subset(s_obj.cancer, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 3)

s_FVF <- FindVariableFeatures(s_sub, selection.method = "vst", nfeatures = 2000)

s_nor <- NormalizeData(s_FVF, normalization.method = "LogNormalize", scale.factor = 10000)

all.genes <- rownames(s_nor)
s_scale <- ScaleData(s_nor, features = all.genes)

s_scale <- RunPCA(s_scale, features = VariableFeatures(object = s_scale))
print(s_scale[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(s_scale, dims = 1:2, reduction = "pca")
ElbowPlot(s_scale)
s_scale <- FindNeighbors(s_scale, dims = 1:10)
s_scale <- FindClusters(s_scale, resolution = 0.5)
s_scale <- RunUMAP(s_scale, dims = 1:10)
DimPlot(s_scale, reduction = "umap")
FeaturePlot(s_scale, features = c("TMEM173", "DDX58", "MB21D1", "IL6"))
