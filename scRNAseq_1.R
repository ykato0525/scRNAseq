library(dplyr)
library(Seurat)
library(patchwork)
set.seed(1234)

cat("single cell RNA-seq解析を始めます\n")
# ディレクトリの設定
# dataを格納したdirectoryを入力する
cat("ディレクトリを入力してください")
input.dir <- readLines(stdin(), n=1)

# Seuratオブジェクトの構築
# 必要なら、データを結合するバージョンも作成する
input.project <- readline("プロジェクト名を入力してください")
data <- Read10X(data.dir="inpur.dir")
# ここはdefaultでこの値を使用して良さそう
seurat.data <- CreateSeuratObject(counts = data, project = input.project, min.cells = 3, min.features = 200)
# 確認でseuratオブジェクトを確認する
length(seurat.data@meta.data$nFeature_RNA)

# mtが入るまえのデータであることの確認
head(seurat.data@meta.data)
seurat.data[["percent_mt"]] <- PercentageFeatureSet(seurat.data, pattern="^mt-")
# mtが入ったことの確認
cat("ミトコンドリア由来の遺伝子の割合を追加しました")
head(seurat.data@meta.data)

# violinplotで
cat("nFeatureRNA, nCountRNA, percentmtをpngで出力します")
png("VlnPlot_QC_before_filtering.png", width=800, height=500)
VlnPlot(A549, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# 各パラメータの入力を求める


# ここまでRdataで保存する
save(seurat.data, file="seurat_data_01.Rdata")