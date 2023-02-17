library(dplyr)
library(Seurat)
library(patchwork)
library(progress)
set.seed(1234)

# プログレスバーのための変数
n <- 100

cat("single cell RNA-seq解析を始めます\n")
# ディレクトリの設定

path <- "~/analysis/ykato/scRNAseq/data"
# ディレクトリの表示
setwd(path)
cat(list.files(), "\n")

# dataを格納したdirectoryを入力する
input.dir <- readline("ディレクトリを入力してください\n")

# Seuratオブジェクトの構築
# 必要なら、データを結合するバージョンも作成する
input.project <- readline("プロジェクト名を入力してください\n")
cat("データを読み込みます\n")
data <- Read10X(data.dir=input.dir)


# seurat objectの作成
cat("Seuratオブジェクトを作成します\n")
seurat.data <- CreateSeuratObject(counts = data, project = input.project, min.cells = 3, min.features = 200)
# 確認でseuratオブジェクトを確認する
cat(length(seurat.data@meta.data$nFeature_RNA))
cat("\n")

# mtが入るまえのデータであることの確認
cat("ミトコンドリア由来の遺伝子を除去します\n")
seurat.data[["percent_mt"]] <- PercentageFeatureSet(seurat.data, pattern="^MT-")
# mtが入ったことの確認
cat("ミトコンドリア由来の遺伝子の割合を追加しました\n")
cat("\n")

# violinplotで
cat("nFeatureRNA, nCountRNA, percentmtをpngで出力します\n")
png("VlnPlot_QC_before_filtering.png", width=800, height=500)
VlnPlot(seurat.data, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3)
dev.off()

# 各パラメータの入力を求める


# ここまでRdataで保存する
#save(seurat.data, file="seurat_data_01.Rdata")