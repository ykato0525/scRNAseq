# scRNAseq
Seurat, Scanpyでsingle cell RNAseq解析を行う方法についてまとめています。
再解析ようなので、データはすでに処理されている前提です。

# Scanpy

## インストール
```
conda install -c conda-forge scanpy python-igraph leidenalg
# pip install scanpy
```

## ライブラリのロードと設定の確認

```
import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
```
## データの読み込み
### resultのファイルの設定
```
results_file = 'write/pbmc3k.h5ad'
```
基本的にh5adファイルにanndataをまとめるのでファイル名を指定しておくと良い

### 保存したファイルを呼び出す場合
```
adata = sc.read_h5ad("ファイル名")
```


### ファイルの読み込み
```
adata = sc.read_10x_mtx(
    'data/filtered_gene_bc_matrices/hg19/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
```

## 前処理
```
sc.pl.highest_expr_genes(adata, n_top=20, )
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
```

### ミトコンドリア由来の遺伝子の割合で除去
```
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
# データの確認
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
# データの除去
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
```
### データの標準化
```
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```


# anndataの構造について
```
adata.X # 遺伝子発現の疎行列
adata.obs # 細胞のメタデータ
```

## 細胞を取り出してプロットしたいとき

