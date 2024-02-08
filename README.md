# scRNAseq
Seurat, Scanpyでsingle cell RNAseq解析を行う方法についてまとめています。
再解析ようなので、データはすでに処理されている前提です。

# Scanpy

## インストール
```
conda install -c conda-forge scanpy python-igraph leidenalg
# pip install scanpy
```

## resultのファイルの設定
```
results_file = 'write/pbmc3k.h5ad'
```
基本的にh5adファイルにanndataをまとめるのでファイル名を指定しておくと良い

### 保存したファイルを呼び出す場合
```
anndata = sc.read_h5ad("ファイル名")
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

# mtxのファイルがあればこちらで
```
adata = sc.read_10x_mtx(
    'data/filtered_gene_bc_matrices/hg19/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
```
