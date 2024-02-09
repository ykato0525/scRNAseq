"""
前処理したデータから細胞株を選んで解析するscriptです。
"""

import numpy as np
import pandas as pd
import scanpy as sc
import sys

file_path = sys.argv[1]
adata = sc.read_h5ad(file_path)

selected_cells = adata[adata.obs['subtype_new'] == 'cancer cells'] # 一般化する場合には、これらも変数にする
gene_name = sys.argv[2]

# グラフの描画
# 発現量については指定するのが大変なので、標準化した値を使用する
# 色もdefaultの設定にして、変更可能にするとよい
sc.pl.umap(adata, color=gene_name,color_map = 'YlGnBu', use_raw=False)
