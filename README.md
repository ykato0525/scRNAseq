# scRNAseq
Scanpyでsingle cell RNAseq解析を行う方法についてまとめています。
再解析用なので、データはすでに処理されている前提です。

# Scanpy

## インストール
```
conda install -c conda-forge scanpy python-igraph leidenalg
# pip install scanpy
```


## コードの使いかた
```
python3 cell_type_gene_plot.py <<ファイルのパス>> <<遺伝子名>>
```
