from sctcm.d4c.TCMscore2cell import TCMscore2cell
import scanpy as sc

# 路径
INPUT_TXT = "d4c_TCMscore2cell-input_template.txt"
SC_TXT = "/scTCM/docs/templates/sc_io-input_template.txt"
ANNO_H5AD = "/script4scTCMmethodologyPaper/02.herbCompa/case1/combined_macrophage_cardiomyocyte.h5ad"
FORM_CONFIG = "d4c_TCMscore2cell-formula_herb_config_template.txt"
OUT = "/script4scTCMmethodologyPaper/02.herbCompa/case1"

# 1. 读取 adata 
adata = sc.read(ANNO_H5AD)
adata.obs['celltype_manual'] = adata.obs['major_celltype'].astype(str)  # 绘图分组用
#adata.obs['Group'] = adata.obs['group']                                # 小提琴图用
if 'X_mde' in adata.obsm:
    adata.obsm['X_umap'] = adata.obsm['X_mde']                         # UMAP 绘图用
adata.write(ANNO_H5AD)  
# =========================================================

if __name__ == "__main__":
    TCMscore2cell(
        input_txt=INPUT_TXT,
        sc_meta_txt=SC_TXT,
        anno_h5ad=ANNO_H5AD,
        formula_config_txt=FORM_CONFIG,
        output_dir=OUT
    )
