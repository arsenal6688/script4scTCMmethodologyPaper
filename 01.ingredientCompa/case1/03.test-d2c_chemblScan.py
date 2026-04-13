# -*- coding: utf-8 -*-
from sctcm.d2c import chemblScan
import warnings
warnings.filterwarnings("ignore")

# ====================== 配置路径 ======================
ANNDATA_PATH    = "/script4scTCMmethodologyPaper/01.ingredientCompa/case1/02.fig2k_light.h5ad"
CHEMBL_PKL      = "/script4scTCMmethodologyPaper/01.ingredientCompa/case1/01.chembl_36_fig2k_FINAL.pkl"
OUTPUT_DIR      = "/script4scTCMmethodologyPaper/01.ingredientCompa/case1"
GROUP_BY        = "cell_type_plot"

# ====================== 执行分析 ======================
if __name__ == "__main__":
    chemblScan.run_chembl_drug2cell(
        h5ad_path=ANNDATA_PATH,
        drug_dict_pkl=CHEMBL_PKL,
        output_dir=OUTPUT_DIR,
        groupby=GROUP_BY
    )
