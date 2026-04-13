# -*- coding: utf-8 -*-
from sctcm.d2c import ingreScanCell
import warnings
warnings.filterwarnings("ignore")

# ====================== 配置路径 ======================
ANNDATA_PATH    = "/script4scTCMmethodologyPaper/01.ingredientCompa/case2/02.fig3i_light.h5ad"
INGRED_PKL      = "/scTCM/sctcm/tests/result/tcm_ing2targetall.pkl"
OUTPUT_DIR      = "/script4scTCMmethodologyPaper/01.ingredientCompa/case2"
GROUP_BY        = "cell_type_plot"

# ====================== 执行分析 ======================
if __name__ == "__main__":
    ingreScanCell.run_chembl_drug2cell(
        h5ad_path=ANNDATA_PATH,
        drug_dict_pkl=INGRED_PKL,
        output_dir=OUTPUT_DIR,
        groupby=GROUP_BY
    )
