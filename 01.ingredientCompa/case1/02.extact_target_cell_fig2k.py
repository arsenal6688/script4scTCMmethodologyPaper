import scanpy as sc
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

# ====================== 路径 ======================
GLOBAL_PATH = "global_object_vasc_atlas_cxg.h5ad"
VASC_PATH   = "vas_for_cxg_240904.h5ad"
OUTPUT_PATH = "/script4scTCMmethodologyPaper/01.ingredientCompa/case1/02.fig2k_light.h5ad"

# ====================== 只保留血管细胞 ======================
print("Loading vascular data...")
adata_vasc = sc.read_h5ad(VASC_PATH)

vascular_target = [
    "aorta_coronary_ec",
    "art_ec_1",
    "art_ec_2",
    "brain_art_ec",
    "kidney_art_ec",
    "spleen_art_ec",
    "pul_art_ec",
    "cap_ec",
    "ven_ec_1",
    "ven_ec_2",
    "pul_ven_ec",
    "brain_ven_ec",
    "endocardial_ec",
    "littoral_EC",
    "lymphatic_EC"
]

# 筛选
adata_vasc = adata_vasc[adata_vasc.obs["ann_final"].isin(vascular_target)].copy()

# ====================== 统一命名 ======================
mapping = {
    "aorta_coronary_ec":"Aorta_coronary_ec",
    "art_ec_1":"Art_ec_1","art_ec_2":"Art_ec_2",
    "brain_art_ec":"Brain_art_ec","kidney_art_ec":"Kidney_art_ec","spleen_art_ec":"Spleen_art_ec",
    "pul_art_ec":"Pul_art_ec","cap_ec":"Capillary_ec",
    "ven_ec_1":"Ven_ec_1","ven_ec_2":"Ven_ec_2","pul_ven_ec":"Pul_ven_ec","brain_ven_ec":"Brain_ven_ec",
    "endocardial_ec":"Endocardial_ec","littoral_EC":"Littoral_EC","lymphatic_EC":"Lymphatic_EC"
}
adata_vasc.obs["cell_type_plot"] = adata_vasc.obs["ann_final"].map(mapping)

# ====================== 直接输出 ======================
adata_final = adata_vasc.copy()

print("\n===== 图2k 轻量版提取完成 =====")
print(adata_final.obs["cell_type_plot"].value_counts())

adata_final.write_h5ad(OUTPUT_PATH)
print(f"\n文件已保存：{OUTPUT_PATH}")
