# -*- coding: utf-8 -*-
import scanpy as sc
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

# ====================== 路径 ======================
GLOBAL_PATH = "final_ec_oct2024_cxg.h5ad"
OUTPUT_PATH = "/script4scTCMmethodologyPaper/01.ingredientCompa/case2/02.fig3i_light.h5ad"

# ====================== 仅保留 图3i + 数据集 共有的内皮细胞 ======================
target_celltypes = [
    # 动脉 EC
    "aorta_coronary_ec",
    "art_ec_1",
    "art_ec_2",
    "brain_art_ec",
    "kidney_art_ec",
    "spleen_art_ec",
    "pul_art_ec",

    # 静脉 EC
    "ven_ec_1",
    "ven_ec_2",
    "pul_ven_ec",
    "brain_ven_ec",
    "endocardial_ec",

    # 毛细血管 EC
    "aerocyte_ec",
    "pul_cap_ec",
    "pul_tip_cap_ec",
    "periportal_cap_ec",
    "pericentral_cap_ec",
    "kidney_cap_ec",
    "glomeruli_ec",
    "cap_ec",
    "endometrium_cap_ec",
    "myo_cap_ec",
    "adip_cap_ec",
    "blood_brain_barrier_ec"
]

# ====================== 加载数据并筛选 ======================
print("[1/3] 加载 final_ec_oct2024_cxg.h5ad ...")
adata = sc.read_h5ad(GLOBAL_PATH)

print(f"\n[2/3] 筛选图3i共有的内皮细胞（共{len(target_celltypes)}种）...")
adata_sub = adata[adata.obs["ann_ec_states"].isin(target_celltypes)].copy()

# ====================== 统一命名 ======================
mapping = {
    "aorta_coronary_ec": "Aorta_coronary_ec",
    "art_ec_1": "Art_ec_1",
    "art_ec_2": "Art_ec_2",
    "brain_art_ec": "Brain_art_ec",
    "kidney_art_ec": "Kidney_art_ec",
    "spleen_art_ec": "Spleen_art_ec",
    "pul_art_ec": "Pul_art_ec",
    "ven_ec_1": "Ven_ec_1",
    "ven_ec_2": "Ven_ec_2",
    "pul_ven_ec": "Pul_ven_ec",
    "brain_ven_ec": "Brain_ven_ec",
    "endocardial_ec": "Endocardial_ec",
    "aerocyte_ec": "Aerocyte_ec",
    "pul_cap_ec": "Pul_cap_ec",
    "pul_tip_cap_ec": "Pul_tip_cap_ec",
    "periportal_cap_ec": "Periportal_cap_ec",
    "pericentral_cap_ec": "Pericentral_cap_ec",
    "kidney_cap_ec": "Kidney_cap_ec",
    "glomeruli_ec": "Glomeruli_ec",
    "cap_ec": "Cap_ec",
    "endometrium_cap_ec": "Endometrium_cap_ec",
    "myo_cap_ec": "Myo_cap_ec",
    "adip_cap_ec": "Adip_cap_ec",
    "blood_brain_barrier_ec": "Blood_brain_barrier_ec"
}

adata_sub.obs["cell_type_plot"] = adata_sub.obs["ann_ec_states"].map(mapping)

# ====================== 标准化 ======================
print("\n[3/3] 数据标准化...")
sc.pp.normalize_total(adata_sub, target_sum=1e4)
sc.pp.log1p(adata_sub)
adata_sub.raw = adata_sub

# ====================== 输出统计（按图3i顺序展示）======================
print("\n===== 图3i 提取完成（仅共有的内皮细胞）=====")
# 按图3i横坐标顺序输出统计
order = [v for v in mapping.values()]
counts = adata_sub.obs["cell_type_plot"].value_counts().reindex(order)
print(counts)
print(f"\n✅ 总细胞数：{adata_sub.n_obs}")
print(f"✅ 总基因数：{adata_sub.n_vars}")

# ====================== 保存文件 ======================
adata_sub.write_h5ad(OUTPUT_PATH)
print(f"\n✅ 文件已保存至：{OUTPUT_PATH}")
