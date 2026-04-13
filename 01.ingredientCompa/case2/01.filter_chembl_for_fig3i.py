# -*- coding: utf-8 -*-
import pickle

# ====================== Figure3i 目标药物 CHEMBL ID ======================
TARGET_IDS = {
    "CHEMBL1492500",
    "CHEMBL2107874",
    "CHEMBL1118",
    "CHEMBL1743029",
    "CHEMBL2181927",
}

# ====================== 文件路径 ======================
INPUT_PKL  = "/scTCM/sctcm/tests/result/ChEMBL_parseChEMBLall-chembl_36_nested_drug_dict.pkl"
OUTPUT_PKL = "/script4scTCMmethodologyPaper/01.ingredientCompa/case2/01.chembl_36_fig3i_FINAL.pkl"

# ====================== 加载数据 ======================
with open(INPUT_PKL, "rb") as f:
    nested_dict = pickle.load(f)

drug_targets_merged = {}
result_dict = {}

# ====================== 第一步：精确匹配 + 合并重复药物 ======================
for atc_code, drug_map in nested_dict.items():
    for drug_key, targets in drug_map.items():
        chembl_id = drug_key.split("|")[0]
        if chembl_id in TARGET_IDS:
            if drug_key in drug_targets_merged:
                drug_targets_merged[drug_key] = list(set(drug_targets_merged[drug_key] + targets))
            else:
                drug_targets_merged[drug_key] = targets
                print(f"✅ 首次匹配: {drug_key}")

# ====================== 第二步：重新构建ATC嵌套字典 ======================
for atc_code, drug_map in nested_dict.items():
    keep = {}
    for drug_key in drug_map.keys():
        if drug_key in drug_targets_merged:
            keep[drug_key] = drug_targets_merged[drug_key]
    if keep:
        result_dict[atc_code] = keep

# ====================== 第三步：最终去重 ======================
final_dict = {}
seen_chembl = set()
for atc_code, drug_map in result_dict.items():
    final_keep = {}
    for drug_key, targets in drug_map.items():
        chembl_id = drug_key.split("|")[0]
        if chembl_id not in seen_chembl:
            seen_chembl.add(chembl_id)
            final_keep[drug_key] = targets
            print(f"✅ 最终保留: {drug_key}")
    if final_keep:
        final_dict[atc_code] = final_keep

# ====================== 保存 ======================
with open(OUTPUT_PKL, "wb") as f:
    pickle.dump(final_dict, f)

total_count = sum(len(drugs) for drugs in final_dict.values())
print(f"\n🎉 去重完成！Figure3i 最终唯一药物数量: {total_count}")
