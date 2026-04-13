import anndata as ad

# 读取两个文件
adata1 = ad.read_h5ad("macrophage_subtypes_manual.h5ad")
adata2 = ad.read_h5ad("cardiomyocyte_subtypes_manual.h5ad")

# 合并
adata_combined = ad.concat([adata1, adata2], merge="same")

# 保存为新文件
adata_combined.write_h5ad("/script4scTCMmethodologyPaper/02.herbCompa/case1/combined_macrophage_cardiomyocyte.h5ad")

print("✅ 合并完成！")
print("总细胞数：", adata_combined.n_obs)
print("细胞类型统计：")
print(adata_combined.obs['celltype_manual'].value_counts())
