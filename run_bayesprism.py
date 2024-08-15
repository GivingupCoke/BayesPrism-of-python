import os
import pandas as pd
import pickle
from pybayesprism import process_input, prism, extract

#error
#Traceback (most recent call last):
#   File "/local/zhouxinyu/R_Bayesprism/workshop/GSE103224/run_bayesprism.py", line 17, in <module>
#     sc_dat_filtered_pc = process_input.select_gene_type(sc_dat_filtered, ["protein_coding"])
#   File "/home/zhouxinyu/.local/lib/python3.8/site-packages/pybayesprism/process_input.py", line 165, in select_gene_type
#     input_filtered = input.loc[:, selected_gene_idx.tolist()]
#IndexError: Boolean index has wrong length: 31737 instead of 31843

import os
import pandas as pd

def select_gene_type(input: pd.DataFrame, gene_type):
    assert all([g in ["protein_coding", "pseudogene", "lincRNA"] for g in gene_type])
    
    input_genes = input.columns
    gene_tab_path = os.path.join(os.path.split(__file__)[0], "txt", "gencode.v22.broad.category.txt")
    gene_list = pd.read_table(gene_tab_path, sep="\t", header=None)
    
    if sum(1 for gene in input_genes if gene[:3] == "ENS") > len(input_genes) * 0.8:
        print("ENSEMBL IDs detected.")
        input_genes = [gene.split(".")[0] for gene in input_genes]
        gene_list_7 = gene_list.iloc[:,7].tolist()
        gene_match = [gene_list_7.index(x) for x in input_genes if x in gene_list_7]
        gene_df = gene_list.iloc[gene_match, [7, 8]]
    else:
        print("Gene symbols detected. Recommend to use ENSEMBL IDs for more unique mapping.")
        gene_list_4 = gene_list.iloc[:,4].tolist()
        gene_match = [gene_list_4.index(x) for x in input_genes if x in gene_list_4]
        gene_df = gene_list.iloc[gene_match, [4, 8]]


    # 以下原代码改动了
    # gene_df.columns = ["gene_name", "category"]
    # selected_gene_idx = gene_df["category"].isin(gene_type)
    
    # print("number of genes retained in each category: ")
    # print(gene_df.loc[selected_gene_idx, "category"].value_counts())
    # input_filtered = input.loc[:, selected_gene_idx.tolist()]
    # return input_filtered


    # 修改后的，由于可能有的input_genes在 “gencode.v22.broad.category.txt” 文件中可能不能找到，
    # gene_df长度比input_genes长度要小
    # 最后input_filtered = input.loc[:, selected_gene_idx.tolist()]报错
    gene_df.columns = ["gene_name", "category"]

    # 创建一个包含所有列的布尔索引，初始值为 False
    selected_gene_idx = [False] * len(input_genes)

    # 创建一个字典用于快速查找 gene_df 中的基因名对应的类别
    gene_to_category = dict(zip(gene_df["gene_name"], gene_df["category"]))

    # 更新布尔索引
    for i, gene in enumerate(input_genes):
        if gene in gene_to_category and gene_to_category[gene] in gene_type:
            selected_gene_idx[i] = True

    # 检查生成的布尔索引长度是否正确
    assert len(selected_gene_idx) == len(input_genes), f"Boolean index has wrong length: {len(selected_gene_idx)} instead of {len(input_genes)}"

    # 打印每个类别中保留的基因数量
    print("Number of genes retained in each category: ")
    retained_genes = [input_genes[i] for i in range(len(input_genes)) if selected_gene_idx[i]]
    retained_categories = [gene_to_category[gene] for gene in retained_genes]
    retained_df = pd.DataFrame({"gene_name": retained_genes, "category": retained_categories})
    print(retained_df["category"].value_counts())

    # 根据布尔索引筛选 input 数据框的列
    input_filtered = input.loc[:, selected_gene_idx]
    return input_filtered




bk_dat = pd.read_csv("bkRNA.csv", sep=",", index_col=0)
sc_dat = pd.read_csv("scRNA.csv", sep=",", index_col=0)


cell_state_labels = pd.read_csv("cell.state.labels.csv", header=None).iloc[:,0].tolist()

cell_type_labels = pd.read_csv("cell.type.labels.csv", header=None).iloc[:,0].tolist()

sc_dat_filtered = process_input.cleanup_genes(sc_dat, "count.matrix", "hs", \
                  ["Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"], 5)
                  
sc_dat_filtered_pc = select_gene_type(sc_dat_filtered, ["protein_coding"])

my_prism = prism.Prism.new(reference = sc_dat_filtered_pc, 
                          mixture = bk_dat, input_type = "count.matrix", 
                          cell_type_labels = cell_type_labels, 
                          cell_state_labels = cell_state_labels, 
                          key = "tumor", 
                          outlier_cut = 0.01, 
                          outlier_fraction = 0.1)

bp_res = my_prism.run(n_cores = 36, update_gibbs = True)     

theta = extract.get_fraction(bp_res, "final", "type")
theta_cv = bp_res.posterior_theta_f.theta_cv
Z_tumor = extract.get_exp(bp_res, "type", "tumor")


with open("bp.pkl", "wb") as f:
    pickle.dump(bp_res, f)



