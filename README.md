This is a Python implementation of [BayesPrism v2.2](https://github.com/Danko-Lab/BayesPrism).

See [tutorial.md](https://github.com/ziluwang829/pyBayesPrism/blob/main/tutorial.md)
for usage. 


However, Something went wrong when I ran the code, on the select_gene_type function in the process_input.py file.

Error display:
Traceback (most recent call last):
    File "/local/zhouxinyu/R_Bayesprism/workshop/GSE103224/run_bayesprism.py", line 17, in <module>
    sc_dat_filtered_pc = process_input.select_gene_type(sc_dat_filtered, ["protein_coding"])
    File "/home/zhouxinyu/.local/lib/python3.8/site-packages/pybayesprism/process_input.py", line 165, in select_gene_type
    input_filtered = input.loc[:, selected_gene_idx.tolist()]
IndexError: Boolean index has wrong length: 31737 instead of 31843


May be due to some input_genes in "gencode. V22. Broad. Category. TXT" file may not be able to find, gene_df length is smaller than input_genes length, Finally, input_filtered = input.loc[:, selected_gene_idx.tolist()] displays an error. So I modified the code slightly, as shown in the run_bayesprism.py file.
