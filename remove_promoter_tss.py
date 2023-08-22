#!/usr/bin/python
import pandas as pd
import sys

crm_name = sys.argv[1] #"A549.active_predict.LR.bed"
df_crm = pd.read_csv(crm_name, header=None, index_col=None, sep="\t")
df_crm = df_crm.set_index([0,1,2])

df_tss = pd.read_csv("../ALL_CRMs_300.pos_site.keep_pos_site.n-1.score.sort_36.clean.closest_tss.csv", header=0, index_col=0)
df_tss_enhancer = df_tss[df_tss['8'] > 1000]
df_tss_enhancer = df_tss_enhancer.set_index(['0', '1', '2'])
df_crm_enhancer = df_crm[df_crm.index.isin(df_tss_enhancer.index)].reset_index()
df_crm_enhancer.to_csv(f"{crm_name}.enhancer", sep="\t", index=False, header=False)


