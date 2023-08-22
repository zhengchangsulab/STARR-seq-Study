#!/usr/bin/python
import pandas as pd
import sys

def select_enhancer(file_name, df_silencer_ctcf):
    df_crm = pd.read_csv(file_name, header=None, index_col=None, sep="\t")
    df_crm = df_crm.set_index([0,1,2])
    df_crm_select = df_crm[~df_crm.index.isin(df_silencer_ctcf.index)]

    df_tss = pd.read_csv("../ALL_CRMs_300.pos_site.keep_pos_site.n-1.score.sort_36.clean.closest_tss.csv", header=0, index_col=0)
    df_tss_enhancer = df_tss[df_tss['8'] > 1000]
    df_tss_enhancer = df_tss_enhancer.set_index(['0', '1', '2'])
    
    df_crm_select = df_crm_select[df_crm_select.index.isin(df_tss_enhancer.index)].reset_index()
    df_crm_select.to_csv(f"{file_name}.enhancer", header=False, index=False, sep="\t")


def main():
    cell_name = sys.argv[1]

    ctcf = f"CTCF_folder/{cell_name}.intersect.ge_2.bed.merge.filter.fall_in_crm.crm.sort.uniq"
    silencer = f"{cell_name}.silencer.bed"

    df_ctcf = pd.read_csv(ctcf,header=None, index_col=None, sep="\t")
    df_silencer = pd.read_csv(silencer, header=None, index_col=None, sep="\t")
    df_silencer_clean = df_silencer.iloc[:,:3]
    df_silencer_ctcf = pd.concat([df_silencer_clean,df_ctcf], axis=0).drop_duplicates()
    df_silencer_ctcf = df_silencer_ctcf.set_index([0,1,2])

    active_crm = f"{cell_name}.active_predict.LR.bed"
    non_active_crm = f"{cell_name}.non_active_predict.LR.bed"

    select_enhancer(active_crm, df_silencer_ctcf)
    select_enhancer(non_active_crm, df_silencer_ctcf)


    
    
    
if __name__=="__main__":
    main()
