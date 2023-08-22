#!/usr/bin/python
from __future__ import print_function

import sys
import re
import pyBigWig
import numpy as np

def extract_binding_site(crm_file):

    binding_sites_set = set()
    with open(crm_file) as fin:
        for line in fin:
            line_split = line.strip().split("\t")
            binding_sites = [binding.split("|")[1] for binding in line_split[4].split(",")]
            binding_sites_set.update(binding_sites)

    return sorted(list(binding_sites_set))


def extract_conservation_score(pos_file, conservation_file, flag):
    bw = pyBigWig.open("/projects/zcsu_research/npy/CO_OCCUR/Network_Anlysis/POST_ANALYSIS/reference/"+conservation_file)

    conservation_score_list = []
    with open(pos_file) as fin:
        for line in fin:
            chrom, start, end = line.strip().split("\t")

            start = int(start)
            end = int(end)

            if flag == "gerp":
                chrom = chrom.replace("chr", "")
                
            values = bw.values(chrom, start, end)
            conservation_score_list.extend(values)

            
    bw.close()
    conservation_score = np.array(conservation_score_list)
    conservation_score = conservation_score[np.logical_not(np.isnan(conservation_score))]
    np.save(pos_file+".remove_nan."+flag, conservation_score)
    
def write_out_binding_site(crm_file, binding_site_set):

    with open(crm_file+".binding_site", "w") as fout:
        for binding_site in binding_site_set:
            binding_site = binding_site.replace(":", "\t").replace("-", "\t")
            new_line = binding_site+"\n"
            fout.write(new_line)
        
def main():
    crm_file = sys.argv[1]
    conservation_file = sys.argv[2]
    flag = sys.argv[3]
    #binding_site_set = extract_binding_site(crm_file)
    extract_conservation_score(crm_file, conservation_file, flag)

if __name__=="__main__":
    main()
