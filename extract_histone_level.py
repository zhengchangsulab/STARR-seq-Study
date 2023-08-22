#!/usr/bin/python
import sys
import re
import pyBigWig
import numpy as np
import glob

def extract_binding_site(crm_file):

    binding_sites_set = set()
    with open(crm_file) as fin:
        for line in fin:
            line_split = line.strip().split("\t")
            binding_sites = [binding.split("|")[1] for binding in line_split[4].split(",")]
            binding_sites_set.update(binding_sites)

    return sorted(list(binding_sites_set))

def extract_conservation_score(pos_file, cell_name, bw_name):

    bw = pyBigWig.open(bw_name)
    marker = bw_name.split(".")[1]
    with open(f"{pos_file}.{marker}.level", "w") as fout:
        with open(pos_file) as fin:
            for line in fin:
                line_strip = line.strip()
                line_split = line_strip.split("\t")
                chrom, start, end = line_split[:3]
                
                start = int(start)
                end = int(end)
                
                values = np.array(bw.values(chrom, start, end))
                values = values[np.logical_not(np.isnan(values))]
                
                mean_level = np.mean(values)
                line_strip = line.strip()
                output = f"{chrom}\t{start}\t{end}\t{mean_level}\n"
                fout.write(output)
    bw.close()

def main():
    cell_name = sys.argv[1]
    pos_file = sys.argv[2]
    histone_name = sys.argv[3]
    extract_conservation_score(pos_file, cell_name, histone_name)    
if __name__=="__main__":
    main()
