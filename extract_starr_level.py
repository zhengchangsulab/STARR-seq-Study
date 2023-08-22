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

def extract_conservation_score(pos_file, cell_name):

    bw_files = glob.glob(f"{cell_name}_STARR*.bigWig")
    if len(bw_files) == 1:
        bw = pyBigWig.open(bw_files[0])
        with open(f"{pos_file}.starr_level", "w") as fout:
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

    else:
        bw1 = pyBigWig.open(bw_files[0])
        bw2 = pyBigWig.open(bw_files[1])
        with open(f"{pos_file}.starr_level", "w") as fout:
            with open(pos_file) as fin:
                for line in fin:
                    line_strip = line.strip()
                    line_split = line_strip.split("\t")
                    chrom, start, end = line_split[:3]
                    
                    start = int(start)
                    end = int(end)
                    values_list = []
                    
                    values1 = np.array(bw1.values(chrom, start, end))
                    values1 = values1[np.logical_not(np.isnan(values1))]

                    values2 = np.array(bw2.values(chrom, start, end))
                    values2 = values2[np.logical_not(np.isnan(values2))]

                    values = np.array([values1, values2])
                    mean_level = np.mean(values)

                    output = f"{chrom}\t{start}\t{end}\t{mean_level}\n"
                    fout.write(output)
        bw1.close()
        bw2.close()

def main():
    cell_name = sys.argv[1]
    pos_file = sys.argv[2]
    extract_conservation_score(pos_file, cell_name)    
if __name__=="__main__":
    main()
