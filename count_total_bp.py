#!/usr/bin/python
from __future__ import print_function
import sys

def compute_total_bp(file_name):

    total_bp = 0
    with open(file_name) as fin:
        for line in fin:
            line_split = line.strip().split("\t")
            length = float(line_split[2]) - float(line_split[1])
            total_bp += length
    with open(file_name+".count_bp", "w") as fout:
        fout.write("total:"+str(total_bp))

def main():
    file_name = sys.argv[1]
    compute_total_bp(file_name)

if __name__=="__main__":
    main()
