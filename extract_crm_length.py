#!/usr/bin/python
from __future__ import print_function
import sys
import numpy as np

def extract_length(crm_name):

    crm_length_list = []
    with open(crm_name) as fin:
        for line in fin:
            line_split = line.strip().split("\t")
            length = int(line_split[2]) - int(line_split[1])
            crm_length_list.append(length)


    crm_length = np.array(crm_length_list)
    print(np.mean(crm_length))
    print(np.median(crm_length))

    np.save(crm_name+".length", crm_length)

            
            
def main():
    crm_name = sys.argv[1]
    extract_length(crm_name)


if __name__=="__main__":
    main()
