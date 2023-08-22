#!/usr/bin/python
import sys
import glob

def main():
    file_name = sys.argv[1]
    chr_index = ["chr"+str(i+1) for i in range(22)] + ["chrX", "chrY"]
    file_out = "{}.filter".format(file_name)
    with open(file_out,"w") as fout:
        with open(file_name) as fin:
            for line in fin:
                chr_name = line.strip().split()[0]
                if chr_name in chr_index:
                    fout.write(line)
                    
if __name__=="__main__":
    main()
