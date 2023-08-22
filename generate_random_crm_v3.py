#!/usr/bin/python
import numpy as np
import sys


def build_dataset_dict():
    chrom_peak_dict = {}
    with open("../TF_all_sort_motifs_unmerge.1000.bed.cat.clean.sort.merge") as fin:
        for line in fin:
            chrom, start, end = line.strip().split("\t")
            try:
                chrom_peak_dict[chrom].append((int(start), int(end), int(end) - int(start)))
            except:
                chrom_peak_dict[chrom] = [(int(start), int(end), int(end)-int(start))]

    return chrom_peak_dict


def get_longest_peak_length_in_dataset_chrom(chrom_peak_dict):

    dataset_chrom_longest_peak = {}
    for chrom in chrom_peak_dict.keys():
        peak_info_list = chrom_peak_dict[chrom]
        longest_peak_length = max([peak_length for _, _, peak_length in peak_info_list])
        dataset_chrom_longest_peak[chrom] = longest_peak_length
    return dataset_chrom_longest_peak
        
def paser_genome_chrom_size():
    genome_size_dict = {}
    with open("/projects/zcsu_research1/npy/cistrom/reference/hg38.chrom.sizes") as fin:
        for line in fin:
            chrom, size = line.strip().split("\t")
            size = int(size)
            genome_size_dict[chrom] = size


    return genome_size_dict

def random_select_in_dataset_peak(chrom, crm_length, chrom_peak_dict):
    chrom_dataset_peaks = chrom_peak_dict[chrom]
    random_peak_index = np.random.randint(len(chrom_dataset_peaks) - 1, size=1)[0]
    
    random_dataset_peak_start, random_dataset_peak_end, random_dataset_peak_length = chrom_dataset_peaks[random_peak_index]
    #random_dataset_peak_length = random_dataset_peak_end - random_dataset_peak_start
    
    # randomly select dataset peak > crm
    
    while random_dataset_peak_length <= crm_length:
        random_peak_index = np.random.randint(len(chrom_dataset_peaks) - 1, size=1)[0]
        random_dataset_peak_start, random_dataset_peak_end, random_dataset_peak_length = chrom_dataset_peaks[random_peak_index]
        
    random_select_start_point = np.random.randint(low=random_dataset_peak_start, high=random_dataset_peak_end-crm_length, size=1)[0]
    random_select_end_point = random_select_start_point + crm_length

    random_peak = "{}\t{}\t{}\n".format(chrom, random_select_start_point, random_select_end_point)
    return random_peak

def random_select_out_dataset_peak(chrom, crm_length,  genome_size_dict):
    size = genome_size_dict[chrom]
    random_select_start_point = np.random.randint(size - crm_length - 1, size=1)[0]
    random_select_end_point = random_select_start_point + crm_length
    random_peak = "{}\t{}\t{}\n".format(chrom, random_select_start_point, random_select_end_point)
    return random_peak


def main():    

    crm_name = sys.argv[1]

    chrom_peak_dict = build_dataset_dict()
    dataset_chrom_longest_peak = get_longest_peak_length_in_dataset_chrom(chrom_peak_dict)
    genome_size_dict = paser_genome_chrom_size()

    random_crm_name = "{}.random_peak".format(crm_name)

    with open(random_crm_name, "w") as fout:
        with open(crm_name) as fin:
            for line in fin:
                line_split = line.strip().split("\t")
                chrom, start, end = line_split[:3]
                crm_length = int(end) - int(start)

                if crm_length > dataset_chrom_longest_peak[chrom]:
                    random_peak = random_select_out_dataset_peak(chrom, crm_length, genome_size_dict)
                else:
                    random_peak = random_select_in_dataset_peak(chrom, crm_length, chrom_peak_dict)

                fout.write(random_peak)


if __name__=="__main__":
    main()

