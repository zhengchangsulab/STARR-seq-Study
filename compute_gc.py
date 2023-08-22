#!/usr/bin/python
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import sys

def compute_gc_content(bed_file):
    genome_file = '/projects/zcsu_research1/npy/cistrom/reference/hg38.fa'  # Replace with the actual path to your reference genome FASTA file
    dict_fasta = SeqIO.to_dict(SeqIO.parse(open(genome_file), 'fasta'))

    with open(f"{bed_file}.gc.bed", "w") as fout:
        with open(bed_file, 'r') as bed:
            for line in bed:
                fields = line.strip().split('\t')
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                
                #Fetch the sequence for the interval
                long_seq_record = dict_fasta[chrom]
                long_seq = long_seq_record.seq
                #alphabet = long_seq.alphabet
                short_seq = str(long_seq)[start-1:end]
                #Compute the GC content for the sequence
                gc_content = gc_fraction(short_seq)
                line_str = line.strip()
                output = f"{line_str}\t{gc_content}\n"
                fout.write(output)

bed_file = sys.argv[1]
compute_gc_content(bed_file)
