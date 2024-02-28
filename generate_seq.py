import utils as utils
import networkx as nx
import argparse
import json
from Node import Node, NodeEncoder
import pickle
import os
import copy

ALPHABET = ['G', 'A', 'T', 'C']

def parse_args():
    parser = argparse.ArgumentParser(description='Description of your program')
    parser.add_argument('--seq_length', type=int, default=10_000, help='Length of the sequence')
    args = parser.parse_args()
    return args


    
if __name__ == "__main__":
    args = parse_args()

    file_name = f"random_seq_l-{args.seq_length}.fna"
    with open(f'./input_sequences/{file_name}', 'w') as f:
        f.write(f"> Random Genomic Sequence - Length {args.seq_length}\n")
        sequence = utils.generate_random_nucleotide_string(args.seq_length)
        for i, bp in enumerate(sequence):
            f.write(bp)
            if(i % 70 == 0 and i != 0):
                f.write("\n")
   
