import utils as utils
import networkx as nx
import argparse
import json
from Node import Node, NodeEncoder
import pickle
import os
import copy

ALPHABET = ['A', 'C', 'G', 'T']
def parse_args():
    parser = argparse.ArgumentParser(description='Description of your program')
    parser.add_argument('--minimizer_size', type=int, default=5, help='Minimizer size parameter')
    parser.add_argument('--coverage', type=int, default=30, help='Coverage to do, number reads is computed in function')
    parser.add_argument('--minimizer_cov', type=float, default=0.25, help='Proportion of possible minimizers to select')
    parser.add_argument('--minspace_overlap', type=int, default=3, help='Minimum overlap between reads')
    parser.add_argument('--input_sequence_path', type=str, default=None, help='File containing the input sequence file')
    parser.add_argument('--min_cluster_encoding', type=int, default=6, help='Minimum number of minimizers to cluster them together')
    parser.add_argument('--read_length', type=int, default=150, help='Length of each read')
    parser.add_argument('--k', type=int, default=41, help='Length of each kmer')
    parser.add_argument('--error_rate', type=float, default=0.03, help='Error rate for reads')
    parser.add_argument('--neighborhood_size', type=int, default=8, help='Size of the neighborhood to consider')

    args = parser.parse_args()
    return args




    
if __name__ == "__main__":
    args = parse_args()
    minimizer_size = args.minimizer_size
    minspace_overlap = args.minspace_overlap
    minimizer_cov = args.minimizer_cov
    neighborhood_size = args.neighborhood_size
    min_cluster_encoding = args.min_cluster_encoding
    input_sequence_path = args.input_sequence_path
    read_length = args.read_length
    k = args.k
    error_rate = args.error_rate
    minimizer_set_size = int((len(ALPHABET)**minimizer_size) * minimizer_cov)
    selected_minimizers = utils.generate_all_possible_sequences(ALPHABET, minimizer_size)[:minimizer_set_size]




    #GOD bless chat GPT
    res_name = input_sequence_path.split('/')[-1].split('.')[0]
    res_name__with_args = f"{res_name}_cover-{args.coverage}x_minspace_overlap-{minspace_overlap}_kmer-{k}_neighborhood-{neighborhood_size}"

    dir_name = utils.create_new_res_dir("components_" + res_name__with_args)
  


    f_logs = open(f"{dir_name}/logs.txt", 'w')

    G, encoding_index, G_neighborhoods, reads, G_neighborhoods_nodes = None, None, None, None, None



    sequence = None
    print("Loading sequence from file...")
    with open(input_sequence_path, "r") as f_in:
        sequence = ''.join([l.strip() for l in f_in.read().splitlines()][1:])
        seq_length = len(sequence)
        f_in.close()

        assert seq_length > 0, "Sequence length is 0"
        assert all([c in ALPHABET for c in sequence]), "Sequence contains invalid characters"

    print(f"Succesfully loaded sequence of length {len(sequence)}")

    num_reads = int(args.coverage * seq_length / read_length)

    print(f"Generating {num_reads} for coverage {args.coverage}x...")
    reads = utils.generate_reads(sequence, k, num_reads, read_length, ALPHABET, with_errors=True, error_rate=error_rate)

    print(f"Number of reads generated: {len(reads)}")



    print("Building graph...")
    G, encoding_index = utils.create_networkx_graph(reads, selected_minimizers, k, directed=False)




    nodes_with_degree_one_out = [node for node, degree in G.out_degree if degree == 1]
    nodes_with_degree_one_in = [node for node, degree in G.in_degree if degree == 1]
    print("Generating neighborhoods...")



    #Nodes that have no outgoing nodes and their neighbors (end of path)
    G_neighborhoods_nodes_out = [(n, list(nx.generators.ego_graph(nx.reverse_view(G), n, radius=neighborhood_size).nodes)) for n in nodes_with_degree_one_out]
    
    #Nodes that have no incoming nodes and their neighbors (start of path)
    G_neighborhoods_nodes_in = [(n, list(nx.generators.ego_graph(G, n, radius=neighborhood_size).nodes)) for n in nodes_with_degree_one_in]




    
    selected_nodes = set()


    #Taking all nodes from the neighborhoods and adding them the root_node

    #Note: We will have duplicated for nodes in two neighborhoods
    for root_node, nodes in G_neighborhoods_nodes_in + G_neighborhoods_nodes_out:
        for node in nodes:
    
            
            attributes = G.nodes[node]
    
            attributes['root_node'] = root_node

            selected_nodes.add(node)

    print(f"Number of selected nodes: {len(selected_nodes)} out of {len(G.nodes)}")



    res_overlap = utils.find_overlapping_matches(G, selected_nodes, reads, minspace_overlap)

    matches = res_overlap['data']


    print("Building assembly...")

    G_enhanced = copy.deepcopy(G)

   
    for m in matches:
        out_node, in_node = m.replace(' ', '').split('|')
        #print(f"Adding edge {out_node} -> {in_node}")
        G_enhanced.add_edge(out_node, in_node, special = True)
        G_enhanced.add_edge(in_node, out_node, special = True)



    print("Computing strongly connected components for unsimplified and simplified graphs...")

    
    strongly_conne_normal_unsimplified = list(nx.strongly_connected_components(G))
    strongly_conne_special_unsimplified = list(nx.strongly_connected_components(G_enhanced))
    print("Components for unsimplified graph done")


    print("Trimming graph...")
    utils.trim_graph(G_enhanced, k, 7, f=f_logs)
    utils.trim_graph(G, k, 7, f=f_logs)



    print("Computing strongly connected components for simplified graphs...")
    strongly_conne_normal_simplified = list(nx.strongly_connected_components(G))
    strongly_conne_special_simplified = list(nx.strongly_connected_components(G_enhanced))


    with open(f"{dir_name}/res.txt", 'w') as f:
        print(f"Saving results...")


        print("Results for unsimplified graph", file=f)
        print(f"# Components normal Graph Unsimplified: {len(strongly_conne_normal_unsimplified)} components", file=f)
        print(f"# Components special Graph Unsimplified: {len(strongly_conne_special_unsimplified)} components", file=f)

        print(f"Top 10 Size Components normal Graph Unsimplified: {sorted(map(len, strongly_conne_normal_unsimplified), reverse=True)[:10]}", file=f)
        print(f"Top 10 Size Components special Graph Unsimplified: {sorted(map(len, strongly_conne_special_unsimplified), reverse=True)[:10]}", file=f)



        print("Results for simplified graph", file=f)
        print(f"# Components normal Graph Simplified: {len(strongly_conne_normal_simplified)} components", file=f)
        print(f"# Components special Graph Simplified: {len(strongly_conne_special_simplified)} components", file=f)
        print(f"Top 10 Components normal Graph Simplified: {sorted(map(len, strongly_conne_normal_simplified), reverse=True)[:10]}", file=f)
        print(f"Top 10 Components special Graph Simplified: {sorted(map(len, strongly_conne_special_simplified), reverse=True)[:10]}", file=f)



