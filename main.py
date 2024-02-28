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
    parser.add_argument('--k', type=int, default=41, help='Size of a kmer')
    parser.add_argument('--save_name', type=str, default = None , help='Save graph to file')
    parser.add_argument('--load_name', type=str, default = None , help='Load graph from file')
    parser.add_argument('--coverage', type=int, default=30, help='Coverage to do, number reads is computed in function')
    parser.add_argument('--read_length', type=int, default=150, help='Length of each read')
    parser.add_argument('--error_rate', type=float, default=0.03, help='Error rate, set to 0 for no errors')
    parser.add_argument('--minimizer_cov', type=float, default=0.25, help='Proportion of possible minimizers to select')
    parser.add_argument('--minspace_overlap', type=int, default=3, help='Minimum overlap between reads')
    parser.add_argument('--input_sequence_path', type=str, default=None, help='File containing the input sequence file')
    parser.add_argument('--min_cluster_encoding', type=int, default=6, help='Minimum number of minimizers to cluster them together')
    parser.add_argument('--bck_size_analysis', type=bool, default=False, help='Whether to perform the bucket size analysis')
    parser.add_argument('--general_statistics', type=bool, default=False, help='Whether to compute statistics about entire graph (VERY EXPENSIVE) ')
    parser.add_argument('--num_threads', type=int, default=1, help='Number of threads to use')
    parser.add_argument('--neighborhood_size', type=int, default=0, help='Size of the neighborhood to consider')
    parser.add_argument('--stub_cleaning_depth', type=int, default=3, help=' Factor of k on which to perform the stub cleaning')
    parser.add_argument('--special_edge_first', action='store_true', help='Whether to add special edges first')
    parser.add_argument('--number_minimizer_sets', type=int, default=1, help='Number of minimizer sets to generate')
    args = parser.parse_args()
    return args




    
if __name__ == "__main__":
    args = parse_args()
    minimizer_size = args.minimizer_size
    minspace_overlap = args.minspace_overlap
    minimizer_cov = args.minimizer_cov
    neighborhood_size = args.neighborhood_size
    k = args.k
    read_length = args.read_length
    min_cluster_encoding = args.min_cluster_encoding
    error_rate = args.error_rate
    input_sequence_path = args.input_sequence_path
    bucket_size_analysis = args.bck_size_analysis

    minimizer_set_size = int((len(ALPHABET)**minimizer_size) * minimizer_cov)

    minimizer_sets = [utils.generate_all_possible_sequences(ALPHABET, minimizer_size)[:minimizer_set_size] for _ in range(args.number_minimizer_sets)]



    #GOD bless chat GPT
    res_name = input_sequence_path.split('/')[-1].split('.')[0]
    res_name__with_args = f"{res_name}_cover-{args.coverage}x_minspace_overlap-{minspace_overlap}_kmer-{k}_neighborhood-{neighborhood_size}_cleaning-depth-{args.stub_cleaning_depth}_special_edge_first-{args.special_edge_first}_er-{error_rate}_nminimizer_sets-{args.number_minimizer_sets}"

    dir_name = utils.create_new_res_dir(res_name__with_args)
  


    f_logs = open(f"{dir_name}/logs.txt", 'w')


    # Unfortunately loading graphs from memory (else part) creates undeterministic behavior

    print("No graph provided, generating graph...")

    sequence = None
    seq_length = None
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
    G, encoding_indexes = utils.create_networkx_graph(reads, minimizer_sets, k)


    G_undirected = G.to_undirected()
    connected_components = list(nx.connected_components(G_undirected))


    #TODO: One could do this only for the nodes where between which we might add an edge ie. nodes_leafs union nodes_source
    for n in G.nodes:
        G.nodes[n]['component'] = utils.find_component(n, connected_components)




    #Set to 1 if graph is directed
    nodes_leaf = [node for node, degree in G.out_degree if degree == 0]
    nodes_source = [node for node, degree in G.in_degree if degree == 0]

    #Sanity check: There should not be any singleton nodes
    assert len(set(nodes_leaf).intersection(set(nodes_source))) == 0, "There should not be any singleton nodes"


    print(f"Number of connected components: {len(connected_components)}")
    print(f"Number of nodes in graph: {len(G.nodes)}")

    

    #Now these are only added between different components
    res_overlap = utils.find_overlapping_matches(G, nodes_leaf, nodes_source, reads, minspace_overlap, connected_components)

    matches = res_overlap['data']


    print("Building assembly...")

    G_enhanced = copy.deepcopy(G)


   #Special edges added (in_node, out_node), weight)
    special_edges = []
    for m in matches:
        out_node, in_node = m.replace(' ', '').split('|')
        minimizer_set_idx = int(matches[m]['minimizer_set_idx'])

        print(f"Adding edge between {out_node} --> {in_node}...", file=f_logs)
        print(f"Encoding: {G.nodes[out_node]['minimizer_encodings']} --> {G.nodes[in_node]['minimizer_encodings']}", file=f_logs)
        print(f"Number predecessors of out_node: {len(list(G.predecessors(out_node)))}", file=f_logs)
        print(f"Number successors of out_node: {len(list(G.successors(out_node)))}", file=f_logs)
        print(f"Number predecessors of in_node: {len(list(G.predecessors(in_node)))}", file=f_logs)
        print(f"Number successors of in_node: {len(list(G.successors(in_node)))}", file=f_logs)
        print()

        edge_weight = matches[m]['edge_weight']


        #We add the edge to the special graph only once even if match from several minimizer sets
        if(not G_enhanced.has_edge(out_node, in_node)):
            G_enhanced.add_edge(out_node, in_node, special = True, minimizer_set_idx=minimizer_set_idx, weight=edge_weight)
        special_edges.append(((out_node, in_node), edge_weight))


    #To control the nodes from which the assembly will start
    special_edges = special_edges if args.special_edge_first else {}

    G_enhanced_components = list(nx.connected_components(G_enhanced.to_undirected()))

    assembly_normal_unsi, paths_normal_unsi = utils.build_assembly(G, set(), minimizer_sets, basic=True)
    assembly_special_unsi, paths_special_unsi = utils.build_assembly(G_enhanced, special_edges, minimizer_sets, basic=False, g_original=G)


    G_enhanced_unsi_size = len(G_enhanced.nodes)
    
    print(f"Trimming graph special, size before simplification: {len(G_enhanced.nodes)}")
    print(f"Trimming special graph.....")
    utils.trim_graph(G_enhanced, int(args.stub_cleaning_depth * k), 7, f=f_logs)
    G_enhanced_simplified_components = list(nx.connected_components(G_enhanced.to_undirected()))
    G_enhanced_si_size = len(G_enhanced.nodes)

    print()

    G_unsi_size = len(G.nodes)
    utils.trim_graph(G, int(args.stub_cleaning_depth * k), 7, f=f_logs)
    G_si_size = len(G.nodes)

    G_simplified_components = list(nx.connected_components(G.to_undirected()))
    


    assembly_normal_simplified, paths_normal_simplified = utils.build_assembly(G, set(), minimizer_sets, basic=True)
    print(f"Assembly normal simplified, size: {len(assembly_normal_simplified)}")
    assembly_special_simplified, paths_special_simplified = utils.build_assembly(G_enhanced, special_edges, minimizer_sets, basic=False, g_original=G)
    print(f"Assembly special simplified, size: {len(assembly_special_simplified)}")

    #TODO: 
    #Should we add nodes in same bucket - i.e. take all keys bigger than some param in index and add edges between all those nodes



    len_n50_normal_unsi, n50_idx_normal_unsi = utils.compute_N50(assembly_normal_unsi, seq_length)
    len_n50_special_unsi, n50_idx_special_unsi = utils.compute_N50(assembly_special_unsi, seq_length)

    len_n50_normal_simplified, n50_idx_normal_simplified = utils.compute_N50(assembly_normal_simplified, seq_length)
    len_n50_special_simplified, n50_idx_special_simplified = utils.compute_N50(assembly_special_simplified, seq_length)





    print(f"Saving statistics...")
    with open(f"{dir_name}/general_stats.txt", 'w') as f:

        #Computing the metrics

        number_valid_matches_found = res_overlap['meta_data']['number_valid_matches_found']
        number_invalid_matches_found = res_overlap['meta_data']['number_invalid_matches_found']
        additional_components_connected = res_overlap['meta_data']['additional_components_connected']
        # number_components = res_overlap['meta_data']['number_connected_components']
        # component_connections = res_overlap['meta_data']['component_connections']
        
        accuracy = None
        if (1.0 * number_valid_matches_found + number_invalid_matches_found) != 0:
            accuracy = number_valid_matches_found /  (1.0 * number_valid_matches_found + number_invalid_matches_found)

        print(f"-----Graph statistics-----", file=f)
        print(f"Number of nodes in graph: {len(G.nodes)}", file=f)
        print(f"Number of selected leaves: {len(nodes_leaf)}", file=f)
        print(f"Number of selected sources: {len(nodes_source)}", file=f)
        print(f"Number of selected nodes: {len(nodes_source) + len(nodes_leaf)}", file=f)


        print(f"Number of valid matches found: {number_valid_matches_found}, number of invalid matches found: {number_invalid_matches_found}", file=f)
        # print(f"Number of additional components connected: {additional_components_connected} among {number_components}", file=f)
        
        if accuracy is not None:
            print(f"Accuracy : {100 * round(accuracy, 3)}%", file=f)

        print(f"Number of nodes in graph normal (unsimplified): {G_unsi_size}", file=f)
        print(f"Number of nodes in graph special (unsimplified): {G_enhanced_unsi_size}", file=f)

        print(f"Number of nodes in graph normal (simplified): {G_si_size}", file=f)
        print(f"Number of nodes in graph special (simplified): {G_enhanced_si_size}", file=f)
    
        print(f"Number of nodes removed in graph normal (unsimplified): {G_unsi_size - G_si_size}", file=f)
        print(f"Number of nodes removed in graph special (unsimplified): {G_enhanced_unsi_size - G_enhanced_si_size}", file=f)
        print("\n", file=f)


        print(f"-----Assembly statistics-----", file=f)
        print(f"-----Unsimplified-----", file=f)
        print(f"Number of contigs in assembly (unsimplified): {len(assembly_normal_unsi)}", file=f)
        print(f"Number of contigs in assembly with additional edges (unsimplified): {len(assembly_special_unsi)}", file=f)

        print(f"N50 for normal assembly (unsimplified): {len_n50_normal_unsi} at index {n50_idx_normal_unsi}", file=f)
        print(f"N50 for special assembly (unsimplified): {len_n50_special_unsi} at index {n50_idx_special_unsi}", file=f)
        print("\n", file=f)

        print(f"-----Simplified-----", file=f)

        print(f"Number of contigs in assembly (simplified): {len(assembly_normal_simplified)}", file=f)
        print(f"Number of contigs in assembly with additional edges (simplified): {len(assembly_special_simplified)}", file=f)
       
        print(f"N50 for normal assembly (simplified): {len_n50_normal_simplified} at index {n50_idx_normal_simplified}", file=f)
        print(f"N50 for special assembly (simplified): {len_n50_special_simplified} at index {n50_idx_special_simplified}", file=f)
        print("\n", file=f)




        print(f"-----Paths statistics-----", file=f)
        print(f"-----Unsimplified-----", file=f)
        paths_len_normal_unsi= sorted([len(p) for p in paths_normal_unsi], reverse=True)
        paths_len_special_unsi = sorted([len(p) for p in paths_special_unsi], reverse=True)
        paths_len_normal_simplified = sorted([len(p) for p in paths_normal_simplified], reverse=True)
        paths_len_special_simplified = sorted([len(p) for p in paths_special_simplified], reverse=True)

        print(f"Number of paths in normal assembly (unsimplified): {len(paths_len_normal_unsi)}", file=f)
        print(f"Number of paths in special assembly (unsimplified): {len(paths_len_special_unsi)}", file=f)

        print(f"Top 50 Longest path in normal assembly (unsimplified): {paths_len_normal_unsi[:50]}", file=f)
        print(f"Top 50 Longest path in special assembly (unsimplified): {paths_len_special_unsi[:50]}", file=f)
        print("\n",file=f)

        print(f"-----Simplified-----", file=f)

        print(f"Number of paths in normal assembly (simplified): {len(paths_len_normal_simplified)}", file=f)
        print(f"Number of paths in special assembly (simplified): {len(paths_len_special_simplified)}", file=f)


        print(f"Top 50 Longest path in normal assembly (simplified): {paths_len_normal_simplified[:50]}", file=f)
        print(f"Top 50 Longest path in special assembly (simplified): {paths_len_special_simplified[:50]}", file=f)




    print(f"Saving connected components sizes...")
    with open(f"{dir_name}/connected_components_sizes.txt", 'w') as f:

        idx_with_size = sorted([(idx, len(c)) for idx, c in enumerate(connected_components)], key=lambda x: x[1], reverse=True)

        # Print all component sizes
        print("--------Original graph--------", file=f)
        for idx, size in idx_with_size:
            print(f"Component {idx} size: {size}", file=f)

        
        idx_with_size = sorted([(idx, len(c)) for idx, c in enumerate(G_enhanced_components)], key=lambda x: x[1], reverse=True)
        print("--------Special graph--------", file=f)
        for idx, size in idx_with_size:
            print(f"Component {idx} size: {size}", file=f)


        print("--------Original graph simplified--------", file=f)
        idx_with_size = sorted([(idx, len(c)) for idx, c in enumerate(G_simplified_components)], key=lambda x: x[1], reverse=True)
        for idx, size in idx_with_size:
            print(f"Component {idx} size: {size}", file=f)

        print("--------Special graph simplified--------", file=f)
        idx_with_size = sorted([(idx, len(c)) for idx, c in enumerate(G_enhanced_simplified_components)], key=lambda x: x[1], reverse=True)
        for idx, size in idx_with_size:
            print(f"Component {idx} size: {size}", file=f)


    


    print(f"Saving assemblies..")

    with open(f"{dir_name}/assembly_normal_unsimplified.fasta", 'w') as f:
        utils.write_assembly(assembly_normal_unsi, f)


    with open(f"{dir_name}/assembly_special_unsimplified.fasta", 'w') as f:
        utils.write_assembly(assembly_special_unsi, f)

    
    with open(f"{dir_name}/assembly_normal_simplified.fasta", 'w') as f:
        utils.write_assembly(assembly_normal_simplified, f)
    
    with open(f"{dir_name}/assembly_special_simplified.fasta", 'w') as f:
        utils.write_assembly(assembly_special_simplified, f)

    
    print(f"Saving additional edges Json details..")
    with open(f"{dir_name}/matches_details.json", 'w') as f:
        json.dump(res_overlap, f, indent=4)

    
    print(f"Saving readset..")
    with open(f"{dir_name}/reads.txt", 'w') as f:
        utils.write_assembly(reads.keys(), f)


 
"""
    if(args.general_statistics):
        print("Computing general statistics...")

        res_all = utils.find_allmatches(G, selected_nodes, reads)

        number_valid_matches = res_all['meta_data']['number_valid_matches']
        number_invalid_matches = res_all['meta_data']['number_invalid_matches']
        total_number_pairs = res_all['meta_data']['total_number_pairs']

        print(f"Number of valid matches: {number_valid_matches}, number of invalid matches: {number_invalid_matches}")


        true_pos_rate = number_valid_matches_found/ (1.0 * number_valid_matches)
        false_pos_rate = number_invalid_matches_found/ (1.0 * number_invalid_matches)

        false_positive = number_invalid_matches_found


        print(f"True positive rate : {100 * round(true_pos_rate, 3)}%", file=f)
        print(f"False negative rate : {100 * round(false_pos_rate, 3)}%", file=f)
"""



        
