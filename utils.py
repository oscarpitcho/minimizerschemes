from itertools import product
import random

import networkx as nx
from typing import List, Set, Tuple, Dict
from Node import Node
from networkx.classes.function import all_neighbors



def intervals_overlap(interval1, interval2):
    r1 = set(range(interval1[0], interval1[1] + 1))
    r2 = set(range(interval2[0], interval2[1] + 1))
    return len(r1.intersection(r2))  


#Let's just store a sequence of nodes first, we can then use those to generate the graph
#We will need to label to which path(s) belongs each edge

def build_paths(reads, k_mersize) -> Dict[str, List[str]]:
    paths = {}
    for read in reads:
        path = []
        if(len(read) >= k_mersize and k_mersize > 0):
            for i in range(0, len(read) - k_mersize + 1):
                path.append(read[i:i+k_mersize])
        paths[read] = path
    return paths



def build_adj_and_coverage(reads, k_mersize) -> Dict[str, List[str]]:
    adjacencies = {}
    read_coverage = {}
    for read in reads:
        for i in range(0, len(read) - k_mersize - 1):
                n = read[i:i+k_mersize]
                n_succ = read[i+1:i+k_mersize+1]
                if n not in adjacencies:
                    adjacencies[n] = {}
                if n_succ not in adjacencies:
                    adjacencies[n_succ] = {}

                if n_succ not in adjacencies[n]:
                    adjacencies[n][n_succ] = 0

                if n not in read_coverage:
                    read_coverage[n] = set()
                if n_succ not in read_coverage:
                    read_coverage[n_succ] = set()
                
                adjacencies[n][n_succ] += 1


                read_coverage[n].add(read)
                read_coverage[n_succ].add(read)

    return adjacencies, read_coverage


def find_minimizers_in_string(s, minimizers)-> List[str]:
    minimizer_lengths = [len(m) for m in minimizers]
    smallest_min = min(minimizer_lengths) if len(minimizer_lengths) > 0 else 0
    encoded_s = []
    i = 0 

    while i < len(s) - smallest_min + 1:
        found = False
        for j in range(len(minimizers)):
            if s[i:i + len(minimizers[j])] == minimizers[j]:
                encoded_s.append(f'M_{j}')
                i += len(minimizers[j])
                found = True
                break

        if not found:
            i += 1
    return encoded_s



def create_networkx_graph(reads, minimizer_sets: List[List[str]], k: int, directed=True):
    """Creates a networkx dBJ from the set of reads. Each node includes the following attributes:
        - corres_reads: the set of reads that walked through this node
        - minimizer_encoding: the encoded node using the provided minimizer set
        
        
        returns:
            - The networkx graph
            - An index of of the graph on the attribute of enconding for fast querying encodings"""
    


    #TODO: Replace reads dictionary with list and index --> store indices in adjacency list
    print("Building adjacency list...")
    adjacencies, read_coverage = build_adj_and_coverage(reads, k)
    print(f"adjacency list length: {len(adjacencies)}")

    #We build one index for each set of minimizers
    indexes = [{} for _ in range(len(minimizer_sets))]


    G = nx.DiGraph()
    j = 0
    
    for node, neighbors in adjacencies.items():
        j += 1


        if j % 10_000 == 0:
            print(f"Bulding Graph, processed {j}/{len(adjacencies)} nodes")


        minimizer_encodings = [find_minimizers_in_string(node, minimizer_set) for minimizer_set in minimizer_sets]
        encoded_strings = [''.join(encoding) for encoding in minimizer_encodings]
        if(not G.has_node(node)):
            G.add_node(node)

        #We add the metadata of the node

        attributes = G.nodes[node]

        attributes['corres_reads'] = read_coverage[node]
        attributes['minimizer_encodings'] = encoded_strings

        #print(G.nodes[node])




        #Add the node to the index, one index for each set of minimizers
        for i in range(len(minimizer_sets)):
            if(encoded_strings[i] not in indexes[i]):
                indexes[i][encoded_strings[i]] = set(node)
            else:
                indexes[i][encoded_strings[i]].add(node)
                
 


        #We then iterate over the neighbors and add them to the graph
        for neigh in neighbors:

            #Meta data added later
            if(not G.has_node(neigh)):
                G.add_node(neigh)

            
            
            G.add_edge(node, neigh, weight=adjacencies[node][neigh])
           
           
            if(not directed):
                G.add_edge(neigh, node, weight=adjacencies[node][neigh])


    return G, indexes


def generate_all_possible_sequences(alphabet, length):
    """returns possible sequences in shuffled order"""
    res = [''.join(x) for x in product(''.join(alphabet), repeat=length)]
    random.shuffle(res)
    return res


def generate_random_nucleotide_string(length):
    nucleotides = ['A', 'T', 'C', 'G']
    return ''.join(random.choice(nucleotides) for _ in range(length))




def generate_reads(nucleotide_string, k, num_reads, read_length,alphabet, with_errors=False, error_rate=0.03):
    """Generates reads from a nucleotide string:
    - Assumes that the reads are of the same length.
    - There are no two identical duplicate reads
    - The read start points are distributed uniformly across the nucleotide string
    
    returns:
    - reads: a dictionary of reads and their start and end positions in the nucleotide string"""
    reads = {}  # Using a set to ensure uniqueness -
    max_length = len(nucleotide_string)
    
    read_conflicts = 0
    while len(reads) < num_reads:
        start_pos = random.randint(0, max_length-1)
        end_pos = start_pos + read_length
    
        read = nucleotide_string[start_pos:start_pos + read_length]
        if(len(read) < read_length):
            continue
        
        if(with_errors):
            read = list(read)
            for i in range(len(read)):
                if random.random() < error_rate:
                    read[i] = random.choice(alphabet)
            read = ''.join(read)

        if read in reads:
            read_conflicts += 1


        #We want start and end indices to be included
        reads[read] = (start_pos,end_pos - 1)
    
    print(f"Number of read conflicts: {read_conflicts}")

    return reads


#Takes two strings and returns the size of overlap (returns the max overlap between the two sides)
#Postitive if end of seq1 overlaps with start of seq2, negative if reverse
def suffix_overlap(seq1, seq2):
    overlap = ""
    for i in reversed(range(min(len(seq1),len(seq2)))):
        if (seq1[-i - 1:] == seq2[:i+1]):
            overlap = seq1[-i - 1:]
        else:
            return overlap
        
    return overlap


#TODO: Node encapsulating class that should be hashable by mapping it to a string
def code_to_string(list):
    return ''.join(list)


"""Takes a string and set of minimizers and returns the string as a list of minimizers it contains"""
def minimizer_string_to_list(string):
    import re
    pattern = r'M_\d+'  # '\d+' matches one or more digits
    result = re.findall(pattern, string)
    return result

def compute_suffix_buckets(G, nodes_leaf, minspace_overlap):
    """
    Takes a set of nodes and returns two dictionaries containing prefixes/suffixes and the nodes starting/ending with them
    """


    overlap_buckets_per_set = [{} for _ in range(len(G.nodes[nodes_leaf[0]]['minimizer_encodings']))]
    print(f"Number of overlap buckets maps: {len(overlap_buckets_per_set)}")
    print(f"Computing overlap buckets...")


    for nl in nodes_leaf:

        encodings = G.nodes[nl]['minimizer_encodings']



        for i, encoding_str in enumerate(encodings):
            encoding = Node.encoded_s_to_list(encoding_str)

            if(len(encoding) >= minspace_overlap):
                #Prefix and suffix should be strings for hashing

                suffix = ''.join(encoding[-minspace_overlap:])

                if suffix not in overlap_buckets_per_set[i]:
                    overlap_buckets_per_set[i][suffix] = set()


                overlap_buckets_per_set[i][suffix].add(nl)


    
    print(f"Finished computing overlap buckets")
    print(f"Number of overlap bucket sets {len(overlap_buckets_per_set)}")
    return overlap_buckets_per_set



def find_overlapping_matches(G, nodes_leaf, nodes_source, reads, minspace_overlap, connected_components):
    def component_connection_exists(c1, c2, component_connections):
        return (f"c{c1} | c{c2}" in component_connections) or f"c{c1} | c{c2}" in component_connections
        


    """
    returns a dictionary containing
    {
        meta_data : {\n
        'number_valid_matches_found' : number_valid_matches_found,\n
        'number_invalid_matches_found' : number_invalid_matches_found
        },\n
        'data' : data
    
    }

    data follows the following template.
    Each entry in the dictionary is equivalent to pairs of nodes that overlap in minimizer space but with different root nodes

    (left_node, right_node) : {\n
        'same_readset': reads1 == reads2,\n
        'overlaps': overlaps,\n
        'root_nodes': (rn1, rn2)
    }
    
    """

    print(f"Computing overlap buckets...")
    suffix_buckets_per_set = compute_suffix_buckets(G, nodes_leaf, minspace_overlap)


    j = 0
    number_valid_matches_found = 0
    number_invalid_matches_found = 0

    data  = {}
    
    #components = list(nx.strongly_connected_components(G))

    component_connections = set()
    additional_components_connected = 0



    #We create connection LEAF_NODE ---> SOURCE_NODE
    #We iterate over each source node to find leaf nodes that share the same suffix

    for i, suffix_buckets in enumerate(suffix_buckets_per_set):
        for sn in nodes_source:
            sn_attr = G.nodes[sn]
            sn_c = sn_attr['component']

            sn_encoding = Node.encoded_s_to_list(sn_attr['minimizer_encodings'][i])

            j += 1
            sn_prefix = ''.join(sn_encoding[:minspace_overlap])
            if j % 10_000 == 0:
                print(f"Finding overlap matches - Processed {j}/{len(nodes_leaf)} leaf nodes")

            #Checking the nodes that end with the start of N1    
            if sn_prefix in suffix_buckets:
                for ln in suffix_buckets[sn_prefix]:
                    ln_attr = G.nodes[ln]

                    ln_c = ln_attr['component']


                    #We did BFS to get the neighborhoods,
                    #so we can assume that if two nodes are in the same neighborhood they are connected (not interesting)
                    #We also check that they are from different components
                    if ln_c != sn_c:

                        reads_ln = ln_attr['corres_reads']
                        reads_sn = sn_attr['corres_reads']
                        edge_weight = len(connected_components[ln_c]) + len(connected_components[sn_c])

                        #Note: Iterating over different data structures, its a mess
                        
                        overlaps = []
                        for r1 in reads_ln:
                            for r2 in reads_sn:
                                overlaps.append(intervals_overlap(reads[r1], reads[r2]))

                        if(reads_ln == reads_sn):
                            overlaps.append(-1)

                        if any(overlaps) == 0:
                            number_invalid_matches_found += 1
                        else:
                            number_valid_matches_found += 1

                        data[f"{ln} | {sn}"] = {
                        'same_readset': reads_ln == reads_sn,
                        'overlaps': overlaps,
                        'components' : (ln_c, sn_c),
                        'edge_weight' : edge_weight,
                        'encodings' : (sn_attr['minimizer_encodings'][i], ln_attr['minimizer_encodings'][i]),
                        'minimizer_set_idx' : i
                        }


                        #Checking if this connection connects two components


        meta_data  = {
            'number_valid_matches_found' : number_valid_matches_found,
            'number_invalid_matches_found' : number_invalid_matches_found,
            #'number_connected_components' : len(components),
            'additional_components_connected' : additional_components_connected,
            'component_connections' : list(component_connections)

        }

    return {'data' : data, 'meta_data' : meta_data}
                    
                    

def find_all_matches(G, selected_nodes, reads): 
    """
    returns:
    {'meta_data' : {
        'total_number_pairs' : total_number_pairs,
        'number_valid_matches' : number_valid_matches,
        'number_invalid_matches' : number_invalid_matches
    }, 'data' : data}
    
    """


     #Getting the total matches that should be found in the graph


    data = {} 



    total_number_pairs  = 0
    number_valid_matches = 0
    number_invalid_matches = 0
    print(f"Processing all matches...")

    i = 0
    for n1 in selected_nodes:
        i += 1
        if i % 100 == 0:
            print(f"Processed {i}/{len(selected_nodes)} nodes for all matches")
        n1_attributes = G.nodes[n1]
        for n2 in selected_nodes:
            n2_attributes = G.nodes[n2]
            if(n1 != n2 and n1_attributes['root_node'] != n2_attributes['root_node']):
                reads1 = n1_attributes['corres_reads']
                reads2 = n2_attributes['corres_reads']

                overlaps = []
                for r1 in reads1:
                    for r2 in reads2:
                        if r1 != r2:
                            overlaps.append(intervals_overlap(reads[r1], reads[r2]))


                #TO be consistent with otherfunction, we count a match only if all reads that walked through the nodes overlapped
                contiguous = all(overlaps) > 0
                
                data[(n1, n2)] = {
                    'overlaps': overlaps,
                    'contiguous' : contiguous
                }

                if(contiguous):
                    number_valid_matches += 1
                else:
                    number_invalid_matches += 1

                total_number_pairs += 1


    
    return {'meta_data' : {
        'total_number_pairs' : total_number_pairs,
        'number_valid_matches' : number_valid_matches,
        'number_invalid_matches' : number_invalid_matches
    }, 'data' : data}
    


"""Returns the index of the component that contains the node"""
def find_component(node, components):
    for i in range(len(components)):
        if node in components[i]:
            return i
    return -1


def build_assembly(graph: nx.DiGraph, special_edges, minimizer_sets: List[List[str]], basic=True, g_original=None):
    """Builds an assembly from a graph. The assembly is a list of nodes that are connected in the graph
    If basic is set to false then we will check for added connections and adjust the assembly accordingly"""



    """
    TODO: Two behaviors to define for the assembly:
    - Case when node is added to the stack by two different nodes before being walked.
    Right now it can occur in two different paths under these circumstances.

    - Should we flag nodes as being walked or edges? I.e. reuse a node as long as we don't reuse the edges.
    
    """

    #We only get the special edges that are present in the graph, might have been removed by the trimming
    special_edges_present = [(edge, w) for edge, w in special_edges if graph.has_edge(edge[0], edge[1])]

    paths = build_walks(graph, special_edges_present, basic=basic)
    #z(F"paths: {paths}")

    assembly = transform_walks_to_contigs(graph, paths, minimizer_sets, basic=basic)

    return assembly, paths

def build_walks(graph: nx.DiGraph, special_edges, basic):


    paths = []
    walked_edges = set()
    stack = []

    #We first launch walks from all special edges

    if not basic:
        decreasing_weight_special_edges = sorted(special_edges, key=lambda x: x[1], reverse=True)

        for ((out_node, in_node), _) in decreasing_weight_special_edges:
            #Checking that the edge has not been walked by another special walk
            if not (out_node, in_node) in walked_edges:
              
                downstream_walk = walk_graph_downstream(graph, in_node, stack, walked_edges)
                upstream_walk = walk_graph_upstream(graph, out_node, walked_edges)

                #print(f"Length of upstream walk: {len(upstream_walk)}")
                #print(f"Length of downstream walk: {len(downstream_walk)}")

                """
                if(len(upstream_walk) == 1):
                    print(f"Upstream walk of length 1, out node: {out_node}, in node: {in_node}")
                    neigh = list(graph.predecessors(out_node))
                    print(f"Number of predecessors of out node: {len(neigh)}")
                    print(f"Number of successors of out node: {len(list(graph.successors(out_node)))}")
                    print(f"Number of neighbors (succ/predec) of in node: {len(list(all_neighbors(graph, in_node)))}")

                
                if(len(downstream_walk) == 1):
                    print(f"Downstream walk of length 1, out node: {out_node}, in node: {in_node}")
                    neigh = list(graph.successors(in_node))
                    print(f"Number of predecessors of in node: {len(list(graph.predecessors(in_node)))}")
                    print(f"Number of successors of in node: {len(neigh)}")
                    print(f"Number of neighbors (succ/predec) of in node: {len(list(all_neighbors(graph, in_node)))}")
                """

                #print(f"Upsream walk: {upstream_walk}")

                walk = list(reversed(upstream_walk)) + downstream_walk

                #print(f"Special walk length: {len(walk)}")
                walked_edges.add((out_node, in_node))
                paths.append(walk)
            else:
                print(f"Edge {out_node} -> {in_node} has already been walked by another special walk")
        
    #We then add the nodes with indegree 0 to the end of the stack
    indegree_0_nodes = [node for node, degree in graph.in_degree if degree == 0]
    indegree_0_nodes = sorted(indegree_0_nodes)
    
    stack = stack + indegree_0_nodes


    curr_path = []
    while len(stack) > 0:
        curr_node = stack.pop()
        curr_path.append(curr_node)

        outgoing_edges = set(graph.out_edges(curr_node))

        unwalked_outgoing_edges = outgoing_edges.difference(walked_edges)

        #Sorting the successors for determinism
        unwalked_outgoing_edges = sorted(unwalked_outgoing_edges, key=lambda x: x[1])

        #Last node in the path
        if(len(unwalked_outgoing_edges) == 0):
            paths.append(curr_path)
            curr_path = []
        else:

            #We add the neighbors of the unwalked outgoing edges to the stack
            for edge in unwalked_outgoing_edges:
                stack.append(edge[1])
                #We mark the edge as walked
                walked_edges.add(edge)
    
    return paths


        

def walk_graph_upstream(graph: nx.DiGraph, starting_node: str, walked_edges: Set):
    """
    Emulates an upstream dfs walk on the graph starting from the starting node.
    Modifies the walked edges set

    returns:
    - A list of nodes walked in upstream order
    """

    walk = []

    curr_node = starting_node

    #print(f"---UPSTREAM WALK STARTING FROM NODE {starting_node}---")

    
    
    walk.append(curr_node)
    #There is an incoming edge that has not been walked
    while not len(set(graph.in_edges(curr_node)).difference(walked_edges)) == 0 :

        #print(f"Walked edges: {walked_edges}")


        #We select the next node to walk at random 
        predecessors = set(graph.predecessors(curr_node))


        #We remove the edges that have already been walked
        predecessors = {p for p in predecessors if (p, curr_node) not in walked_edges}




        #We select the next node to walk using lexigraphic order
        next_node = sorted(predecessors)[0]
        #print(f"Node has {len(predecessors)} unwalked predecessors, selected node: {next_node}")

        #Walking upstream so the edge is reversed
        walked_edges.add((next_node, curr_node))

        walk.append(next_node)

        curr_node = next_node

    predecessors = set(graph.predecessors(curr_node))
    accessible_predecessors = {p for p in predecessors if (p, curr_node) not in walked_edges}

    #print(f"last node number of predecessors: {len(predecessors)}")
    #print(f"last node number of accessible predecessors: {len(accessible_predecessors)}")
    #print(f"---UPSTREAM WALK ENDING FROM NODE {starting_node}---")

    return walk

def walk_graph_downstream(graph: nx.DiGraph, starting_node: str, stack: List, walked_edges: Set):
    walk = []

    curr_node = starting_node
    walk.append(curr_node)



    while not len(set(graph.out_edges(curr_node)).difference(walked_edges)) == 0:

        #We select the next node to walk at random 
        successors = set(graph.successors(curr_node))

        #We remove the edges that have already been walked
        successors = {s for s in successors if (curr_node, s) not in walked_edges}

        #We select the first node in alphabetic order- this is to ensure determinism
        next_node = sorted(successors)[0]
        walked_edges.add((curr_node, next_node))


        #We add the neighbors of the unwalked outgoing edges to the stack
        for edge in graph.out_edges(curr_node):
            stack.append(edge[1])
            #We mark the edge as walked
            walked_edges.add(edge)

        
        walk.append(next_node)
        curr_node = next_node


    return walk


def transform_walks_to_contigs(G, paths, minimizer_sets: List[List[str]], basic=True):
    """Walks the paths and returns a list of contigs


    Args:
        G (nx.DiGraph): The graph over which the paths were computed
        paths (List[List[str]]): The paths, each is a sequence of nodes
        minimizer_sets (List(List[str])): 
        basic (bool, optional): If true, we assume all edges in the paths are normal. If false, we check for special edges. 
    
    """

    assembly = []

    path_is_special = False

    for path in paths:
        
        curr_contig  = []
        for n_idx in range(len(path)):
            curr_node = path[n_idx]


            if(n_idx == 0):
                curr_contig = curr_contig + list(curr_node)

            elif(n_idx > 0):

                #We assume all paths are continueous
                if(basic):
                    # Append the last character
                    curr_contig.append(curr_node[-1])

                else:
                    prev_node = path[n_idx - 1]

                    #Check this syntax
                    #print(f"Prev node: {prev_node}, curr node: {curr_node}, {prev_node in graph.nodes}, {curr_node in graph.nodes}")

                    edges_attr = G.edges[prev_node, curr_node]

                    if(edges_attr.get('special', False)):

                        minimizer_set_idx = edges_attr['minimizer_set_idx']
                        minimizer_set = minimizer_sets[minimizer_set_idx]


                        prev_node_encoding = G.nodes[prev_node]['minimizer_encodings'][minimizer_set_idx]

                      
           
                        final_minimizer = minimizer_string_to_list(prev_node_encoding)[-1]
                        
                        #M_XXX
                        #print(f"Final minimizer: {final_minimizer}")
                        minimizer_idx = int(final_minimizer.split('_')[1])

                        minimizer_base_sequence = minimizer_set[minimizer_idx]

                        #print(minimizer_base_sequence)


                        #This piece of code needs to be tested
                        import re
                        idx_prev_node_last_occ = [m.start() for m in re.finditer(minimizer_base_sequence, prev_node)][-1]
                        idx_curr_node_first_occ = [m.start() for m in re.finditer(minimizer_base_sequence, curr_node)][0]

                        prev_node_trailing_chars = len(prev_node[idx_prev_node_last_occ + len(minimizer_base_sequence):])

                        curr_node_trailing_chars_idx = min(idx_curr_node_first_occ + len(minimizer_base_sequence) + prev_node_trailing_chars, len(curr_node))
                        
                        trailing_bases = curr_node[curr_node_trailing_chars_idx:]


                        curr_contig = curr_contig + list(trailing_bases)


                    # Normal edge, we append the last character
                    else:
                        curr_contig.append(curr_node[-1])
                    

        curr_contig = ''.join(curr_contig)

        assembly.append(curr_contig)



    return assembly

        
def compute_N50(contigs, seq_length):
    print(f"Computing N50 for {len(contigs)} contigs")
    contigs = sorted(contigs, key=len, reverse=True)
    total_length = 0
    n50_contig_idx = 0

    for contig in contigs:
        print(f"Current contig length: {len(contig)}")
        n50_contig_idx += 1
        total_length += len(contig)
        if(total_length >= seq_length/2):
            return len(contig), n50_contig_idx
        
    return -1, -1


def write_assembly(assembly, out_f):
    assembly_sorted = sorted(assembly, key=len, reverse=True)
    for i, contig in enumerate(assembly_sorted):
        out_f.write(f">CONTIG_{i}_length_{len(contig)}\n")
        for i in range(0, len(contig), 60):
            out_f.write(contig[i:i+60] + "\n")


import sys 
"IMPORTANT, THIS FUNCTION MODIFIES THE INPUT GRAPH"     
def trim_graph(graph: nx.DiGraph, stub_length: int, iterations: int, f=sys.stdout):
    def investigate_stub(node, graph_to_main, i):

        curr_node = node
        nodes_walked = set(curr_node)


        for j in range(stub_length):
            #print(f"Current node: {curr_node} with stub length {j} in iteration {i}")
            neighbors = list(graph_to_main.neighbors(curr_node))

            #If the node has several neighbors, this implies we stop (bubble like or incoming stub)
            if(len(neighbors) > 1):
                #print(f"Iteration {i} - Found several incoming nodes when walking up a stub, # nodes walked: {len(nodes_walked)}, no removal", file=f)
                break


            #We were walking a small isolated path, remove all nodes
            if(len(neighbors) == 0):
                nodes_walked.add(curr_node)
                graph.remove_nodes_from(nodes_walked)
                print(f"Iteration {i} - Isolated path - Removed {len(nodes_walked)} nodes", file=f)
                break


            #If we have already walked this node
            #It can happen but the starting node should not be in the set
            if(curr_node in nodes_walked):
                print(f"Iteration {i} - Found a cycle - Not removing nodes", file=f)
                break

            
            #Note we are checking in the orginal graph
            #This node is the start of a fork
            #We stop and remove all nodes walked
            if(graph.out_degree[curr_node] > 1):

                #We do not remove the current node as it is the fork
                print(f"Iteration {i} - Found a fork - Removed {len(nodes_walked)} nodes", file=f)
                graph.remove_nodes_from(nodes_walked)
                break

            nodes_walked.add(curr_node)

            curr_node = neighbors[0]
            

            


    #We repeat the process for several steps to clear nested stubs
    for i in range(iterations):
        print(f"Iteration {i} - Number of nodes: {len(graph.nodes)}", file=f)
        degree_0_nodes = [node for node, degree in graph.degree if degree == 0]
        graph.remove_nodes_from(degree_0_nodes)


        out_degree_0_nodes = [node for node, degree in graph.out_degree if degree == 0]


        reversed = graph.reverse(copy=False)

        #We go up a stub until we find a node with out degree > 1 (i.e. in degree 1 in the reversed graph)

        for node in out_degree_0_nodes:
            #For nodes with out degreee 0 we walk up the stub
            investigate_stub(node, reversed, i)


        in_degree_0_nodes = [node for node, degree in graph.in_degree if degree == 0]
        for node in in_degree_0_nodes:
            #For nodes with in degreee 0 we walk down the stub
            investigate_stub(node, graph, i)
           


def create_new_res_dir(name):
    import os
    results_dir = "./results"
    dir_name = f"{results_dir}/{name}"

    # Check existing directories
    existing_dirs = [d for d in os.listdir(results_dir) if os.path.isdir(os.path.join(results_dir, d))]
    count = sum(name in dir for dir in existing_dirs)

    # Append count to dir_name if necessary
    if count > 0:
        dir_name = f"{dir_name}_{count}"

    # Create directory
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    else:
        # If the directory already exists, increment the count and create a new one
        while os.path.exists(dir_name):
            count += 1
            dir_name = f"{count}_{dir_name}"
        os.makedirs(dir_name)

    return dir_name
