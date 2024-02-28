import unittest
from Node import Node
from utils import build_paths, create_networkx_graph, find_minimizers_in_string, build_assembly, build_walks, walk_graph_upstream
import networkx as nx
import copy


class TestNodes(unittest.TestCase):

    def test_from_encoded_string(self):
        node = Node.from_encoded_string("M_12M_17M_238")
        self.assertEqual(str(node), 'M_12M_17M_238')

        node = Node.from_encoded_string("M_1")
        self.assertEqual(str(node), 'M_1')

        node = Node.from_encoded_string("M_99M_100")
        self.assertEqual(str(node), 'M_99M_100')

        node = Node.from_encoded_string("")
        self.assertEqual(str(node), '')

        node = Node.from_encoded_string("M_")
        self.assertEqual(str(node), '')


class TestBuildPaths(unittest.TestCase):

    def test_empty_reads(self):
        self.assertEqual(build_paths([], 3), {})

    def test_zero_k_mersize(self):
        self.assertEqual(build_paths(["ACTG"], 0), {"ACTG": []})

    def test_k_mersize_gt_read_length(self):
        self.assertEqual(build_paths(["ACT", "CGT"], 4), {"ACT": [], "CGT": []})

    def test_varied_read_lengths(self):
        self.assertEqual(build_paths(["ACT", "CGTA"], 2), {"ACT": ["AC", "CT"], "CGTA": ["CG", "GT", "TA"]})

    def test_typical_case(self):
        self.assertEqual(build_paths(["ACTG", "CGTA"], 3), {"ACTG": ["ACT", "CTG"], "CGTA": ["CGT", "GTA"]})


class TestCreateNetworkxGraph(unittest.TestCase):

    def test_empty_reads(self):
        G = create_networkx_graph([], 3)
        self.assertEqual(G.number_of_nodes(), 0)
        self.assertEqual(G.number_of_edges(), 0)

    def test_single_read_k2(self):
        G = create_networkx_graph(["ACTG"], 2)
        self.assertEqual(set(G.nodes), set(["AC", "CT", "TG"]))
        self.assertEqual(set(G.edges), set([("AC", "CT"), ("CT", "TG")]))
        self.assertEqual(G.nodes["AC"]["node_object"].original_reads, set(["ACTG"]))

    def test_multiple_reads_k2(self):
        G = create_networkx_graph(["ACTG", "GCTA"], 2)
        self.assertEqual(set(G.nodes), set(["AC", "CT", "TG", "GC", "TA"]))
        self.assertEqual(set(G.edges), set([('AC', 'CT'), ('CT', 'TG'), ('GC', 'CT'), ('CT', 'TA')]))
        self.assertEqual(G.nodes["CT"]["node_object"].original_reads, set(["ACTG", "GCTA"]))

    def test_kmer_eq_read_length(self):
        G = create_networkx_graph(["ACTG"], 4)
        self.assertEqual(set(G.nodes), set(["ACTG"]))
        self.assertEqual(set(G.edges), set([]))
        self.assertEqual(G.nodes["ACTG"]["node_object"].original_reads, set(["ACTG"]))

    def test_kmer_gt_read_length(self):
        G = create_networkx_graph(["ACTG"], 5)
        self.assertEqual(G.number_of_nodes(), 0)
        self.assertEqual(G.number_of_edges(), 0)

    def test_multiple_reads_diff_lengths(self):
        G = create_networkx_graph(["ACT", "GCTA"], 2)
        self.assertEqual(set(G.nodes), set(["AC", "CT", "GC", "TA"]))
        self.assertEqual(set(G.edges), set([('AC', 'CT'), ('GC', 'CT'), ('CT', 'TA')]))
        self.assertEqual(G.nodes["CT"]["node_object"].original_reads, set(["ACT", "GCTA"]))

    def test_reads_attributes_all_nodes(self):
        G = create_networkx_graph(["ACT", "GCTA", "TAA"], 2)
        self.assertEqual(G.nodes["CT"]["node_object"].original_reads, set(["ACT", "GCTA"]))
        self.assertEqual(G.nodes["GC"]["node_object"].original_reads, set(["GCTA"]))
        self.assertEqual(G.nodes["TA"]["node_object"].original_reads, set(["GCTA", "TAA"]))
        self.assertEqual(G.nodes["AA"]["node_object"].original_reads, set(["TAA"]))


class TestFindMinimizersInString(unittest.TestCase):

    def test_basic_functionality(self):
        self.assertEqual(find_minimizers_in_string('AAAATTBBBB', ['AAAA', 'BBBB']), ['M_0', 'M_1'])

    def test_no_matches(self):
        self.assertEqual(find_minimizers_in_string('CCCCCDDDDD', ['AAAA', 'BBBB']), [])

    def test_partial_matches(self):
        self.assertEqual(find_minimizers_in_string('AAAACCCC', ['AAAA', 'BBBB']), ['M_0'])

    def test_repeated_minimizers(self):
        self.assertEqual(find_minimizers_in_string('AAAAAAAABBBB', ['AAAA', 'BBBB']), ['M_0', 'M_0', 'M_1'])
    def test_repeated_minimizers_2(self):
        self.assertEqual(find_minimizers_in_string('AAAAAAABBBB', ['AAAA', 'BBBB']), ['M_0', 'M_1'])

    def test_empty_string(self):
        self.assertEqual(find_minimizers_in_string('', ['AAAA', 'BBBB']), [])

    def test_empty_minimizers_set(self):
        self.assertEqual(find_minimizers_in_string('AAAABBBB', []), [])

    def test_overlapping_minimizers(self):
        self.assertEqual(find_minimizers_in_string('AAAABBBB', ['AAA', 'BBB']), ['M_0', 'M_1'])

    def test_special_characters(self):
        self.assertEqual(find_minimizers_in_string('123@#$$%456', ['123', '@#$']), ['M_0', 'M_1'])

    def test_case_sensitivity(self):
        self.assertEqual(find_minimizers_in_string('aaabbBb', ['AAA', 'BBB']), [])






class TestBuildAssembly(unittest.TestCase):

    def test_basic_assembly_no_special_edges(self):
        # Create a simple graph
        g_line = nx.DiGraph()

        minimizer_set_line= ['AB', 'BC', 'CD', 'DE', 'EF', 'FG', 'GH']  # Example minimizer set

        g_line.add_node('ABC', minimizer_encoding=['M_0'])
        g_line.add_node('BCD', minimizer_encoding=['M_1'])
        g_line.add_edge('ABC', 'BCD', special=False)

      
        # Call build_assembly
        assembly = build_assembly(g_line, None, minimizer_set_line, basic=True)

        print(assembly)
        # Define the expected assembly
        expected_assembly = ['ABCD']  # This is an assumption
        self.assertEqual(assembly, expected_assembly)

    def test_basic_assembly_with_special_edges(self):
        g_line = nx.DiGraph()

        minimizer_set_line= ['AB', 'BC', 'CD', 'DE', 'EF', 'FG', 'GH']  # Example minimizer set

        g_line.add_node('ABC', minimizer_encoding='M_0')
        g_line.add_node('BCD', minimizer_encoding='M_1')
        g_line.add_edge('ABC', 'BCD', special=False)

        g_line.add_node('DBCXX', minimizer_encoding='M_1')
        g_line.add_edge('BCD', 'DBCXX', special=True)

        # Call build_assembly
        assembly = build_assembly(g_line, None, minimizer_set_line, basic=False)

        print(assembly)
        # Define the expected assembly
        expected_assembly = ['ABCDX']  # This is an assumption
        self.assertEqual(assembly, expected_assembly)

class TestGraphWalk(unittest.TestCase):

    def test_walks_on_graph_line(self):

        """
        a -> b -> c -> d -[special]-> e -> f -> g -> h
        
        """

        g_line_normal = nx.DiGraph()

        g_line_normal.add_node('A')
        g_line_normal.add_node('B')
        g_line_normal.add_node('C')
        g_line_normal.add_node('D')
        g_line_normal.add_node('E')
        g_line_normal.add_node('F')

        g_line_normal.add_edge('A', 'B', special=False)
        g_line_normal.add_edge('B', 'C', special=False)
        g_line_normal.add_edge('C', 'D', special=False)
        g_line_normal.add_edge('E', 'F', special=False)
        g_line_normal.add_edge('F', 'G', special=False)
        g_line_normal.add_edge('G', 'H', special=False)


        g_line_special = copy.deepcopy(g_line_normal)
        g_line_special.add_edge('D', 'E', special=True)

        # Call build_walks
        walks_special = set([''.join(x) for x in build_walks(g_line_special, set(('D', 'E')), True)])
        walks_no_special = set([''.join(x) for x in build_walks(g_line_normal, set(), False)])

        print(f"walks_special: {walks_special}")
        print(f"walks_no_special: {walks_no_special}")
        self.assertEqual(walks_special, {''.join(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'])})
        self.assertEqual(walks_no_special, {''.join(['A', 'B', 'C', 'D']), ''.join(['E', 'F', 'G', 'H'])})


    def test_walk_upstream_1(self):

        g = nx.DiGraph()

        """
        
       a -> b -> c -> d -> i --[special]--> j
                           ̂̉/
                          /   
        e -> f -> g -> h 
        
        """

        g.add_node('A')
        g.add_node('B')
        g.add_node('C')
        g.add_node('D')
        g.add_node('E')
        g.add_node('F')
        g.add_node('G')
        g.add_node('H')
        g.add_node('I')
        g.add_node('J')

        g.add_edge('A', 'B', special=False)
        g.add_edge('B', 'C', special=False)
        g.add_edge('C', 'D', special=False)
        g.add_edge('E', 'F', special=False)
        g.add_edge('F', 'G', special=False)
        g.add_edge('G', 'H', special=False)
        g.add_edge('D', 'I', special=False)
        g.add_edge('I', 'J', special=True)
        g.add_edge('H', 'I', special=False)

        # Call build_walks
        walk = walk_graph_upstream(g, 'I', set())
        self.assertEqual(walk, ['I', 'D', 'C', 'B', 'A'])


if __name__ == "__main__":
    unittest.main()