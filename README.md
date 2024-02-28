# Minimizer Schemes For Graph Assembly 

The following codebase is used to test the effects of minimizers encoding on DBJ graph assembly.

All experiments can be run by executing main.py with the correct arguments

example command: ```python3 main.py --input_sequence_path input_sequences/random_seq_l-100000.fna --k 61 --error_rate 0.03 --coverage 30 --number_minimizer_sets 1 --special_edge_first```

The resulting assemblies are all saved to disk along with metrics. For detailed comparison we recommend using Quast.
