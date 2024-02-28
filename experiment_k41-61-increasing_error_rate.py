import os
import sys
import subprocess

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

error_rates = [0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10 , 0.11, 0.12, 0.13, 0.14, 0.15]
number_runs_per_experiment = 1


def command(k, error_rate):
    return f"python3 main.py --input_sequence_path input_sequences/random_seq_l-100000.fna --k {k} --error_rate {error_rate} --coverage 30"

print("Running experiments for k=41 and k=61 with increasing error rate")
for k in [41, 61]:   
    print(f"Running experiments for k={k}") 
    for error_rate in error_rates:
        print(f"Running experiments for error_rate={error_rate}")
        for run in range(number_runs_per_experiment):
            print(f"Running experiment {run+1}/{number_runs_per_experiment}")

            sh = subprocess.run(command(k, error_rate), shell=True, capture_output=True)