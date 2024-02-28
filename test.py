import utils

def parse_fasta(file_path):
    with open(file_path, 'r') as file:
        contigs = []
        current_contig = []
        for line in file:
            line = line.strip()

            #Hitting a new contig
            if line.startswith('>'):
                contigs.append(''.join(current_contig))
                current_contig = []
            else:
                current_contig.append(line)
        return contigs


contigs = parse_fasta('QuastSequences/k61_simplified_special_assembly.fasta')
print(f"Number of contigs: {len(contigs)}")


print(utils.compute_N50(parse_fasta('QuastSequences/k61_simplified_special_assembly.fasta'), 100000))