from Bio import SeqIO

input_file = "data/raw/gene.fna"
output_file = "data/cleaned_sequences.fasta"

with open(output_file, "w") as fout:
    for record in SeqIO.parse(input_file, "fasta"):
        seq_str = str(record.seq).upper()
        # Replace any 'N' or invalid chars with 'A' (arbitrary choice)
        seq_str = ''.join([base if base in 'ACGT' else 'A' for base in seq_str])
        fout.write(f">{record.id}\n{seq_str}\n")
