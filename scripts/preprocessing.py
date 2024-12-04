from Bio import SeqIO

input_file = "data/raw/human_promoters.fasta"
output_file = "data/cleaned_sequences.fasta"

with open(output_file, "w") as fout:
    for record in SeqIO.parse(input_file, "fasta"):
        seq_str = str(record.seq).upper().replace('N','A')  # Replace 'N' with 'A' arbitrarily
        fout.write(f">{record.id}\n{seq_str}\n")
