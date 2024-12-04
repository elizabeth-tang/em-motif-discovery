import random

alphabet = ['A','C','G','T']
motif = "ACGTA"  # known motif for synthetic test

with open("data/raw/synthetic_sequences.fasta", "w") as f:
    for i in range(50):
        seq = ''.join(random.choices(alphabet, k=100))
        # Embed motif at random position
        start = random.randint(0, 95)
        seq = seq[:start] + motif + seq[start+5:]
        f.write(f">synthetic_seq_{i}\n{seq}\n")
