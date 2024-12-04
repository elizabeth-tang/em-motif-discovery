"""
Run this file using command:
```
python -c "
from src.model import initialize_model
from src.em import run_em_algorithm
import numpy as np
import os

# Ensure results directory exists
if not os.path.exists('results'):
    os.makedirs('results')

model = initialize_model()
final_model = run_em_algorithm(model, 'data/cleaned_sequences.fasta', W=8, max_iterations=80)

# After EM converges, save the discovered PWM
np.save('results/discovered_pwm.npy', final_model.pwm)
print('Discovered PWM saved to results/discovered_pwm.npy')
"
```
"""
import numpy as np
from Bio import SeqIO
from math import log

def seq_to_numeric(seq, alphabet=['A','C','G','T']):
    alpha_dict = {c:i for i,c in enumerate(alphabet)}
    return np.array([alpha_dict.get(base, -1) for base in seq], dtype=int)

def calc_substring_probability(substring_indices, pwm):
    if np.any(substring_indices < 0):
        return 1e-12
    probs = [pwm[pos, base_idx] for pos, base_idx in enumerate(substring_indices)]
    return np.prod(probs)

def calc_substring_background_prob(substring_indices, background_probs):
    if np.any(substring_indices < 0):
        return 1e-12
    probs = [background_probs[base_idx] for base_idx in substring_indices]
    return np.prod(probs)

def expectation_step(model, sequences, W):
    posteriors = []
    for seq_record in sequences:
        seq_str = str(seq_record.seq)
        seq_idx = seq_to_numeric(seq_str)
        seq_length = len(seq_idx)
        
        position_posteriors = []
        for i in range(seq_length - W + 1):
            motif_like = calc_substring_probability(seq_idx[i:i+W], model.pwm)
            bg_like = calc_substring_background_prob(seq_idx[i:i+W], model.background_probs)
            numerator = motif_like * model.motif_prior
            denominator = numerator + bg_like * (1 - model.motif_prior)
            posterior = numerator / (denominator + 1e-15)
            position_posteriors.append(posterior)
        posteriors.append(np.array(position_posteriors))
    return posteriors

def maximization_step(model, sequences, posteriors, W):
    alphabet_len = model.pwm.shape[1]
    counts = np.zeros((W, alphabet_len))
    total_posterior = 0.0
    total_positions = 0

    for seq_record, pos_array in zip(sequences, posteriors):
        seq_str = str(seq_record.seq)
        seq_idx = seq_to_numeric(seq_str)
        seq_length = len(seq_idx)
        
        for i in range(seq_length - W + 1):
            p = pos_array[i]
            if p > 1e-15:
                substring = seq_idx[i:i+W]
                for pos in range(W):
                    base_idx = substring[pos]
                    if base_idx >= 0:
                        counts[pos, base_idx] += p
                total_posterior += p
        total_positions += (seq_length - W + 1)

    pwm = (counts + 0.001) / (np.sum(counts, axis=1, keepdims=True) + 0.001*alphabet_len)
    new_motif_prior = total_posterior / (total_positions + 1e-15)
    
    model.pwm = pwm
    model.motif_prior = new_motif_prior
    return model

def compute_log_likelihood(model, sequences, W):
    total_ll = 0.0
    for seq_record in sequences:
        seq_str = str(seq_record.seq)
        seq_idx = seq_to_numeric(seq_str)
        seq_length = len(seq_idx)
        for i in range(seq_length - W + 1):
            motif_like = calc_substring_probability(seq_idx[i:i+W], model.pwm)
            bg_like = calc_substring_background_prob(seq_idx[i:i+W], model.background_probs)
            p_data = model.motif_prior * motif_like + (1 - model.motif_prior) * bg_like
            total_ll += np.log(p_data+1e-15)
    return total_ll

def run_em_algorithm(model, sequence_file, W, max_iterations=50, convergence_threshold=1e-5):
    sequences = list(SeqIO.parse(sequence_file, "fasta"))
    prev_ll = None
    for iteration in range(max_iterations):
        posteriors = expectation_step(model, sequences, W)
        model = maximization_step(model, sequences, posteriors, W)
        ll = compute_log_likelihood(model, sequences, W)
        
        print(f"Iteration {iteration}, Log-Likelihood: {ll}")
        
        # Check for convergence
        if prev_ll is not None:
            ll_diff = abs(ll - prev_ll)
            if ll_diff < convergence_threshold:
                print("Convergence reached.")
                break
                
        prev_ll = ll
    return model

