"""
# Run and check the KL divergence between discovered and known PWM:

## First step:
for i in {1..5}; do
    python src/em.py
    cp results/discovered_pwm.npy results/discovered_pwm_run$i.npy
done

## Second step:
for i in {1..5}; do
    python scripts/evaluate_performance.py --input results/discovered_pwm_run$i.npy
done
"""

import argparse
import numpy as np
import pandas as pd
import logomaker
import os

def load_pwm_from_npy(file_path):
    """
    Load a PWM from a .npy file.
    The PWM should be a numpy array of shape (W,4).
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"PWM file not found at {file_path}")
    pwm = np.load(file_path)
    if pwm.shape[1] != 4:
        raise ValueError("Expected PWM with 4 columns (A,C,G,T).")
    return pwm

def parse_pfm(file_path):
    """
    Parse a PFM (Position Frequency Matrix) file in the JASPAR-like format.
    Rows correspond to A, C, G, T in that order.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"PFM file not found at {file_path}")

    counts = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                continue
            if line:
                row = list(map(float, line.split()))
                counts.append(row)
    pfm = np.array(counts)  # shape: (4, W)
    return pfm

def pfm_to_pwm(pfm):
    """
    Convert a PFM to a PWM by normalizing each column to sum to 1.
    """
    col_sums = pfm.sum(axis=0)
    pwm = pfm / (col_sums + 1e-15)
    return pwm

def kl_divergence(p, q):
    """
    Compute KL divergence between distributions p and q.
    """
    return np.sum(p * np.log((p+1e-15)/(q+1e-15)))

def compare_pwms(discovered_pwm, known_pwm):
    """
    Compute total KL divergence column by column.
    Ensure that both PWMs have shape (W,4).
    """
    if discovered_pwm.shape != known_pwm.shape:
        raise ValueError("Discovered PWM and known PWM must have the same shape.")

    total_div = 0.0
    for i in range(discovered_pwm.shape[0]):
        p_col = discovered_pwm[i]
        q_col = known_pwm[i]
        total_div += kl_divergence(p_col, q_col)
    return total_div

def generate_logo(pwm, output_path="results/discovered_motif_logo.png"):
    """
    Generate a sequence logo from a PWM and save as PNG.
    pwm: shape (W,4), columns = A,C,G,T
    """
    df = pd.DataFrame(pwm, columns=['A','C','G','T'])
    # Using a commonly available font:
    logo = logomaker.Logo(df)
    logo.ax.set_title("Discovered Motif Logo", fontsize=14)
    logo.ax.set_ylabel("Information (bits)", fontsize=12)
    logo.ax.set_xlabel("Position", fontsize=12)
    logo.fig.tight_layout()
    logo.fig.savefig(output_path, dpi=150)
    print(f"Motif logo saved to {output_path}")

if __name__ == "__main__":
    # Argument parsing for input PWM file
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, required=True, help='Path to discovered PWM file')
    args = parser.parse_args()

    discovered_pwm_path = args.input

    # Known PFM path (adjust if needed)
    known_pfm_path = "data/raw/MA0011.1.pfm"

    # Load discovered PWM from given input path
    discovered_pwm = load_pwm_from_npy(discovered_pwm_path)

    # Parse known PFM and convert to PWM, then transpose to (W,4)
    known_pfm = parse_pfm(known_pfm_path)
    known_pwm = pfm_to_pwm(known_pfm)
    known_pwm = known_pwm.T  # Transpose so shape matches discovered_pwm

    # Compute KL divergence
    divergence = compare_pwms(discovered_pwm, known_pwm)
    print("KL divergence between discovered and known PWM:", divergence)

    # Generate motif logo for the discovered PWM
    generate_logo(discovered_pwm, "results/discovered_motif_logo.png")
