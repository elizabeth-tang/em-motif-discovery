import numpy as np
import pandas as pd
import logomaker
import os

def load_pwm_from_npy(file_path):
    """
    Load a PWM from a .npy file.
    The PWM should be a numpy array of shape (W,4), corresponding to A,C,G,T.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Discovered PWM file not found at {file_path}")
    pwm = np.load(file_path)
    if pwm.shape[1] != 4:
        raise ValueError("Expected PWM with 4 columns (A,C,G,T).")
    return pwm

def parse_pfm(file_path):
    """
    Parse a PFM (Position Frequency Matrix) file in a JASPAR-like format.
    Expected format:
    >MA0011.1 br
    3.00 5.00 ...
    1.00 2.00 ...
    1.00 1.00 ...
    7.00 4.00 ...
    
    Rows correspond to [A, C, G, T] in order.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"PFM file not found at {file_path}")

    counts = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Header line, skip
                continue
            if line:
                row = list(map(float, line.split()))
                counts.append(row)
    pfm = np.array(counts)  # shape should be (4, W)
    return pfm

def pfm_to_pwm(pfm):
    """
    Convert a Position Frequency Matrix (PFM) to a Position Weight Matrix (PWM)
    by normalizing each column so that sum of (A,C,G,T) = 1.
    """
    col_sums = pfm.sum(axis=0)
    pwm = pfm / (col_sums + 1e-15)
    return pwm

def kl_divergence(p, q):
    """
    Compute the KL divergence between two probability distributions p and q.
    p and q must be arrays of the same length and represent valid probability distributions.
    """
    return np.sum(p * np.log((p + 1e-15) / (q + 1e-15)))

def compare_pwms(discovered_pwm, known_pwm):
    """
    Compute the total KL divergence between discovered and known PWMs.
    This will sum the KL divergence column-by-column across all motif positions.
    Both PWMs must be the same shape.
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
    Generate and save a sequence logo from a PWM using logomaker.
    pwm shape: (W,4) with order A,C,G,T.
    """
    df = pd.DataFrame(pwm, columns=['A','C','G','T'])
    logo = logomaker.Logo(df, shade_below=.5, fade_below=.5, font_name='Arial Rounded MT Bold')
    logo.ax.set_title("Discovered Motif Logo", fontsize=14)
    logo.ax.set_ylabel("Information (bits)", fontsize=12)
    logo.ax.set_xlabel("Position", fontsize=12)
    logo.fig.tight_layout()
    logo.fig.savefig(output_path, dpi=150)
    print(f"Motif logo saved to {output_path}")

if __name__ == "__main__":
    # File paths
    discovered_pwm_path = "results/discovered_pwm.npy"    # The PWM discovered by your EM algorithm
    known_pfm_path = "data/raw/MA0011.1.pfm"              # The known motif PFM from JASPAR
    logo_output_path = "results/discovered_motif_logo.png"

    # Load discovered PWM
    discovered_pwm = load_pwm_from_npy(discovered_pwm_path)

    # Parse the known PFM and convert it to PWM
    known_pfm = parse_pfm(known_pfm_path)
    known_pwm = pfm_to_pwm(known_pfm)

    # Compare the discovered motif with the known motif using KL divergence
    divergence = compare_pwms(discovered_pwm, known_pwm)
    print("KL divergence between discovered and known PWM:", divergence)

    # Generate a sequence logo for the discovered motif
    generate_logo(discovered_pwm, logo_output_path)
