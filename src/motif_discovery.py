import numpy as np
from matplotlib import pyplot as plt


# Helper function to initialize parameters
def initialize_params(alphabet_size, W):
    """
    Initialize the motif and background model parameters and mixing coefficient.
    """
    motif_probs = np.random.dirichlet([1] * alphabet_size, W)
    bg_probs = [0.25, 0.25, 0.25, 0.25]
    lambda_param = 0.5  # Initial guess for the mixing coefficient
    return motif_probs, bg_probs, lambda_param


# Expectation step
def expectation_step(X, motif_probs, bg_probs, lambda_param, seq_lens):
    """
    Compute the expected membership probabilities (E-step).
    """
    n = len(X)
    Z = np.zeros((n, 2))  # Initialize Z as a 2D array

    # Compute raw Z values for motif and background
    for i, subseq in enumerate(X):
        P_motif = lambda_param * np.prod(
            [motif_probs[j, char] for j, char in enumerate(subseq)]
        )
        P_bg = (1 - lambda_param) * np.prod([bg_probs[char] for char in subseq])

        # Normalize probabilities to sum to 1
        total = P_motif + P_bg
        if total > 0:
            Z[i, 0] = P_motif / total  # Probability for motif
            Z[i, 1] = P_bg / total  # Probability for background
        else:
            Z[i, 0] = 0
            Z[i, 1] = 1  # This case should not normally occur

    # Ensure global normalization across rows (each subsequence)
    row_sums = np.sum(Z, axis=1)
    Z /= row_sums[:, None]  # Normalize each row to sum to 1

    assert np.allclose(np.sum(Z, axis=1), 1), "Probabilities not normalized correctly."

    return Z


# Maximization step
def maximization_step(X, Z, alphabet_size, W):
    """
    Update the motif and background model parameters and mixing coefficient (M-step).
    """
    lambda_param = np.mean(Z[:, 0])
    motif_probs = np.zeros((W, alphabet_size))
    bg_probs = np.zeros(alphabet_size)

    # Update motif probabilities
    for j in range(W):
        for i in range(alphabet_size):
            motif_probs[j, i] = sum(
                Z[k, 0] for k, subseq in enumerate(X) if subseq[j] == i
            )
    motif_probs = motif_probs / motif_probs.sum(axis=1, keepdims=True)

    # Update background probabilities
    for i in range(alphabet_size):
        bg_probs[i] = sum(
            (Z[k, 1]) * sum(1 for char in subseq if char == i)
            for k, subseq in enumerate(X)
        )
    bg_probs /= bg_probs.sum()

    return motif_probs, bg_probs, lambda_param


# Log-likelihood calculation
def compute_log_likelihood(X, Z, motif_probs, bg_probs, lambda_param, eps=1e-6):
    """
    Compute the log-likelihood of the data given the current parameters.
    """
    log_likelihood = 0
    for i, subseq in enumerate(X):
        P_motif = lambda_param * np.prod(
            [motif_probs[j, char] for j, char in enumerate(subseq)]
        )
        P_bg = (1 - lambda_param) * np.prod([bg_probs[char] for char in subseq])

        # Add contributions from motif and background components
        log_likelihood += Z[i, 0] * np.log(P_motif + eps) + Z[i, 1] * np.log(P_bg + eps)

    return log_likelihood


# Main EM algorithm
def run_em(X, alphabet_size, W, seq_lens, max_iters=500, tol=1e-4):
    """
    Run the EM algorithm to find the motif and background models.
    """
    # Initialize parameters
    motif_probs, bg_probs, lambda_param = initialize_params(alphabet_size, W)
    log_likelihoods = []
    prev_log_likelihood = -np.inf

    for iteration in range(max_iters):
        # E-step
        Z = expectation_step(X, motif_probs, bg_probs, lambda_param, seq_lens)

        # M-step
        motif_probs, bg_probs, lambda_param = maximization_step(X, Z, alphabet_size, W)

        # Log-likelihood
        log_likelihood = compute_log_likelihood(
            X, Z, motif_probs, bg_probs, lambda_param
        )
        print(f"Iteration {iteration + 1}, Log-Likelihood: {log_likelihood}")
        log_likelihoods.append(log_likelihood)

        # Check for convergence
        if np.abs(log_likelihood - prev_log_likelihood) < tol:
            break
        prev_log_likelihood = log_likelihood
        print(extract_motif_from_probs(motif_probs, alphabet))

    return motif_probs, bg_probs, lambda_param, Z, log_likelihoods


# Motif discovery process
def discover_motifs(sequences, W, alphabet):
    """
    Discover motifs iteratively using the mixture model.
    """
    # Map sequences to numerical representation
    char_to_index = {char: i for i, char in enumerate(alphabet)}
    sequences_num = [[char_to_index[char] for char in seq] for seq in sequences]
    X = []
    seq_lens = []
    for seq in sequences_num:
        for i in range(len(seq) - W + 1):
            X.append(seq[i : i + W])
        seq_lens.append(len(seq) - W + 1)
    alphabet_size = len(alphabet)

    # Run EM to discover motifs
    motif_probs, bg_probs, lambda_param, Z, log_likelihoods = run_em(
        X, alphabet_size, W, seq_lens
    )
    motif = extract_motif_from_probs(motif_probs, alphabet)

    # Return motif probabilities and significant subsequences
    return motif_probs, bg_probs, lambda_param, log_likelihoods


def extract_motif_from_probs(motif_probs, alphabet):
    """
    Extract the consensus motif from the motif probabilities matrix.

    Parameters:
        motif_probs (numpy.ndarray): The motif probability matrix of shape (W, alphabet_size).
        alphabet (list): List of characters in the alphabet (e.g., ['A', 'C', 'G', 'T']).

    Returns:
        str: The consensus motif as a string.
    """
    consensus_motif = ""
    for position_probs in motif_probs:
        # Find the index of the maximum probability for this position
        max_index = np.argmax(position_probs)
        # Map the index back to the corresponding character
        consensus_motif += alphabet[max_index]
    return consensus_motif


def read_fasta(filename):
    with open(filename, "r") as f:
        output = []
        s = ""
        for l in f.readlines():
            if l.strip()[0] == ">":
                # skip the line that begins with ">"
                if s == "":
                    continue
                output.append(s)
                s = ""
            # keep appending to the current sequence when the line doesn't begin
            # with ">"
            else:
                s += l.strip()
        output.append(s)
        return output


# Example usage
if __name__ == "__main__":
    # Input sequences

    # Test/Artificial Sequences
    # sequences = ["ACGTTGAC", "TGCACGTT", "ACGACGTT", "GTTACGTT"]  # W = 5
    # sequences = read_fasta("data/simple_data/test.fasta") # W = 4
    # sequences = read_fasta("data/simple_data/motif.fasta") # W = 10

    # Real Sequences
    sequences = read_fasta("data/real_data/U10081.1.fasta")  # W = 11
    # sequences = read_fasta("data/real_data/MA0011.1.fasta") # W = 6, 8
    W = 11  # Width of the motif
    trials = 10
    alphabet = ["A", "C", "G", "T"]

    # Discover motifs
    discovered_motifs = []

    for i in range(trials):
        print(f"Running Trial 1:")
        discovered_motifs.append(discover_motifs(sequences, W, alphabet))

    motif_probs, bg_probs, lambda_param, log_likelihoods = max(
        discovered_motifs, key=(lambda x: x[3][-1])
    )

    # Print results
    print("\nMotif Position Weight Matrix:")
    for row in motif_probs:
        print(" ".join(f"{val:.2f}" for val in row))
    print("\nBackground Frequencies:")
    print(" ".join(f"{val:.2f}" for val in bg_probs))
    print(f"\nMixing Parameter (lambda): {lambda_param:.2f}")
    print(extract_motif_from_probs(motif_probs, alphabet))

    np.save("results/discovered_pwm.npy", motif_probs)
    print("Discovered PWM saved to results/discovered_pwm.npy")

    plt.plot(log_likelihoods)
    plt.title("Convergence")
    plt.xlabel("Iterations")
    plt.ylabel("Log Likelihood")
    plt.savefig("results/log_likelihood.png")
