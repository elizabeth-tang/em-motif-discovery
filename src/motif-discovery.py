import numpy as np
import statistics


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
def expectation_step(X, motif_probs, bg_probs, lambda_param):
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
        Z[i, 0] = P_motif / total  # Probability for motif
        Z[i, 1] = P_bg / total  # Probability for background

    for i in range(n - W + 1):  # Iterate over all windows of size W
        window_sum_motif = np.sum(Z[i : i + W, 0])  # Sum over motif probabilities
        if window_sum_motif > 1:
            Z[i : i + W, 0] *= 1 / window_sum_motif  # Scale motif probabilities

        window_sum_bg = np.sum(Z[i : i + W, 1])  # Sum over background probabilities
        if window_sum_bg > 1:
            Z[i : i + W, 1] *= 1 / window_sum_bg

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
def run_em(X, alphabet_size, W, max_iters=1000, tol=1e-4):
    """
    Run the EM algorithm to find the motif and background models.
    """
    # Initialize parameters
    motif_probs, bg_probs, lambda_param = initialize_params(alphabet_size, W)
    prev_log_likelihood = -np.inf

    for iteration in range(max_iters):
        # E-step
        Z = expectation_step(X, motif_probs, bg_probs, lambda_param)

        # M-step
        motif_probs, bg_probs, lambda_param = maximization_step(X, Z, alphabet_size, W)

        # Log-likelihood
        log_likelihood = compute_log_likelihood(
            X, Z, motif_probs, bg_probs, lambda_param
        )
        print(f"Iteration {iteration + 1}, Log-Likelihood: {log_likelihood}")

        # Check for convergence
        if np.abs(log_likelihood - prev_log_likelihood) < tol:
            break
        prev_log_likelihood = log_likelihood

    return motif_probs, bg_probs, lambda_param, Z


# Motif discovery process
def discover_motifs(sequences, W, alphabet):
    """
    Discover motifs iteratively using the mixture model.
    """
    # Map sequences to numerical representation
    char_to_index = {char: i for i, char in enumerate(alphabet)}
    sequences_num = [[char_to_index[char] for char in seq] for seq in sequences]
    X = [seq[i : i + W] for seq in sequences_num for i in range(len(seq) - W + 1)]
    alphabet_size = len(alphabet)

    # Run EM to discover motifs
    motif_probs, bg_probs, lambda_param, Z = run_em(X, alphabet_size, W)
    motif = extract_motif_from_probs(motif_probs, alphabet)

    # Return motif probabilities and significant subsequences
    return motif_probs, bg_probs, lambda_param


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
    # sequences = ["ACGTTGAC", "TGCACGTT", "ACGACGTT", "GTTACGTC"]
    sequences = read_fasta("data/cleaned_sequences.fasta")
    # sequences = read_fasta("data/simple_data/test.fasta")
    W = 6  # Width of the motif
    alphabet = ["A", "C", "G", "T"]

    # Discover motifs
    motif_probs, bg_probs, lambda_param = discover_motifs(sequences, W, alphabet)

    # Print results
    print("\nMotif Position Weight Matrix:")
    for row in motif_probs:
        print(" ".join(f"{val:.2f}" for val in row))
    print("\nBackground Frequencies:")
    print(" ".join(f"{val:.2f}" for val in bg_probs))
    print(f"\nMixing Parameter (lambda): {lambda_param:.2f}")
    print(extract_motif_from_probs(motif_probs, alphabet))

    # print(f"Consensus Motif: {motif}")
    # print(read_fasta("data/cleaned_sequences.fasta"))
