# em-motif-discovery


## Research Papers

### Bailey & Elkan (1994)
- Introduces MEME algorithm using EM for motif discovery.
- Finite mixture model: sequences considered as mixture of motif vs background.
- Uses position weight matrices (PWMs) and updates them iteratively.
- Initial guess of motif and updates parameters to maximize likelihood.

### Lawrence & Reilly (1990)
- Original formulation of EM-based motif discovery.
- EM steps:
  E-step: Compute posterior probability of each position being a motif instance.
  M-step: Re-estimate PWM and mixture parameters based on these probabilities.
- Convergence based on stabilizing likelihood.


## Introduction to the Project

This project implements a motif discovery algorithm using the Expectation-Maximization (EM) approach with a finite mixture model. The goal is to identify a biologically relevant motif in a given set of DNA sequences and then compare the discovered motif to a known reference motif from databases like JASPAR. The discovered motif is represented as a Position Weight Matrix (PWM) and visualized using sequence logos.

Throughout the development and execution of this project, we encountered various challenges. This README documents the process, the problems we faced, the debugging steps taken, and the solutions implemented. It serves as a detailed record of the coding and troubleshooting journey.

## How to Run:

1. **Setup:**
  - Created a virtual environment and installed necessary Python packages (`numpy`, `pandas`, `scipy`, `biopython`, `matplotlib`, `logomaker`, `seaborn`).
  - In `src/motif_discovery`, modify the below fields to the appropriate data/settings that you want to test on:
    - sequences: the DNA sequences to extract the motif from. A few options are already available
    - width: the width of the motif to find
    - trials: the number of trials to run
2. **Running Expectation-Maximization Algorithm:**
  - Run `python src/motif_discovery.py`
3. **Visualization:**
  - Copy and paste the Motif Position Weight Matrix into the data_str field of `src/visualization.ipynb`
  - Run the below cells to create the motif logo and a heatmap of the PWM


## Project Steps and Implementation Details

1. **Initial Setup:**
  - Created a virtual environment and installed necessary Python packages (`numpy`, `pandas`, `scipy`, `biopython`, `matplotlib`, `logomaker`, `seaborn`).
  - Structured the repository into folders (`src/`, `data/`, `results/`, `scripts/`, `configs/`, `docs/`, `tests/`).
  - Implemented the EM algorithm based on references from literature (e.g., Bailey & Elkan 1994).

2. **Data Preparation:**
  - Cleaned the input sequences (`data/cleaned_sequences.fasta`) ensuring uppercase and no invalid characters.
  - Used either a real dataset or a synthetic dataset with a known motif planted for testing.
  - Downloaded a known motif PFM from JASPAR (e.g., `MA0011.1.pfm`) for comparison.

3. **Running the EM Algorithm:**
  - Implemented `src/motif_discovery.py` which contains the EM algorithm
  - After running EM, saved the discovered PWM to `results/discovered_pwm.npy`.
  - Compared discovered PWM with the known PWM by converting the known PFM to a PWM and computing KL divergence.

5. **Evaluation and Visualization:**
  - Developed `scripts/evaluate_performance.py` to:
    - Parse command-line arguments to specify the discovered PWM file.
    - Compute KL divergence between discovered and known PWMs.
  - Implemented `src/visualization.ipynb` to:
    - Display a heatmap for any given PWM
    - Create a motif logo for any given PWM

## Tips and Best Practices

- **Parameter Tuning:**  
  Start with known or literature-supported motif widths. Use background probabilities derived from the data to give a more realistic model. Try increasing `max_iterations` or adjusting the motif_prior if the motif fails to emerge.




