# em-motif-discovery


## Researched Paper Summaries

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

## Project Steps and Implementation Details

1. **Initial Setup:**
   - Created a virtual environment and installed necessary Python packages (`numpy`, `pandas`, `scipy`, `biopython`, `matplotlib`, `logomaker`, `seaborn`).
   - Structured the repository into folders (`src/`, `data/`, `results/`, `scripts/`, `configs/`, `docs/`, `tests/`).
   - Implemented the EM algorithm based on references from literature (e.g., Bailey & Elkan 1994).

2. **Model Configuration:**
   - Used a `configs/model_params.json` file to store parameters like motif width, maximum EM iterations, and background nucleotide probabilities.
   - Initially set uniform background probabilities `[0.25, 0.25, 0.25, 0.25]` and a fixed motif width (e.g., 8 bp).
   - Initialized the PWM uniformly or with simple random distributions.

3. **Data Preparation:**
   - Cleaned the input sequences (`data/cleaned_sequences.fasta`) ensuring uppercase and no invalid characters.
   - Used either a real dataset or a synthetic dataset with a known motif planted for testing.
   - Downloaded a known motif PFM from JASPAR (e.g., `MA0011.1.pfm`) for comparison.

4. **Running the EM Algorithm:**
   - Implemented `src/em.py` containing `expectation_step`, `maximization_step`, and `run_em_algorithm`.
   - After running EM, saved the discovered PWM to `results/discovered_pwm.npy`.
   - Compared discovered PWM with the known PWM by converting the known PFM to a PWM and computing KL divergence.

5. **Evaluation and Visualization:**
   - Developed `scripts/evaluate_performance.py` to:
     - Parse command-line arguments to specify the discovered PWM file.
     - Compute KL divergence between discovered and known PWMs.
     - Generate a motif logo with `logomaker`.

## Tips and Best Practices

- **Parameter Tuning:**  
  Start with known or literature-supported motif widths. Use background probabilities derived from the data to give a more realistic model. Try increasing `max_iterations` or adjusting the motif_prior if the motif fails to emerge.


## Conclusion

Through careful debugging, parameter adjustments, data quality checks, and code refinements, we addressed numerous issues:

- We ensured that discovered and known PWMs were correctly oriented.
- We verified that the evaluation script used the correct input PWMs.
- We introduced variability in initialization and tested multiple runs.
- We refined parameters and verified that the motif logo is rendered from a properly shaped PWM.

This systematic approach to debugging and refining the code and analysis pipelines helped us move from incorrect results to a more stable and meaningful motif discovery process.




