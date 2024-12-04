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
   - Created a virtual environment and installed necessary Python packages (`numpy`, `pandas`, `scipy`, `biopython`, `matplotlib`, `logomaker`).
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

## Problems Encountered and Their Solutions

1. **Shape Mismatch Between Known and Discovered PWMs:**
   - **Issue:** We initially got a `ValueError` due to mismatched PWM shapes. The discovered PWM was `(W,4)` while the known PWM was `(4,W)`.
   - **Solution:** Transposed the known PWM to `(W,4)` so both PWMs align correctly. Updated the code to confirm correct shapes before comparison.

2. **Same KL Divergence and Same Results After Multiple Runs:**
   - **Issue:** Even after running the EM algorithm multiple times, the KL divergence remained identical, and the discovered motif did not change.
   - **Potential Causes and Fixes:**
     - **Lack of True Random Initialization:** Ensured the PWM initialization in `src/model.py` uses `np.random.random` and no fixed seed.
     - **Argument Parsing in `evaluate_performance.py`:** Confirmed that `--input` parameter properly selects different PWM files. Updated the code to parse arguments using `argparse` and tested multiple runs again.

3. **Strange Motif Logo (Rows of Single Nucleotides):**
   - **Issue:** The motif logo displayed a pattern where each position seemed dominated by a single nucleotide row (e.g., top row all A’s, next row all C’s), not a realistic motif distribution.
   - **Possible Causes:**
     - Incorrect matrix orientation when generating logos.
     - The EM algorithm converging to a trivial or non-informative solution due to non-ideal parameters or input data.
   - **Solutions:**
     - Verified that the PWM has shape `(W,4)` before passing it to `logomaker`.
     - Printed intermediate PWMs during the EM steps to ensure they update meaningfully.
     - Tuned parameters (motif width, motif_prior, background_probs) to better reflect the underlying biology.
     - Checked the dataset to ensure it contains a strong motif signal. Considered testing on a synthetic dataset with a known motif.
     - Ensured true randomness in initialization and possibly tested different subsets of the sequence data.

4. **Font Warnings for Sequence Logo:**
   - **Issue:** Multiple warnings about `Arial Rounded MT Bold` font not found.
   - **Solution:** Changed the font to a commonly available default font (e.g., `DejaVu Sans`), or omitted `font_name` entirely in `logomaker.Logo()` calls.

## Tips and Best Practices

- **Debugging EM Steps:**  
  Add print statements after each iteration of EM to track changes in log-likelihood, PWM values, and posterior probabilities. This helps identify if parameters are updating or if the algorithm is stuck.

- **Parameter Tuning:**  
  Start with known or literature-supported motif widths. Use background probabilities derived from the data to give a more realistic model. Try increasing `max_iterations` or adjusting the motif_prior if the motif fails to emerge.

- **Data Quality Checks:**  
  Ensure the input dataset contains the motif you’re trying to discover. Test on a synthetic dataset where you know the motif and see if the EM algorithm can recover it. If not, there may be a coding or conceptual issue.

- **Multiple Runs and Random Seeds:**
  Run the algorithm multiple times with random initializations to avoid converging on a single local maximum. Ensure no fixed random seed is forcing the same initial conditions each run.

## Conclusion

Through careful debugging, parameter adjustments, data quality checks, and code refinements, we addressed numerous issues:

- We ensured that discovered and known PWMs were correctly oriented.
- We verified that the evaluation script used the correct input PWMs.
- We introduced variability in initialization and tested multiple runs.
- We refined parameters and verified that the motif logo is rendered from a properly shaped PWM.

This systematic approach to debugging and refining the code and analysis pipelines helped us move from incorrect results to a more stable and meaningful motif discovery process.




