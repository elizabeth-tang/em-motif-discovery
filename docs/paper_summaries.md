## Paper Summaries

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


