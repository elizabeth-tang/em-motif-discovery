"""
RUn using `python src/model.py
"""
import json
import numpy as np

class FiniteMixtureModel:
    def __init__(self, pwm, background_probs, motif_prior):
        self.pwm = pwm
        self.background_probs = background_probs
        self.motif_prior = motif_prior

def initialize_model(config_path="configs/model_params.json"):
    with open(config_path) as f:
        params = json.load(f)
    W = params["motif_width"]
    alph_len = len(params["alphabet"])

    pwm = np.random.random((W, alph_len))
    pwm = pwm / pwm.sum(axis=1, keepdims=True)

    background_probs = np.array(params["background_probs"])
    motif_prior = 0.05
    model = FiniteMixtureModel(pwm, background_probs, motif_prior)
    return model

if __name__ == "__main__":
    model = initialize_model()
    print("Initial PWM:\n", model.pwm)
    print("Background:", model.background_probs)
    print("Motif prior:", model.motif_prior)
