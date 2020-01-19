import random
from sys import exit
from time import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb
from pyDOE import lhs  # pip install pyDOE

from utils import eval_signal, oscilating, find_peaks


def pairplot(samples_number, m_range, K_range, base_amplituda, base_perioda):
    ms = np.linspace(m_range[0], m_range[1], m_range[1] - m_range[0] + 1, dtype=int)
    Ks = np.linspace(10 ** K_range[0], 10 ** K_range[1], samples_number)
    oscilating_samples = []
    print(base_amplituda, base_perioda)
    while len(oscilating_samples) < 10:
        new_params = list(params)
        new_params[4] = random.choice(ms)
        new_params[5] = random.choice(ms)
        new_params[6] = random.choice(ms)
        new_params[7] = random.choice(Ks)
        new_params[8] = random.choice(Ks)
        new_params[9] = random.choice(Ks)

        A, B, C = simulate(model, tuple(new_params))
        peaks_A, minimums_A, peaks_vals_A, minimums_vals_A = find_peaks(A)

        if oscilating(peaks_vals_A, 0.01):
            amplituda, perioda = eval_signal(peaks_A, minimums_vals_A, peaks_vals_A)
            if amplituda >= base_amplituda and perioda <= base_perioda:
                oscilating_samples.append(list(new_params[4:]))
                print(len(oscilating_samples), amplituda, perioda)

    df = pd.DataFrame(np.array(oscilating_samples), columns=["m1", "m2", "m3", "K1", "K2", "K3"])
    sb.pairplot(df, markers="+")


# util function
def denormalize(normalized_value, min, max, discrete=False):
    if discrete:
        values = np.linspace(min, max, max - min + 1, dtype=int)
        splitter = 1 / (max - min + 1)
        start_splitter = splitter
        counter = 0
        while counter < (max - min):
            if normalized_value < start_splitter:
                return values[counter]
            counter += 1
            start_splitter += splitter
        return values[-1]
    else:
        value = ((max - min) * normalized_value) + min
        return value


def pairplot_lhs(samples_number, m_range, K_range, base_amplituda, base_perioda):
    samples = lhs(6, samples=samples_number, criterion='center')
    oscilating_samples = []
    new_params = list(params)
    better_both_counter = 0
    better_amp_counter = 0
    better_per_counter = 0
    sample_counter = 0

    for sample in samples:

        # denormalize
        new_params[4] = denormalize(sample[0], m_range[0], m_range[1], discrete=True)  # m1
        new_params[5] = denormalize(sample[1], m_range[0], m_range[1], discrete=True)  # m2
        new_params[6] = denormalize(sample[2], m_range[0], m_range[1], discrete=True)  # m3
        new_params[7] = denormalize(sample[3], K_range[0], K_range[1], discrete=False)  # K1
        new_params[8] = denormalize(sample[4], K_range[0], K_range[1], discrete=False)  # K2
        new_params[9] = denormalize(sample[5], K_range[0], K_range[1], discrete=False)  # K3

        A, B, C = simulate(model, tuple(new_params))
        peaks_A, minimums_A, peaks_vals_A, minimums_vals_A = find_peaks(A)

        if oscilating(peaks_vals_A, 0.01):
            amplituda, perioda = eval_signal(peaks_A, minimums_vals_A, peaks_vals_A)

            if amplituda > 0 and perioda > 0:
                if amplituda >= base_amplituda and perioda <= base_perioda:
                    oscilating_samples.append(list(new_params[4:]) + ["Boljše oboje"])
                    better_both_counter += 1
                elif amplituda > base_amplituda:
                    oscilating_samples.append(list(new_params[4:]) + ["Boljša amplituda"])
                    better_amp_counter += 1
                elif perioda < base_perioda:
                    oscilating_samples.append(list(new_params[4:]) + ["Boljša perioda"])
                    better_per_counter += 1

        sample_counter += 1
        print("Progress: " + str(sample_counter / samples_number))

    print("Base Values:", base_amplituda, base_perioda)
    print("Samples: ", samples_number)
    print("Better amp: ", better_amp_counter, "Better per: ", better_per_counter, "Better both: ", better_both_counter)

    df = pd.DataFrame(oscilating_samples, columns=["m1", "m2", "m3", "K1", "K2", "K3", "Primerjava z osnovo"])
    if better_both_counter > 0:
        sb.pairplot(df, markers="+", hue="Primerjava z osnovo",
                    hue_order=["Boljša perioda", "Boljša amplituda", "Boljše oboje"])
    else:
        sb.pairplot(df, markers="+", hue="Primerjava z osnovo", hue_order=["Boljša perioda", "Boljša amplituda"])
    plt.show()


# USER PARAMETERS
model_type = 'positive'
# model_type = 'negative'
# END OF USER PARAMETERS

# SETUP
start_time = time()

_temp = __import__(model_type + '_core', globals(), locals(), ['get_model_params', 'model', 'simulate'], 0)
get_model_params = _temp.get_model_params
model = _temp.model
simulate = _temp.simulate

params = get_model_params()
# END OF SETUP

# base model vrednosti
A, B, C = simulate(model, params)
peaks_A, minimums_A, peaks_vals_A, minimums_vals_A = find_peaks(A)
base_amplituda, base_perioda = eval_signal(peaks_A, minimums_vals_A, peaks_vals_A)

# pairplot(10000, [1, 4], [-1, 3], base_amplituda, base_perioda)
pairplot_lhs(200000, [1, 4], [0.1, 1000], base_amplituda, base_perioda)

print('Job completed in ', time() - start_time, ' seconds.')
exit(0)
