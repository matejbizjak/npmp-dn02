from datetime import datetime
from sys import exit
from time import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb

from utils import eval_signal, oscilating, find_peaks


def preveri_obmocje(samples_number):
    # params = (alpha, alpha0, beta, n, m1, m2, m3, K1, K2, K3)
    new_params = list(params)

    # zgradimo prazne matrike, kamor bomo shranjevali rezultate
    oscilacije = [[np.nan for x in range(samples_number)] for x in range(4)]
    amplitude = [[np.nan for x in range(samples_number)] for x in range(4)]
    periode = [[np.nan for x in range(samples_number)] for x in range(4)]

    # sampling
    m_samples = [4, 3, 2, 1]
    K_samples = np.linspace(0.1, 1000, samples_number)

    # za vse sample pozenemo simulacijo in rezultate shranimo v prej narejene matrike
    for x in range(len(m_samples)):
        for y in range(len(K_samples)):
            new_params[4] = m_samples[x]
            new_params[5] = m_samples[x]
            new_params[6] = m_samples[x]
            new_params[7] = K_samples[y]
            new_params[8] = K_samples[y]
            new_params[9] = K_samples[y]

            A, B, C = simulate(model, tuple(new_params))

            peaks_A, minimums_A, peaks_vals_A, minimums_vals_A = find_peaks(A)

            if oscilating(peaks_vals_A, 0.01):
                oscilacije[x][y] = 1

                amplituda, perioda = eval_signal(peaks_A, minimums_vals_A, peaks_vals_A)
                amplitude[x][y] = amplituda
                periode[x][y] = perioda

    # narisemo grafe
    plot_graph(oscilacije, m_samples, K_samples, fig_osc, 'Oscilacije')
    plot_graph(amplitude, m_samples, K_samples, fig_amp, 'Amplitude')
    plot_graph(periode, m_samples, K_samples, fig_per, 'Periode')


def plot_graph(data, m_samples, K_samples, figure, figure_title):
    data_frame = pd.DataFrame(data, columns=list(map(lambda k: int(k), K_samples)), index=m_samples)
    subplot = figure.add_subplot(1, 1, 1)
    figure.suptitle(figure_title, fontsize=30)
    sb.heatmap(data_frame, ax=subplot, cbar=True, square=False, mask=data_frame.isnull(), cmap="YlGnBu")
    subplot.set_xlabel("K1 = K2 = K3")
    subplot.set_ylabel("m1 = m2 = m3")


# USER PARAMETERS
model_type = 'positive'
# model_type = 'negative'

figure_size_x = 15
figure_size_y = 15

sb.set(font_scale=2)

samples_number = 250
# END OF USER PARAMETERS

# SETUP
start_time = time()

_temp = __import__(model_type + '_core', globals(), locals(), ['get_model_params', 'model', 'simulate'], 0)
get_model_params = _temp.get_model_params
model = _temp.model
simulate = _temp.simulate

params = get_model_params()

fig_osc = plt.figure(figsize=(figure_size_x, figure_size_y))
fig_amp = plt.figure(figsize=(figure_size_x, figure_size_y))
fig_per = plt.figure(figsize=(figure_size_x, figure_size_y))
# END OF SETUP

preveri_obmocje(samples_number)

# SAVE GRAPHS
print('Saving graphs...')
current_time = datetime.now().strftime("%m-%d-%Y_%H-%M-%S")
# eps
fig_osc.savefig('plots/temp/eps/' + model_type[0] + '_' + str(current_time) + '_osc' + str(samples_number) + '.eps',
                format='eps')
fig_amp.savefig('plots/temp/eps/' + model_type[0] + '_' + str(current_time) + '_amp' + str(samples_number) + '.eps',
                format='eps')
fig_per.savefig('plots/temp/eps/' + model_type[0] + '_' + str(current_time) + '_per' + str(samples_number) + '.eps',
                format='eps')
# png
fig_osc.savefig('plots/temp/png/' + model_type[0] + '_' + str(current_time) + '_osc' + str(samples_number) + '.png',
                format='png')
fig_amp.savefig('plots/temp/png/' + model_type[0] + '_' + str(current_time) + '_amp' + str(samples_number) + '.png',
                format='png')
fig_per.savefig('plots/temp/png/' + model_type[0] + '_' + str(current_time) + '_per' + str(samples_number) + '.png',
                format='png')

# plt.show()
# END OF SAVE GRAPHS

print('Job completed in ', time() - start_time, ' seconds.')
exit(0)
