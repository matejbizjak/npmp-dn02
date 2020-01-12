from datetime import datetime
from sys import exit

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb

from utils import eval_signal, oscilating, find_peaks


def preveri_obmocje(samples_number):
    # vse kombinacije m-jev in K-jev
    graph_index = 1
    for i in range(1, 4):
        for j in range(1, 4):
            print(graph_index, 'out of ', 9)
            sampling(i, j, samples_number, graph_index)
            graph_index += 1


def sampling(m_number, K_number, samples_number, graph_index):
    # params = (alpha, alpha0, beta, n, m1, m2, m3, K1, K2, K3)
    new_params = list(params)

    # zgradimo prazne matrike, kamor bomo shranjevali rezultate
    oscilacije = [[np.nan for x in range(samples_number)] for x in range(4)]
    amplitude = [[np.nan for x in range(samples_number)] for x in range(4)]
    periode = [[np.nan for x in range(samples_number)] for x in range(4)]

    # sampling
    m_samples = [4, 3, 2, 1]
    K_samples = np.logspace(-1, 3, samples_number)

    # za vse sample pozenemo simulacijo in rezultate shranimo v prej narejene matrike
    for x in range(len(m_samples)):
        print('\t', x + 1, 'out of ', len(m_samples))
        for y in range(len(K_samples)):
            new_params[3 + m_number] = m_samples[x]
            new_params[6 + K_number] = K_samples[y]

            A, B, C = simulate(model, tuple(new_params))

            peaks_A, minimums_A, peaks_vals_A, minimums_vals_A = find_peaks(A)

            if oscilating(peaks_vals_A, 0.01):
                oscilacije[x][y] = 1

                amplituda, perioda = eval_signal(peaks_A, minimums_vals_A, peaks_vals_A)
                amplitude[x][y] = amplituda
                periode[x][y] = perioda

    # narisemo grafe
    plot_graph(oscilacije, m_number, K_number, m_samples, K_samples, graph_index, fig_osc, 'Oscilacije')
    plot_graph(amplitude, m_number, K_number, m_samples, K_samples, graph_index, fig_amp, 'Amplitude')
    plot_graph(periode, m_number, K_number, m_samples, K_samples, graph_index, fig_per, 'Periode')


def plot_graph(data, m_number, K_number, m_samples, K_samples, graph_index, figure, figure_title):
    data_frame = pd.DataFrame(data, columns=list(map(lambda k: round(k, 3), K_samples)), index=m_samples)
    subplot = figure.add_subplot(3, 3, graph_index)
    figure.suptitle(figure_title, fontsize=30)
    sb.heatmap(data_frame, ax=subplot, cbar=True, square=False, mask=data_frame.isnull(), cmap="YlGnBu")
    var_list = ["alpha", "alpha0", "beta", "n", "m1", "m2", "m3", "K1", "K2", "K3"]
    subplot.set_xlabel(var_list[6 + K_number])
    subplot.set_ylabel(var_list[3 + m_number])


# SETTINGS
# model_type = 'positive'
model_type = 'negative'

# figure_size_x = 100
# figure_size_y = 10
figure_size_x = 30
figure_size_y = 15

samples_number = 100
# END OF SETTINGS

_temp = __import__(model_type + '_core', globals(), locals(), ['get_model_params', 'model', 'simulate'], 0)
get_model_params = _temp.get_model_params
model = _temp.model
simulate = _temp.simulate

params = get_model_params()

fig_osc = plt.figure(figsize=(figure_size_x, figure_size_y))
fig_amp = plt.figure(figsize=(figure_size_x, figure_size_y))
fig_per = plt.figure(figsize=(figure_size_x, figure_size_y))

preveri_obmocje(samples_number)

# SHRANJEVANJE GRAFOV
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

exit(0)
