import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb

from negative_core import get_model_params, model, simulate
# from positive_core import get_model_params, model, simulate
from utils import eval_signal, oscilating, find_peaks

from genetics import genetic_algorithm


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
        print('\t', x)
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

    # graf za oscilacije
    data_frame1 = pd.DataFrame(oscilacije, columns=list(map(lambda k: round(k, 3), K_samples)), index=m_samples)
    var_list = ["alpha", "alpha0", "beta", "n", "m1", "m2", "m3", "K1", "K2", "K3"]
    ax1 = fig1.add_subplot(3, 3, graph_index)
    fig1.suptitle('Oscilacije', fontsize=30)
    sb.heatmap(data_frame1, ax=ax1, cbar=True, square=False, mask=data_frame1.isnull(), cmap="YlGnBu")
    ax1.set_xlabel(var_list[6 + K_number])
    ax1.set_ylabel(var_list[3 + m_number])

    # graf za amplitude
    data_frame2 = pd.DataFrame(amplitude, columns=list(map(lambda k: round(k, 3), K_samples)), index=m_samples)
    var_list = ["alpha", "alpha0", "beta", "n", "m1", "m2", "m3", "K1", "K2", "K3"]
    ax2 = fig2.add_subplot(3, 3, graph_index)
    fig2.suptitle('Amplitude', fontsize=30)
    sb.heatmap(data_frame2, ax=ax2, cbar=True, square=False, mask=data_frame2.isnull(), cmap="YlGnBu")
    ax2.set_xlabel(var_list[6 + K_number])
    ax2.set_ylabel(var_list[3 + m_number])

    # graf za periode
    data_frame3 = pd.DataFrame(periode, columns=list(map(lambda k: round(k, 3), K_samples)), index=m_samples)
    var_list = ["alpha", "alpha0", "beta", "n", "m1", "m2", "m3", "K1", "K2", "K3"]
    ax3 = fig3.add_subplot(3, 3, graph_index)
    fig3.suptitle('Periode', fontsize=30)
    sb.heatmap(data_frame3, ax=ax3, cbar=True, square=False, mask=data_frame3.isnull(), cmap="YlGnBu")
    ax3.set_xlabel(var_list[6 + K_number])
    ax3.set_ylabel(var_list[3 + m_number])


params = get_model_params()
# velikost koncne slike v px*10 (x, y)
# fig = plt.figure(figsize=(100, 10))
# fig2 = plt.figure(figsize=(100, 10))
# fig3 = plt.figure(figsize=(100, 10))

fig1 = plt.figure(figsize=(25, 10))
fig2 = plt.figure(figsize=(25, 10))
fig3 = plt.figure(figsize=(25, 10))

limiting_params = genetic_algorithm(mode="all_params")
# limiting_params = genetic_algorithm()
print(limiting_params)

preveri_obmocje(50)
# ce zelis spremeniti tip modela (poz. v neg. povratno zanko samo spremeni importe)

# SHRANJEVANJE GRAFOV
# plt.savefig('testplot.eps', format='eps')

plt.show()
