import random
from datetime import datetime
from sys import exit
from time import time

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
    """
    :param m_number:
    :param K_number:
    :param samples_number:
    :param graph_index:
    :return:
    """

    # params = (alpha, alpha0, beta, n, m1, m2, m3, K1, K2, K3)
    new_params = list(params)

    # zgradimo prazne matrike, kamor bomo shranjevali rezultate
    oscilacije = [[np.nan for x in range(samples_number)] for x in range(4)]
    amplitude = [[np.nan for x in range(samples_number)] for x in range(4)]
    periode = [[np.nan for x in range(samples_number)] for x in range(4)]

    # sampling
    m_samples = [4, 3, 2, 1]
    # K_samples = np.logspace(-1, 3, samples_number)
    K_samples = np.linspace(0.1, 1000, samples_number)

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


def pairplot2(m_number, K_number, samples_number, graph_index):
    # in development
    # params = (alpha, alpha0, beta, n, m1, m2, m3, K1, K2, K3)
    new_params = list(params)

    # zgradimo prazne matrike, kamor bomo shranjevali rezultate
    oscilacije = []
    amplitude = []
    periode = []

    # sampling
    m_samples = [4, 3, 2, 1]
    # K_samples = np.logspace(-1, 3, samples_number)
    K_samples = np.linspace(0.1, 1000, samples_number)

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

def pairplot(samples_number, m_range, K_range, base_amplituda, base_perioda):
    ms = np.linspace(m_range[0], m_range[1], m_range[1]- m_range[0] + 1, dtype=int)
    Ks = np.linspace(10**K_range[0], 10**K_range[1],samples_number)
    oscilating_samples = []
    print (base_amplituda, base_perioda)
    while len(oscilating_samples)<10:
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
            if (amplituda >= base_amplituda and perioda<=base_perioda):
                oscilating_samples.append(list(new_params[4:]))
                print(len(oscilating_samples), amplituda, perioda)

    df = pd.DataFrame(np.array(oscilating_samples), columns=["m1", "m2", "m3", "K1", "K2", "K3"])
    sb.pairplot(df, markers="+")



# util function
def denormalize(normalized_value, min, max, discrete=False):
    if discrete:
        values = np.linspace(min, max, max-min+1, dtype=int)
        splitter = 1/(max-min+1)
        start_splitter = splitter
        counter = 0
        while counter<(max-min):
            if (normalized_value<start_splitter):
                return values[counter]
            counter+=1
            start_splitter +=splitter;
        return values[-1]
    else:
        value = ((max-min) * normalized_value) + min 
        return value

from pyDOE import lhs # pip install pyDOE
def pairplot_lhs(samples_number, m_range, K_range, base_amplituda, base_perioda):
    samples = lhs(6, samples=samples_number, criterion='center')
    oscilating_samples = []
    new_params = list(params)
    better_both_counter=0
    better_amp_counter=0
    better_per_counter=0
    sample_counter = 0

    for sample in samples:

        # denormalize
        new_params[4] = denormalize(sample[0], m_range[0], m_range[1], discrete=True)   # m1
        new_params[5] = denormalize(sample[1], m_range[0], m_range[1], discrete=True)   # m2
        new_params[6] = denormalize(sample[2], m_range[0], m_range[1], discrete=True)   # m3
        new_params[7] = denormalize(sample[3], K_range[0], K_range[1], discrete=False)  # K1
        new_params[8] = denormalize(sample[4], K_range[0], K_range[1], discrete=False)  # K2
        new_params[9] = denormalize(sample[5], K_range[0], K_range[1], discrete=False)  # K3

        A, B, C = simulate(model, tuple(new_params))
        peaks_A, minimums_A, peaks_vals_A, minimums_vals_A = find_peaks(A)

        if oscilating(peaks_vals_A, 0.01):
            amplituda, perioda = eval_signal(peaks_A, minimums_vals_A, peaks_vals_A)
            
            if (amplituda>0 and perioda>0):
                if (amplituda >= base_amplituda and perioda<=base_perioda):
                    oscilating_samples.append(list(new_params[4:]) + ["Boljše oboje"])
                    better_both_counter+=1
                elif (amplituda>base_amplituda):
                    oscilating_samples.append(list(new_params[4:]) + ["Boljša amplituda"])
                    better_amp_counter+=1
                elif (perioda<base_perioda):
                    oscilating_samples.append(list(new_params[4:]) + ["Boljša perioda"])
                    better_per_counter+=1
                
        sample_counter +=1
        print("Progress: " + str(sample_counter/samples_number))

    print ("Base Values:", base_amplituda, base_perioda)
    print("Samples: ", samples_number)
    print("Better amp: ", better_amp_counter, "Better per: ", better_per_counter, "Better both: ", better_both_counter)

    df = pd.DataFrame(oscilating_samples, columns=["m1", "m2", "m3", "K1", "K2", "K3", "Primerjava z osnovo"])
    if better_both_counter>0:
        sb.pairplot(df, markers="+", hue="Primerjava z osnovo", hue_order=["Boljša perioda", "Boljša amplituda", "Boljše oboje"])
    else:
        sb.pairplot(df, markers="+", hue="Primerjava z osnovo", hue_order=["Boljša perioda", "Boljša amplituda"])
    plt.show()


# USER PARAMETERS
model_type = 'positive'
#model_type = 'negative'

figure_size_x = 100
figure_size_y = 15
# figure_size_x = 30
# figure_size_y = 15

samples_number = 100
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

# RUN
#preveri_obmocje(samples_number)
# END OF RUN

############## BENJA
# base model vrednosti
A, B, C = simulate(model, params)
peaks_A, minimums_A, peaks_vals_A, minimums_vals_A = find_peaks(A)
base_amplituda, base_perioda = eval_signal(peaks_A, minimums_vals_A, peaks_vals_A)

#pairplot(10000, [1, 4], [-1, 3], base_amplituda, base_perioda)
pairplot_lhs(200000, [1, 4], [0.1, 1000], base_amplituda, base_perioda)
############## BENJA

# SAVE GRAPHS
print('Saving graphs...')
current_time = datetime.now().strftime("%m-%d-%Y_%H-%M-%S")
# eps
#fig_osc.savefig('plots/temp/eps/' + model_type[0] + '_' + str(current_time) + '_osc' + str(samples_number) + '.eps',
#                format='eps')
#fig_amp.savefig('plots/temp/eps/' + model_type[0] + '_' + str(current_time) + '_amp' + str(samples_number) + '.eps',
#                format='eps')
#fig_per.savefig('plots/temp/eps/' + model_type[0] + '_' + str(current_time) + '_per' + str(samples_number) + '.eps',
#                format='eps')
# png
#fig_osc.savefig('plots/temp/png/' + model_type[0] + '_' + str(current_time) + '_osc' + str(samples_number) + '.png',
#                format='png')
#fig_amp.savefig('plots/temp/png/' + model_type[0] + '_' + str(current_time) + '_amp' + str(samples_number) + '.png',
#                format='png')
#fig_per.savefig('plots/temp/png/' + model_type[0] + '_' + str(current_time) + '_per' + str(samples_number) + '.png',
#                format='png')

# plt.show()
# END OF SAVE GRAPHS

print('Job completed in ', time() - start_time, ' seconds.')
exit(0)
