# POZITIVNA POVRATNA ZANKA

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.signal as scipy
import seaborn as sns
from scipy.integrate import odeint

# SETUP
# initial condition
Z0 = [7.02028482, 1.57288111, 58.50876295, 8.69333201, 1.40155783, 53.49915358]

# number of time points
t_end = 100
n = t_end * 10

# time points
t = np.linspace(0, t_end, n)

# parameters
alpha = 216
alpha0 = 0.001 * alpha
n = 2
beta = 5
m1 = 1
m2 = 1
m3 = 1
K1 = 4
K2 = 4
K3 = 4

params = (alpha, alpha0, beta, n, m1, m2, m3, K1, K2, K3)


# function that returns dZ/dt
def model(Z, t, alpha, alpha0, beta, n, m1, m2, m3, K1, K2, K3):
    a = Z[0]
    b = Z[1]
    c = Z[2]
    A = Z[3]
    B = Z[4]
    C = Z[5]

    # mRNAs
    dadt = (alpha * ((A / K1) ** m1)) / (1 + C ** n + (A / K1) ** m1) + alpha0 - a
    dbdt = (alpha * ((B / K2) ** m2)) / (1 + A ** n + (B / K2) ** m2) + alpha0 - b
    dcdt = (alpha * ((C / K3) ** m3)) / (1 + B ** n + (C / K3) ** m3) + alpha0 - c

    # proteins
    dAdt = - beta * (A - a)
    dBdt = - beta * (B - b)
    dCdt = - beta * (C - c)

    dzdt = [dadt, dbdt, dcdt, dAdt, dBdt, dCdt]

    return dzdt


# simulation
def simulate(model, Z0, t, params):
    Z = odeint(model, Z0, t, args=params)
    A = Z[:, 3]
    B = Z[:, 4]
    C = Z[:, 5]
    return A, B, C


# peaks and valleys detection
def findPeaks(simResults):
    peaks = scipy.find_peaks(simResults)[0]  # polozaji
    minimums = scipy.find_peaks(list(map(lambda x: x * (-1), simResults)))[0]  # polozaji
    peaks_vals = list(map(lambda x: simResults[x], peaks))  # vrednosti
    minimums_vals = list(map(lambda x: simResults[x], minimums))  # vrednosti

    return peaks, minimums, peaks_vals, minimums_vals


# check if oscilating
def oscilating(peakVals, tolerance):
    try:
        return abs(peakVals[-2] - peakVals[-3]) < tolerance
    except:
        return False


# find amplitude and periode
def evalSignal(peaks, minVals, peakVals):
    amplitude = (peakVals[-2] - minVals[-2]) / 2
    periode = (peaks[-2] - peaks[-3]) / 10  # 10 = n/t_end
    return amplitude, periode


# 2d analiza

def preveri_obmocje(samples_number):
    graph_index = 1
    for i in range(1, 4):
        for j in range(1, 4):
            print(graph_index, 'out of ', 9)
            sampling(i, j, samples_number, graph_index)
            graph_index += 1


def sampling(m_number, K_number, samples_number, graph_index):
    # params = (alpha, alpha0, beta, n, m1, m2, m3, K1, K2, K3)
    newParams = list(params)

    oscilacije = [[np.nan for x in range(samples_number)] for x in range(4)]
    amplitude = [[np.nan for x in range(samples_number)] for x in range(4)]
    periode = [[np.nan for x in range(samples_number)] for x in range(4)]

    m_samples = [4, 3, 2, 1]
    K_samples = np.logspace(-1, 3, samples_number)

    for x in range(len(m_samples)):
        print('\t', x)
        for y in range(len(K_samples)):
            newParams[3 + m_number] = m_samples[x]
            newParams[6 + K_number] = K_samples[y]

            A, B, C = simulate(model, Z0, t, tuple(newParams))

            peaks_A, minimums_A, peaks_vals_A, minimums_vals_A = findPeaks(A)

            if oscilating(peaks_vals_A, 0.01):
                oscilacije[x][y] = 1

                amplituda, perioda = evalSignal(peaks_A, minimums_vals_A, peaks_vals_A)
                amplitude[x][y] = amplituda
                periode[x][y] = perioda

    # graf za oscilacije
    dataFrame1 = pd.DataFrame(oscilacije, columns=list(map(lambda k: round(k, 3), K_samples)), index=m_samples)
    varList = ["alpha", "alpha0", "beta", "n", "m1", "m2", "m3", "K1", "K2", "K3"]
    ax1 = fig1.add_subplot(3, 3, graph_index)
    fig1.suptitle('Oscilacije', fontsize=30)
    sns.heatmap(dataFrame1, ax=ax1, cbar=True, square=False, mask=dataFrame1.isnull(), cmap="YlGnBu")
    ax1.set_xlabel(varList[6 + K_number])
    ax1.set_ylabel(varList[3 + m_number])

    # graf za amplitude
    dataFrame2 = pd.DataFrame(amplitude, columns=list(map(lambda k: round(k, 3), K_samples)), index=m_samples)
    varList = ["alpha", "alpha0", "beta", "n", "m1", "m2", "m3", "K1", "K2", "K3"]
    ax2 = fig2.add_subplot(3, 3, graph_index)
    fig2.suptitle('Amplitude', fontsize=30)
    sns.heatmap(dataFrame2, ax=ax2, cbar=True, square=False, mask=dataFrame2.isnull(), cmap="YlGnBu")
    ax2.set_xlabel(varList[6 + K_number])
    ax2.set_ylabel(varList[3 + m_number])

    # graf za periode
    dataFrame3 = pd.DataFrame(periode, columns=list(map(lambda k: round(k, 3), K_samples)), index=m_samples)
    varList = ["alpha", "alpha0", "beta", "n", "m1", "m2", "m3", "K1", "K2", "K3"]
    ax3 = fig3.add_subplot(3, 3, graph_index)
    fig3.suptitle('Periode', fontsize=30)
    sns.heatmap(dataFrame3, ax=ax3, cbar=True, square=False, mask=dataFrame3.isnull(), cmap="YlGnBu")
    ax3.set_xlabel(varList[6 + K_number])
    ax3.set_ylabel(varList[3 + m_number])


# ANALIZA

# output figure
# fig = plt.figure(figsize=(100, 10))
# fig2 = plt.figure(figsize=(100, 10))
# fig3 = plt.figure(figsize=(100, 10))

fig1 = plt.figure(figsize=(25, 10))
fig2 = plt.figure(figsize=(25, 10))
fig3 = plt.figure(figsize=(25, 10))

preveri_obmocje(50)

# SHRANJEVANJE GRAFOV
# plt.savefig('testplot.eps', format='eps')

plt.show()
