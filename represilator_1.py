# POZITIVNA POVRATNA ZANKA

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.signal as scipy
import seaborn as sns
from scipy.integrate import odeint


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


# initial condition
# Z0 = [10,0,0,0,0,0]
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


# output figure
fig = plt.figure(figsize=(8, 4))  # 3 grafi - 8,8, sicer 8,4


# # 2d analiza
# def preveriObmocje(var1, var2, range1, range2, n1=100,
#                    n2=100):  # var1 in var2 - indeks spremenljivk (alpha, alpha0, beta, n, m1, m2, m3, K1, K2, K3);
#     newParams = list(params)  # range1, range2 - array [spodnja, zgornja];
#     amplitudes = []  # n1, n2 - število vzorcev
#     periodes = []
#     i = 0
#     for x in range(range1[0], range1[1] + 1, int(abs(range1[1] + 1 - range1[0]) / n1)):
#         amplitudes.append([])
#         periodes.append([])
#         for y in range(range2[0], range2[1] + 1, int(abs(range2[1] + 1 - range2[0]) / n2)):
#             newParams[var1] = x
#             newParams[var2] = y
#
#             A, B, C = simulate(model, Z0, t, tuple(newParams))
#
#             peaks_A, minimums_A, peaks_vals_A, minimums_vals_A = findPeaks(A)
#             peaks_B, minimums_B, peaks_vals_B, minimums_vals_B = findPeaks(B)
#             peaks_C, minimums_C, peaks_vals_C, minimums_vals_C = findPeaks(C)
#
#             if oscilating(peaks_vals_A, 0.01) and oscilating(peaks_vals_A, 0.01) and oscilating(peaks_vals_A,
#                                                                                                 0.01):  # ali moramo preveriti vse tri?
#                 amplitude, periode = evalSignal(peaks_A, minimums_vals_A, peaks_vals_A)
#                 # je dovolj če samo A?
#                 if (amplitude < 0):
#                     amplitude = 0
#
#                 if (periode < 0):
#                     periode = 0
#                 # recimo da samo amplitudo gledamo
#             else:
#                 amplitude = 0
#                 periode = 0
#
#             amplitudes[i].append(amplitude)
#             periodes[i].append(periode)
#         i = i + 1
#
#     dataFrame1 = pd.DataFrame(amplitudes)
#     dataFrame2 = pd.DataFrame(periodes)
#     # print(dataFrame)
#
#     varList = ["alpha", "alpha0", "beta", "n", "m1", "m2", "m3", "K1", "K2", "K3"]
#     ax1 = fig.add_subplot(121)  # če imamo 3 grafe  - 221, sicer 121
#     ax2 = fig.add_subplot(122)  # če imamo 3 grafe  - 222, sicer 122
#     ax1.set_title("Amplitude")
#     ax2.set_title("Periode")
#     sns.heatmap(dataFrame1, xticklabels=list(range(range1[0], range1[1] + 1, int(abs(range1[1] + 1 - range1[0]) / n1))),
#                 yticklabels=list(range(range2[0], range2[1] + 1, int(abs(range2[1] + 1 - range2[0]) / n2))), ax=ax1)
#     ax1.set_xlabel(varList[var2])
#     ax1.set_ylabel(varList[var1])
#     sns.heatmap(dataFrame2, xticklabels=list(range(range1[0], range1[1] + 1, int(abs(range1[1] + 1 - range1[0]) / n1))),
#                 yticklabels=list(range(range2[0], range2[1] + 1, int(abs(range2[1] + 1 - range2[0]) / n2))), ax=ax2)
#     ax2.set_xlabel(varList[var2])
#     ax2.set_ylabel(varList[var1])


# analiza 2.0 (m1=m2=m3 in K1=K2=K3 ampak gledamo vse kombinacije - 2 for zanki, časovno zahtevno)
def preveriObmocjeFiksnoKomb(range1, range2, n1=100,
                             n2=100):
    # indeks spremenljivk (alpha, alpha0, beta, n, m1, m2, m3, K1, K2, K3);
    newParams = list(params)  # range1, range2 - array [spodnja, zgornja];
    amplitudes = []  # n1, n2 - število vzorcev
    periodes = []
    i = 0
    for x in range(range1[0], range1[1] + 1, int(abs(range1[1] + 1 - range1[0]) / n1)):
        amplitudes.append([])
        periodes.append([])
        print("Iteracija: ", x)
        for y in range(range2[0], range2[1] + 1, int(abs(range2[1] + 1 - range2[0]) / n2)):
            newParams[4] = x  # m1
            newParams[5] = x  # m2
            newParams[6] = x  # m3
            newParams[7] = y  # K1
            newParams[8] = y  # K2
            newParams[9] = y  # K3

            A, B, C = simulate(model, Z0, t, tuple(newParams))

            peaks_A, minimums_A, peaks_vals_A, minimums_vals_A = findPeaks(A)
            peaks_B, minimums_B, peaks_vals_B, minimums_vals_B = findPeaks(B)
            peaks_C, minimums_C, peaks_vals_C, minimums_vals_C = findPeaks(C)

            if oscilating(peaks_vals_A, 0.01) and oscilating(peaks_vals_A, 0.01) and oscilating(peaks_vals_A,
                                                                                                0.01):  # ali moramo
                # preveriti vse tri? Dobim podobne rezultate
                amplitude, periode = evalSignal(peaks_A, minimums_vals_A, peaks_vals_A)
                # je dovolj če samo A?
                if (amplitude < 0):
                    amplitude = 0

                if (periode < 0):
                    periode = 0
                # recimo da samo amplitudo gledamo
            else:
                amplitude = 0
                periode = 0

            amplitudes[i].append(amplitude)
            periodes[i].append(periode)
        i = i + 1

    dataFrame1 = pd.DataFrame(amplitudes)
    dataFrame2 = pd.DataFrame(periodes)

    # print(dataFrame)

    varList = ["alpha", "alpha0", "beta", "n", "m1", "m2", "m3", "K1", "K2", "K3"]
    ax1 = fig.add_subplot(121)  # če imamo 3 grafe  - 221, sicer 121
    ax2 = fig.add_subplot(122)  # če imamo 3 grafe  - 222, sicer 122
    ax1.set_title("Amplitude")
    ax2.set_title("Periode")
    sns.heatmap(dataFrame1, xticklabels=list(range(range1[0], range1[1] + 1, int(abs(range1[1] + 1 - range1[0]) / n1))),
                yticklabels=list(range(range2[0], range2[1] + 1, int(abs(range2[1] + 1 - range2[0]) / n2))), ax=ax1)
    ax1.set_xlabel("K")
    ax1.set_ylabel("m")
    sns.heatmap(dataFrame2, xticklabels=list(range(range1[0], range1[1] + 1, int(abs(range1[1] + 1 - range1[0]) / n1))),
                yticklabels=list(range(range2[0], range2[1] + 1, int(abs(range2[1] + 1 - range2[0]) / n2))), ax=ax2)
    ax2.set_xlabel("K")
    ax2.set_ylabel("m")


# analiza 3.0 (m1=m2=m3=K1=K2=K3 le 1 for zanka)
def preveriObmocjeFiksno(range1, n1=100):
    # indeks spremenljivk (alpha, alpha0, beta, n, m1, m2, m3, K1, K2, K3);
    newParams = list(params)  # range1, range2 - array [spodnja, zgornja];
    amplitudes = []  # n1, n2 - število vzorcev
    periodes = []

    for x in range(range1[0], range1[1] + 1, int(abs(range1[1] + 1 - range1[0]) / n1)):
        newParams[4] = x  # m1
        newParams[5] = x  # m2
        newParams[6] = x  # m3
        newParams[7] = x  # K1
        newParams[8] = x  # K2
        newParams[9] = x  # K3

        A, B, C = simulate(model, Z0, t, tuple(newParams))

        peaks_A, minimums_A, peaks_vals_A, minimums_vals_A = findPeaks(A)
        peaks_B, minimums_B, peaks_vals_B, minimums_vals_B = findPeaks(B)
        peaks_C, minimums_C, peaks_vals_C, minimums_vals_C = findPeaks(C)

        if oscilating(peaks_vals_A, 0.01) and oscilating(peaks_vals_A, 0.01) and oscilating(peaks_vals_A,
                                                                                            0.01):  # ali moramo
            # preveriti vse tri?
            amplitude, periode = evalSignal(peaks_A, minimums_vals_A, peaks_vals_A)
            # je dovolj če samo A?
            if (amplitude < 0):
                amplitude = 0

            if (periode < 0):
                periode = 0
            # recimo da samo amplitudo gledamo
        else:
            amplitude = 0
            periode = 0

        amplitudes.append(amplitude)
        periodes.append(periode)

    dataFrame1 = pd.DataFrame(amplitudes)
    dataFrame2 = pd.DataFrame(periodes)

    # print(dataFrame)

    varList = ["alpha", "alpha0", "beta", "n", "m1", "m2", "m3", "K1", "K2", "K3"]
    ax1 = fig.add_subplot(121)  # če imamo 3 grafe  - 221, sicer 121
    ax2 = fig.add_subplot(122)  # če imamo 3 grafe  - 222, sicer 122
    ax1.set_title("Amplitude")
    ax2.set_title("Periode")

    ax1.plot(amplitudes)
    ax1.set_xlabel("K")
    ax1.set_ylabel("m")

    ax2.plot(periodes)
    ax2.set_xlabel("K")
    ax2.set_ylabel("m")


# plot results

def plotResults(results):
    A = results[0]
    peaks_A, minimums_A, peaks_vals_A, minimums_vals_A = findPeaks(A)
    B = results[1]
    C = results[2]
    ax = fig.add_subplot(212)
    ax.plot(t, A, 'g:', label='A(t)')
    ax.plot(t, B, 'b-', label='B(t)')
    ax.plot(t, C, 'r--', label='C(t)')
    ax.set_ylabel('Koncentracije')
    ax.set_xlabel('Čas')
    ax.legend(loc='right')


##################
##################
##################

# ANALIZA

# parametri: alpha, alpha0, beta, n, m1, m2, m3, K1, K2, K3
# indeksi:     0       1      2   3   4   5   6   7   8   9
# preveriObmocje(4, 7, [1, 4], [1, 4], 4, 4)
preveriObmocjeFiksnoKomb([1, 1000], [1, 1000], 1000, 1000)
# preveriObmocjeFiksno([1, 75], 75)

# PRIKAZ GRAFA PRI DOLOČENIH PARAMETRIH

# m1=3
# K1=1
# params= (alpha, alpha0, beta, n, m1, m2, m3, K1, K2, K3)
# plotResults(simulate(model, Z0, t, params))


# SHRANJEVANJE GRAFOV
# plt.savefig('testplot.eps', format='eps')


plt.show()
