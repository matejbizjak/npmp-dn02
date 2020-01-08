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
fig = plt.figure(figsize=(25, 25))  # 3 grafi - 8,8, sicer 8,4


# 2d analiza

# uniformno vzorčenje, dve spremenljivki
def preveriObmocje(var1, var2, range1, range2, n1=100,
                   n2=100):  # var1 in var2 - indeks spremenljivk (alpha, alpha0, beta, n, m1, m2, m3, K1, K2, K3);
    newParams = list(params)  # range1, range2 - array [spodnja, zgornja];
    amplitudes = []  # n1, n2 - število vzorcev
    periodes = []
    i = 0
    for x in range(range1[0], range1[1] + 1, int(abs(range1[1] + 1 - range1[0]) / n1)):
        amplitudes.append([])
        periodes.append([])
        for y in range(range2[0], range2[1] + 1, int(abs(range2[1] + 1 - range2[0]) / n2)):
            newParams[var1] = x
            newParams[var2] = y

            A, B, C = simulate(model, Z0, t, tuple(newParams))

            peaks_A, minimums_A, peaks_vals_A, minimums_vals_A = findPeaks(A)
            peaks_B, minimums_B, peaks_vals_B, minimums_vals_B = findPeaks(B)
            peaks_C, minimums_C, peaks_vals_C, minimums_vals_C = findPeaks(C)

            if oscilating(peaks_vals_A, 0.01) and oscilating(peaks_vals_A, 0.01) and oscilating(peaks_vals_A,
                                                                                                0.01):  # ali moramo preveriti vse tri?
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
    ax1.set_xlabel(varList[var2])
    ax1.set_ylabel(varList[var1])
    sns.heatmap(dataFrame2, xticklabels=list(range(range1[0], range1[1] + 1, int(abs(range1[1] + 1 - range1[0]) / n1))),
                yticklabels=list(range(range2[0], range2[1] + 1, int(abs(range2[1] + 1 - range2[0]) / n2))), ax=ax2)
    ax2.set_xlabel(varList[var2])
    ax2.set_ylabel(varList[var1])


# heatmap K1=K2=K3 / m1=m2=m3
def heatmap_mK(range1, range2, n1=100, n2=100):  # 1 - m, 2 - K
    newParams = list(params)
    amplitudes = []
    periodes = []
    i = 0
    x = range1[0]
    xs_labels = []
    ys_labels = []

    while x <= range1[1]:
        amplitudes.append([])
        periodes.append([])
        xs_labels.append(x)

        y = range2[0]
        while y <= range2[1]:
            if len(ys_labels) < n2:
                ys_labels.append(y)
            newParams[4] = x
            newParams[5] = x
            newParams[6] = x
            newParams[7] = y
            newParams[8] = y
            newParams[9] = y

            A, B, C = simulate(model, Z0, t, tuple(newParams))

            peaks_A, minimums_A, peaks_vals_A, minimums_vals_A = findPeaks(A)
            peaks_B, minimums_B, peaks_vals_B, minimums_vals_B = findPeaks(B)
            peaks_C, minimums_C, peaks_vals_C, minimums_vals_C = findPeaks(C)

            if oscilating(peaks_vals_A, 0.01) and oscilating(peaks_vals_A, 0.01) and oscilating(peaks_vals_A,
                                                                                                0.01):  # ali moramo preveriti vse tri?
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
            y += (abs(range2[1] - range2[0]) / n2)

        i = i + 1
        x += abs(range1[1] - range1[0]) / n1

    dataFrame1 = pd.DataFrame(amplitudes)  # , columns=[1,2,3,4,5,6], index=["a","b","a","b","a","b"])
    dataFrame2 = pd.DataFrame(periodes)
    ax1 = fig.add_subplot(121)  # če imamo 3 grafe  - 221, sicer 121
    ax2 = fig.add_subplot(122)  # če imamo 3 grafe  - 222, sicer 122
    ax1.set_title("Amplitude")
    ax2.set_title("Periode")
    sns.heatmap(dataFrame1, xticklabels=xs_labels, yticklabels=ys_labels, ax=ax1)
    ax1.set_xlabel("K")
    ax1.set_ylabel("m")
    sns.heatmap(dataFrame2, xticklabels=xs_labels, yticklabels=ys_labels, ax=ax2)
    ax2.set_xlabel("K")
    ax2.set_ylabel("m")


# preveri območje, implementirano s sudoku_lhs - število razdelitev območij (n1 in n2) mora bit enako za obe spremenljivki in deljivo s 3
# ns - število uporabljenih vzorcev - tudi mora biti deljivo s 3
def preveriObmocje_LHS(var1, var2, range1, range2, ns=300, n1=999, n2=999):
    if (n1 != n2):
        print("n1 not equal n2")
        return 0

    newParams = list(params)
    amplitudes = [[0 for x in range(n1)] for x in range(n1)]
    periodes = [[0 for x in range(n1)] for x in range(n1)]
    xs = np.linspace(range1[0], range1[1], n1)
    ys = np.linspace(range2[0], range2[1], n2)

    # sudoku lhs
    import sudoku_lhs
    S, _ = sudoku_lhs.sudoku.sample(2, 3, int(ns / 9), visualize=False, verbose=False)
    # S - array arrayev, stolpci - prvi element, vrstica - drugi element

    for [x, y] in S:

        val_x = xs[x]
        val_y = ys[y]

        newParams[var1] = val_x
        newParams[var2] = val_y

        A, B, C = simulate(model, Z0, t, tuple(newParams))

        peaks_A, minimums_A, peaks_vals_A, minimums_vals_A = findPeaks(A)
        peaks_B, minimums_B, peaks_vals_B, minimums_vals_B = findPeaks(B)
        peaks_C, minimums_C, peaks_vals_C, minimums_vals_C = findPeaks(C)

        if oscilating(peaks_vals_A, 0.01) and oscilating(peaks_vals_A, 0.01) and oscilating(peaks_vals_A,
                                                                                            0.01):  # ali moramo preveriti vse tri?
            amplitude, periode = evalSignal(peaks_A, minimums_vals_A, peaks_vals_A)
            # je dovolj če samo A?
            if (amplitude < 0):
                amplitude = 0

            if (periode < 0):
                periode = 0

        else:
            amplitude = 0
            periode = 0

        amplitudes[x][y] = amplitude
        periodes[x][y] = periode

    dataFrame1 = pd.DataFrame(amplitudes, columns=xs, index=ys)
    dataFrame2 = pd.DataFrame(periodes, columns=xs, index=ys)
    # print(dataFrame)

    varList = ["alpha", "alpha0", "beta", "n", "m1", "m2", "m3", "K1", "K2", "K3"]
    ax1 = fig.add_subplot(121)  # če imamo 3 grafe  - 221, sicer 121
    ax2 = fig.add_subplot(122)  # če imamo 3 grafe  - 222, sicer 122
    ax1.set_title("Amplitude")
    ax2.set_title("Periode")
    sns.heatmap(dataFrame1,
                ax=ax1)  # , xticklabels=list(range(range1[0], range1[1] + 1, int(abs(range1[1] + 1 -range1[0])/n1))), yticklabels=list(range(range2[0], range2[1] + 1, int(abs(range2[1] +1 -range2[0])/n2))), ax=ax1)
    ax1.set_xlabel(varList[var2])
    ax1.set_ylabel(varList[var1])
    sns.heatmap(dataFrame2,
                ax=ax2)  # xticklabels=list(range(range1[0], range1[1] + 1, int(abs(range1[1] + 1 -range1[0])/n1))), yticklabels=list(range(range2[0], range2[1] + 1, int(abs(range2[1] +1 -range2[0])/n2))), ax=ax2)
    ax2.set_xlabel(varList[var2])
    ax2.set_ylabel(varList[var1])


# isto kot preveriObmocje_LHS, le da dobimo graf, ki prikezuje le ali oscilira, ali ne oscilira
# preveri območje, implementirano s sudoku_lhs - število razdelitev območij (n1 in n2) mora bit enako za obe spremenljivki in deljivo s 3
# ns - število uporabljenih vzorcev - tudi mora biti deljivo s 3
# multi(Bool) - pove, ali spreminjamo vse K-je in m-je
def preveriOscilacije_LHS(var1, var2, range1, range2, current_index, ns=300, n1=999, n2=999, multi=False):
    if (n1 != n2):
        print("n1 not equal n2")
        return 0

    newParams = list(params)
    oscilacije = [[1 for x in range(n1)] for x in range(n1)]
    xs = np.logspace(range1[0], range1[1], n1)
    ys = np.logspace(range2[0], range2[1], n2)

    # sudoku lhs
    import sudoku_lhs
    S, _ = sudoku_lhs.sudoku.sample(2, 3, int(ns / 9), visualize=False, verbose=False)
    # S - array arrayev, stolpci - prvi element, vrstica - drugi element

    for [x, y] in S:

        val_x = xs[x]
        val_y = ys[y]

        if multi:
            newParams[4] = val_x
            newParams[5] = val_x
            newParams[6] = val_x
            newParams[7] = val_y
            newParams[8] = val_y
            newParams[9] = val_y
        else:
            newParams[var1] = val_x
            newParams[var2] = val_y

        A, B, C = simulate(model, Z0, t, tuple(newParams))

        peaks_A, minimums_A, peaks_vals_A, minimums_vals_A = findPeaks(A)
        # peaks_B, minimums_B, peaks_vals_B, minimums_vals_B = findPeaks(B)
        # peaks_C, minimums_C, peaks_vals_C, minimums_vals_C = findPeaks(C)

        if oscilating(peaks_vals_A, 0.01):
            oscilira = -0.4
            amplitude, periode = evalSignal(peaks_A, minimums_vals_A, peaks_vals_A)
            if (amplitude < 0):
                oscilira = -1

            if (periode < 0):
                oscilira = -1

        else:
            oscilira = -1

        oscilacije[x][y] = oscilira

        dataFrame1 = pd.DataFrame(oscilacije, columns=list(map(lambda x: round(x, 3), xs)),
                                  index=list(map(lambda y: round(y, 3), ys)))

    # print(dataFrame)

    varList = ["alpha", "alpha0", "beta", "n", "m1", "m2", "m3", "K1", "K2", "K3"]
    ax1 = fig.add_subplot(6, 6, current_index)  # če imamo 3 grafe  - 221, sicer 121
    # ax2=fig.add_subplot(122)   #če imamo 3 grafe  - 222, sicer 122
    # ax1.set_title("Oscilacije")
    # ax2.set_title("Periode")

    if current_index == 31:
        sns.heatmap(dataFrame1, ax=ax1, cmap="Pastel1",
                    cbar=False)
        ax1.set_xlabel(varList[var2])
        ax1.set_ylabel(varList[var1])
    elif current_index % 6 == 1:
        sns.heatmap(dataFrame1, ax=ax1, cmap="Pastel1",
                    cbar=False, xticklabels=0)
        ax1.set_ylabel(varList[var1])
    else:
        if current_index >= 31:
            sns.heatmap(dataFrame1, ax=ax1, cmap="Pastel1",
                        cbar=False, yticklabels=0)
            ax1.set_xlabel(varList[var2])
        else:
            sns.heatmap(dataFrame1, ax=ax1, cmap="Pastel1",
                        cbar=False, yticklabels=0,
                        xticklabels=0)

    # ax1.set_xlabel(varList[var2])
    # ax1.set_ylabel(varList[var1])
    # sns.heatmap(dataFrame2, ax=ax2)#xticklabels=list(range(range1[0], range1[1] + 1, int(abs(range1[1] + 1 -range1[0])/n1))), yticklabels=list(range(range2[0], range2[1] + 1, int(abs(range2[1] +1 -range2[0])/n2))), ax=ax2)
    # ax2.set_xlabel(varList[var2])
    # ax2.set_ylabel(varList[var1])


# plot results


##################
##################
##################

# ANALIZA

# parametri: alpha, alpha0, beta, n, m1, m2, m3, K1, K2, K3
# indeksi:     0       1      2   3   4   5   6   7   8   9
# preveriObmocje_LHS(4,7, [1,4], [5,6], 30 , 30, 30)

current_index = 1
for i in range(4, 10):
    for j in range(4, 10):
        preveriOscilacije_LHS(i, j, [-1, 4], [-1, 4], current_index, 30, 30, 30, multi=True)
        print(current_index)
        current_index += 1

# heatmap_mK([10**(-3),1],[10**(-3),1], 5, 5)

# PRIKAZ GRAFA PRI DOLOČENIH PARAMETRIH

m1 = 0.9
K1 = 0.001
m2 = 0.9
K2 = 0.001
m3 = 0.9
K3 = 0.001

params = (alpha, alpha0, beta, n, m1, m2, m3, K1, K2, K3)
# plotResults(simulate(model, Z0, t, params))


# SHRANJEVANJE GRAFOV
# plt.savefig('testplot.eps', format='eps')


plt.show()
