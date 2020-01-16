# POZITIVNA POVRATNA ZANKA

import numpy as np
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
m1 = 0
m2 = 0
m3 = 0
K1 = 1
K2 = 1
K3 = 1


def simulate(model, params):
    Z = odeint(model, Z0, t, args=params)
    A = Z[:, 3]
    B = Z[:, 4]
    C = Z[:, 5]
    return A, B, C


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


def get_model_params():
    params = (alpha, alpha0, beta, n, m1, m2, m3, K1, K2, K3)
    return params
