import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint


# function that returns dZ/dt
def model(Z, t, alpha, alpha0, beta, n):
    a = Z[0]
    b = Z[1]
    c = Z[2]
    A = Z[3]
    B = Z[4]
    C = Z[5]

    # mRNAs
    dadt = alpha / (1 + C ** n) + alpha0 - a
    dbdt = alpha / (1 + A ** n) + alpha0 - b
    dcdt = alpha / (1 + B ** n) + alpha0 - c

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
t_end = 25
n = t_end * 10

# time points
t = np.linspace(0, t_end, n)

alpha = 216
alpha0 = 0.001 * alpha
n = 2
beta = 5

params = (alpha, alpha0, beta, n)

Z = odeint(model, Z0, t, args=params)
A = Z[:, 3]
B = Z[:, 4]
C = Z[:, 5]

# plot results
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(t, A, 'g:', label='A(t)')
ax.plot(t, B, 'b-', label='B(t)')
ax.plot(t, C, 'r--', label='C(t)')
ax.set_ylabel('concentrations')
ax.set_xlabel('time')
ax.legend(loc='best')
plt.show()

##################
##################
#################
