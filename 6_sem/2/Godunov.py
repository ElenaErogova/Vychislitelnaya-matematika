import numpy as np
import math
import matplotlib.pyplot as plt

h = 0.05
tau = h/2
T = 2

u_star = 0
F_plus = 0
F_minus = 0
x = np.zeros(math.ceil(2/h)+1)
u = np.zeros((math.ceil(2/h)+1, math.ceil(T/tau+1)))

for i in range(len(x)):
    x[i] = i * h - 1
    u[i, 0] = x[i]

for n in range(1, math.ceil(T/tau)+1):
    u[0, n] = u[0, n-1]
    u[len(x)-1, n] = u[len(x)-1, n-1]
    for i in range(1, len(x)-1):

        if (u[i-1, n] + u[i, n]) / 2 < 0:
            u_star = u[i, n]
        else:
            u_star = u[i-1, n]

        F_minus = u_star**2 / 2

        if (u[i, n] + u[i+1, n]) / 2 < 0:
            u_star = u[i+1, n]
        else:
            u_star = u[i, n]

        F_plus = u_star**2 / 2

        u[i, n] = u[i, n-1] + tau * (F_minus - F_plus) / h

plt.plot(x, u[:, 1], label='Numerical')
plt.plot(x, x/(1+1*0.025), label='Exact')
plt.legend()
plt.savefig('Graph.png')
plt.show()
