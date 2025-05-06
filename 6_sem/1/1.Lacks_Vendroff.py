import numpy as np
import math
import matplotlib.pyplot as plt

CFL = 1

L = 20
T = 18
h = 0.5

tau = CFL * h
NX = math.ceil(L/h + 1)
Ntau = math.floor(T/tau)+1
u = np.zeros((NX+1, Ntau))
x = np.zeros(NX+1)
t = np.zeros(Ntau)

for i in range(1, len(x)):
    x[i] = (i-1)*h
    u[i, 0] = math.sin(4 * math.pi * x[i] / L)

for time in range(1, Ntau):
    for i in range(2, NX):
        u[i, time] = (CFL**2 + CFL) * u[i-1, time-1] / 2 + (1 - CFL**2) * u[i, time-1] + (CFL**2 - CFL) * u[i+1, time-1]/2
    u[1, time] = (CFL**2 + CFL) * u[NX, time-1] / 2 + (1 - CFL**2) * u[1, time-1] + (CFL**2 - CFL) * u[2, time-1] / 2
    u[NX, time] = (CFL**2 + CFL) * u[NX-1, time-1] / 2 + (1 - CFL**2) * u[NX, time-1] + (CFL**2 - CFL) * u[1, time-1]/2

    if time % 5 == 0:
        print('Time is:', time * tau)
        print(u[:, time])
        plt.plot(x, u[:, time], label=f'T = {time * tau}')

plt.legend()
plt.savefig('LV, Co = 1.01.png')

