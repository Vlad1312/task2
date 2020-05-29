import matplotlib.pyplot as plt

import numpy as np
import os
from numpy import linspace

import requests as rqst
import re
import scipy.special as sp



url = 'http://www.jenyay.net/uploads/Student/Modelling/task_02.txt'
task = rqst.get(url)
z = re.search(r'^6\..+', task.text, flags=re.M)
z1 = (z.group().split(';'))
D = float(z1[0].split('=')[1])
fmin = float(z1[1].split('=')[1])
fmax = float(z1[2].split('=')[1])


def task2():
 os.mkdir("results")
 res = open('results/task_02_4O-506C_Pripechenkov_6.txt', 'w')
 print('f[ГГц]\t    sigma[м^2]', file=res)

 c = 3e8
 r = D / 2
 f = np.linspace(fmin, fmax + 1, 500)
 sigma1 = []
 for fx in f:
    lamda = c / fx
    k = 2 * np.pi / lamda
    h = []
    b = [0]
    a = []
    sigmaN = []
    n = 0
    while n < 70:
        h.append(
            sp.spherical_jn(n, k * r) + 1j * sp.spherical_yn(n, k * r))
        a.append(sp.spherical_jn(n, k * r) / h[n])
        n += 1
    n = 1
    while n < 70:
        b.append((k * r * sp.spherical_jn(n - 1, k * r) - n *
                  sp.spherical_jn(n, k * r)) / (k * r * h[n - 1] - n * h[n]))
        n += 1
    n = 1
    while n < 70:
        sigmaN.append(((-1) ** n) * (n + 1 / 2) * (b[n] - a[n]))
        n += 1
    sigma = (lamda ** 2 / np.pi) * (abs(np.sum(sigmaN))) ** 2
    print(str(fx * 1e-9) + '\t ' + str(sigma), file=res)
    sigma1.append(sigma / (np.pi * (r ** 2)))


 fig, ax = plt.subplots()

 ax.plot (2 * np.pi * f * r / c, sigma1)

 plt.show()

if os.path.exists('results'):
     os.remove('results/task_02_4O-506C_Pripechenkov_6.txt')
     os.rmdir("results")
     task2()
else:
     task2()