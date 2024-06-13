import os
import re
import numpy as np
from matplotlib import pyplot as plt


def read_tec(filename: str):
    with open(filename, 'r') as f:
        lines = f.readlines()
    time_value = 0.0
    for line in lines:
        match = re.search(r'SOLUTIONTIME=(\d+\.?\d*)', line)
        if match:
            time_value = float(match.group(1))
    return time_value, np.loadtxt(lines[2:], dtype=float).T


files = os.listdir("./result")

solution_time = []
Rho = []
U = []
T = []
q = []

for fn in files:
    t, data = read_tec(f"./result/{fn}")
    solution_time.append(t)
    X = data[0]
    Rho.append(data[1])
    U.append(data[2])
    T.append(data[3])
    q.append(data[4])
Rho = np.asarray(Rho)
U = np.asarray(U)
T = np.asarray(T)
q = np.asarray(q)

solution_time = np.asarray(solution_time)
sorted_indices = np.argsort(solution_time)
solution_time.sort()
Rho = Rho[sorted_indices]
U = U[sorted_indices]
T = T[sorted_indices]
q = q[sorted_indices]

print(Rho.shape)

plt.figure(figsize=(15, 5), dpi=100)
plt.pcolormesh(solution_time, X, Rho.T, cmap='coolwarm')
plt.colorbar(label='Rho')
plt.xlabel('Time t')
plt.ylabel('X coordinate')
plt.title('Rho(t, X)')
plt.xlim(min(solution_time), max(solution_time))
plt.show()
