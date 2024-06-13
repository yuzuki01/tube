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

T = []
Data = []

for fn in files:
    t, data = read_tec(f"./result/{fn}")
    T.append(t)
    X = data[0]
    Data.append(data[1])
Data = np.asarray(Data)

sorted_indices = np.argsort(T)
T.sort()
Data = Data[sorted_indices]

plt.figure(figsize=(15, 5), dpi=100)
plt.pcolormesh(T, X, Data.T, cmap='coolwarm')
plt.colorbar(label='Rho')
plt.xlabel('Time t')
plt.ylabel('X coordinate')
plt.title('Rho(t, X)')
plt.xlim(min(T), max(T))
plt.show()
