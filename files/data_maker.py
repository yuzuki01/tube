import os
import re
import pickle
import numpy as np
from matplotlib import pyplot as plt


def read_tec(filename: str):
    with open(filename, 'r') as f:
        lines = f.readlines()
    return np.loadtxt(lines[2:], dtype=float).T


data_path = [os.path.join("./", item) for item in os.listdir("./") if os.path.isdir(item) and item.startswith('kn')]
data_set = []
for folder in data_path:
    config_file = os.path.join(folder, "config.txt")
    step_files = [os.path.join(folder, file) for file in os.listdir(folder) if re.match(r"step-[\d.e-]+.dat.plt", file)]
    case_data = []
    for step_file in step_files:
        step_data = read_tec(step_file)
        case_data.append(step_data)
    data_set.append(case_data)

data_set = np.asarray(data_set)

# 转化为 (X, t, Kn, 7)
# 7 个量为表示 "X", "time", "Kn", "Rho", "U", "T", "q"
data_set = np.transpose(data_set, (3, 1, 0, 2))


# 特征处理
def process_kn(kn):
    return (np.tanh(0.8 * (np.log10(kn) + 2.5)) + 1) / 2


data_set_processed = np.copy(data_set)
kn_data = data_set_processed[:, :, :, 2]
kn_data_processed = process_kn(kn_data)
data_set_processed[:, :, :, 2] = kn_data_processed

print("Data-set shape: ", data_set_processed.shape)
print("Get flattened data.")

X_num, t_num, Kn_num, *args = data_set_processed.shape
data_set_processed = data_set_processed.reshape((-1, 7))
inputs = data_set_processed[:, :3]
outputs = data_set_processed[:, 3:]

data_num = X_num * t_num * Kn_num
print(f"Predicted shape should be: ({data_num}, 3) and ({data_num}, 4)")
print("Inputs shape:", inputs.shape)  # 应该是 (X_num * t_num * Kn_num, 3)
print("Outputs shape:", outputs.shape)  # 应该是 (X_num * t_num * Kn_num, 4)

with open("data.pkl", "wb") as file:
    pickle.dump({
        "input_names": ["X", "t", "f(Kn)"],
        "output_names": ["Rho", "U", "T", "q"],
        "inputs": inputs,
        "outputs": outputs,
    }, file)
    print("Save to pickle file.")
