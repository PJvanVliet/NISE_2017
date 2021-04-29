import sys
import numpy as np 
import matplotlib.pyplot as plt 

def plot(filename):
    """
    Reads the data files and labels, plots each of them.
    """
    data = np.loadtxt(filename)
    time = data[:, 0]

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    lineshapes = ["--", "-", "-."]
    labels = ["NISE", "Prezhdo", "Alternative"]
    for i in range(data.shape[1] - 1):
        pop = data[:, i+1]
        ax.plot(time, pop, lineshapes[i], label=labels[i])
    ax.set_xlabel(r"$t$ (fs)")
    ax.set_xlim(0, time[-1])
    ax.set_ylabel(r"$\bar{\rho}_{00}(t)$")
    plt.legend()
    plt.grid()
    plt.savefig(f"{filename[:-4]}.png")
    plt.close()

fn = sys.argv[1]
plot(fn)
