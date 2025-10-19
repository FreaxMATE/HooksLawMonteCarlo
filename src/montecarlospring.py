import math
import numpy as np
import matplotlib.pyplot as plt

rng = np.random.default_rng(12345)

class Link:
    """
        dir: 0 -> left, 1 -> right
    """
    def __init__(self, dir=None):
        if not dir:
            self.dir = 1 if rng.random() > 0.5 else 0
        else:
            self.dir = dir
    
    def get_dir(self):
        return self.dir

class RubberBand:
    def __init__(self, N, a=1):
        self.N = N
        self.a = a
        self.links = [Link() for _ in range(N)]

    def length(self):
        return self.a*(2*np.sum([l.dir for l in self.links]) - self.N)


def probability_length(N, a=1):
    """
        Calculates the probability of the chain having the length L.
    """
    prob = np.zeros(N)
    lengths = np.zeros(N)
    for n in range(N):
        L = a*(2*n - N)
        # omega_Nn = np.exp(N*np.log(2) - L**2/(2*N*a**2))
        omega_Nn = math.comb(N, n)
        prob[n] = omega_Nn/2**N
        lengths[n] = L
    return prob, lengths

def histogram_plot(lengths, probs_theory, lengths_mc):
    """
        Plots the histogram and return the chi^2 error
    """
    fig, ax = plt.subplots(2, 1)

    dl = 1
    edges = [l - dl for l in lengths]
    edges.append(lengths[-1] + dl)

    counts, edges_np = np.histogram(lengths_mc, bins=edges)
    lengths_mc_normed = counts / counts.sum()

    ax[0].bar(lengths, lengths_mc_normed, width=2*dl, color='tab:orange', edgecolor='black', alpha=0.5, label='Length for n_bands='+str(n_bands))
    ax[0].plot(lengths, probs_theory)
    ax[0].set_ylabel('Pobability')
    ax[0].legend()
    ax[1].plot(lengths, probs_theory/lengths_mc_normed)
    ax[1].axhline(y=1, color='k', linewidth=1)
    ax[1].set_ylabel('Ratio Theory/MC')
    ax[1].set_xlabel('Length L/a')

    return np.sum((lengths_mc_normed - probs_theory)**2/probs_theory)

if __name__ == "__main__":

    N = 100

    probs_theory, lengths = probability_length(N=N)

    n_bands = 10000

    bands = [RubberBand(N=N) for _ in range(n_bands)]
    lengths_mc = [b.length() for b in bands]

    chi_squared = histogram_plot(lengths, probs_theory, lengths_mc)
    print(chi_squared)

    plt.savefig('out/histogram_v1.png', dpi=300)
    plt.show()
