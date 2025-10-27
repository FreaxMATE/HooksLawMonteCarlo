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

def histogram_plot(lengths, probs_theory, lengths_mc, n_bands):
    """
        Plots the histogram and return the chi^2 error
    """
    fig, ax = plt.subplots(2, 1, height_ratios=[2, 1], sharex=True, figsize=(12, 8))

    dl = 1
    edges = [l - dl for l in lengths]
    edges.append(lengths[-1] + dl)

    counts, edges_np = np.histogram(lengths_mc, bins=edges)
    lengths_mc_normed = counts / counts.sum()

    print(counts)

    ax[0].step(lengths, lengths_mc_normed, where='mid', label='Length for n_bands='+str(n_bands))
    ax[0].errorbar(lengths, lengths_mc_normed, fmt=',', color='lawngreen', yerr=np.sqrt(counts)/np.sum(counts), capsize=3, label=r'$\sigma = \sqrt{counts}$')
    ax[0].plot(lengths, probs_theory, label='Theory')
    ax[0].set_ylabel('Probability')
    ax[0].legend()
    ax[1].step(lengths, lengths_mc_normed/probs_theory, where='mid')
    ax[1].axhline(y=1, color='k', linewidth=1)
    ax[1].set_ylabel('Ratio MC/Theory')
    ax[1].set_xlabel('Length L/a')
    ax[1].set_xlim(-70, 70)

    return np.sum((counts - n_bands*probs_theory)**2/(n_bands*probs_theory))/(len(lengths)-1)

if __name__ == "__main__":

    N = 100

    probs_theory, lengths = probability_length(N=N)

    n_bands = 10000

    bands = [RubberBand(N=N) for _ in range(n_bands)]
    lengths_mc = [b.length() for b in bands]

    chi_squared = histogram_plot(lengths, probs_theory, lengths_mc, n_bands)
    print(chi_squared)

    plt.savefig('out/histogram_v1.png', dpi=300)
    plt.show()
