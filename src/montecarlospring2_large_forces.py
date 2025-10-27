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

        # Boltzmann weight
        self.w = 0

        self.kB = 1
        self.T = 1
        self.beta = 1 / (self.kB*self.T)

        self.links = [Link() for _ in range(N)]

    def get_beta(self):
        return self.beta

    def boltzmann_weight(self, force):
        self.w = np.exp(self.beta*self.length()*force)
        return self.w

    def length(self):
        return self.a*(2*np.sum([l.dir for l in self.links]) - self.N)


def probability_length_force(N, rubberbands, a=1, force=0.0):
    """
        Calculates the probability of the chain having the length L.
    """
    prob = np.zeros(N)
    lengths = np.zeros(N)
    prob_l_numerator = np.zeros(N)

    Z = 0
    for n in range(N):
        lengths[n] = a*(2*n - N)
        omega_Nn = math.comb(N, n)
        prob_l_numerator[n] = omega_Nn*np.exp(rubberbands[n].get_beta()*force*lengths[n])
        Z += prob_l_numerator[n]

    for n in range(N):
        prob[n] = prob_l_numerator[n]/Z

    return prob, lengths

def weighted_histogram_plot(lengths, probs_theory, lengths_mc, boltzmann_weights, n_bands, force, color_index):
    """
        Plots the histogram and return the chi^2 error
    """

    colors=['tab:blue', 'tab:orange', 'tab:green']
    colors_dark=['steelblue', 'darkorange', 'darkgreen']

    dl = 1
    edges = [l - dl for l in lengths]
    edges.append(lengths[-1] + dl)

    print('Lengths: ', np.shape(lengths_mc))
    print('Boltzmann weights: ', np.shape(boltzmann_weights))

    counts, edges_np = np.histogram(lengths_mc, bins=edges, weights=boltzmann_weights)
    lengths_mc_normed = counts / counts.sum()

    ax[0].step(lengths, lengths_mc_normed, where='mid', color=colors[color_index], label='n_bands='+str(n_bands)+', Force f = '+str(force))
    ax[0].plot(lengths, probs_theory, color=colors_dark[color_index], label='Theory force f = '+str(force))
    ax[0].set_ylabel('Probability')
    ax[0].legend()

    print('Force')
    print('  Theory: ', probs_theory)
    print('  Length: ', lengths_mc_normed)
    print('  Ratio: ', lengths_mc_normed/probs_theory)
    ax[1].step(lengths, lengths_mc_normed/probs_theory, where='mid', color=colors[color_index])
    ax[1].axhline(y=1, color='k', linewidth=1)
    ax[1].set_ylabel('Ratio MC/Theory')
    ax[1].set_xlabel('Length L/a')
    # ax[1].set_xlim(-70, 70)
    ax[1].set_ylim(-0.5, 3)
    ax[0].set_ylim(-0.1, 0.2)

    return np.sum((counts - n_bands*probs_theory)**2/(n_bands*probs_theory))/(len(lengths)-1)

if __name__ == "__main__":

    fig, ax = plt.subplots(2, 1, height_ratios=[2, 1], sharex=True, figsize=(12, 8))

    N = 100
    n_bands = 10000

    bands = [RubberBand(N=N) for _ in range(n_bands)]


    for color_index, force in enumerate([0.1, 0.5, 1]):
        lengths_mc = [b.length() for b in bands]
        boltzmann_weights = [b.boltzmann_weight(force) for b in bands]

        probs_theory, lengths = probability_length_force(N=N, rubberbands=bands, force=force)

        chi_squared = weighted_histogram_plot(lengths, probs_theory, lengths_mc, boltzmann_weights, n_bands, force, color_index)
        print(chi_squared)

    plt.savefig('out/weighted_histogram_2_large_forces.png', dpi=300)
    plt.show()
