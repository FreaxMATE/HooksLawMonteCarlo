import math
import numpy as np
import matplotlib.pyplot as plt

rng = np.random.default_rng(12345)

class Link:
    """
        dir: 0 -> left, 1 -> right
    """
    def __init__(self, a, beta, dir=None, biased=False, force=0):
        if not dir:
            if biased:
                prob_right = 0.5*(1 + math.tanh(beta*force*a))
                self.dir = 1 if rng.random() < prob_right else 0
            else:
                self.dir = 1 if rng.random() > 0.5 else 0
        else:
            self.dir = dir
    
    def get_dir(self):
        return self.dir

class RubberBand:
    def __init__(self, N, a=1, biased=False, force=0):
        self.N = N
        self.a = a

        # Boltzmann weight
        self.w = 0

        self.kB = 1
        self.T = 1
        self.beta = 1 / (self.kB*self.T)

        self.links = [Link(a=a, beta=self.beta, biased=biased, force=force) for _ in range(N)]

    def get_beta(self):
        return self.beta

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

def histogram_plot(lengths, probs_theory, lengths_mc, ax):
    """
        Plots the histogram and return the chi^2 error
    """

    dl = 1
    edges = [l - dl for l in lengths]
    edges.append(lengths[-1] + dl)

    counts, edges_np = np.histogram(lengths_mc, bins=edges)
    lengths_mc_normed = counts / counts.sum()

    ax[0].bar(lengths, lengths_mc_normed, width=2*dl, edgecolor='black', alpha=0.5, label='n_bands='+str(n_bands)+', Force f = '+str(force))
    ax[0].plot(lengths, probs_theory)
    ax[0].set_ylabel('Pobability')
    ax[0].legend()
    ax[1].plot(lengths, probs_theory/lengths_mc_normed)
    ax[1].axhline(y=1, color='k', linewidth=1)
    ax[1].set_ylabel('Ratio Theory/MC')
    ax[1].set_xlabel('Length L/a')

    return np.sum((lengths_mc_normed - probs_theory)**2/probs_theory)

if __name__ == "__main__":

    fig, ax = plt.subplots(2, 1)

    forces = np.linspace(0.001, 5, 100)

    lengths_mean = np.zeros(len(forces))
    lengths_std = np.zeros(len(forces))

    N = 100
    n_bands = 10000
    a = 1

    for i, force in enumerate(forces):
        bands = [RubberBand(N=N, a=a, biased=True, force=force) for _ in range(n_bands)]

        lengths_mc = [b.length() for b in bands]

        probs_theory, lengths = probability_length_force(N=N, rubberbands=bands, force=force)

        # fig_temp, ax_temp = plt.subplots()
        # ax_temp.plot(lengths_mc)
        # plt.show()

        # lengths_mean[i] = 1/len(lengths_mc) * np.mean(lengths_mc)
        lengths_mean[i] = np.mean(lengths_mc)
        lengths_std[i] = np.std(lengths_mc)

        print('Force: ', force)
        print('Avg length: ', lengths_mean[i])
        print('Probability: ', 0.5*(1 + math.tanh(force)))

        chi_squared = histogram_plot(lengths, probs_theory, lengths_mc, ax)

    plt.savefig('out/histogram_3.png', dpi=300)


    plt.close()
    fig, ax = plt.subplots()

    forces_indiscrete = np.linspace(forces[0], forces[-1], 1000)
    forces_indiscrete_begin = forces_indiscrete[:int(len(forces_indiscrete)/10)]

    ax.errorbar(forces, lengths_mean, fmt='.', yerr=lengths_std, capsize=2, label='MC')

    lf_14 = N*a*np.tanh(1*forces_indiscrete*1)
    lf_hooks = N*a**2*forces_indiscrete_begin / (1*1)

    first_range = int(len(forces)/16)
    forces_fit = forces[:first_range]
    lengths_mean_fit = lengths_mean[:first_range]
    coeff, resid, _, _, _ = np.polyfit(forces_fit, lengths_mean_fit, deg=1, full=True)
    lengths_mean_fit_values = coeff[0]*forces_fit + coeff[1]
    ax.plot(forces_fit, lengths_mean_fit_values, label='Fit')
    print('Fit: ', coeff[0], '* forces + ', coeff[1])
    print(resid)

    ax.plot(forces_indiscrete_begin, lf_hooks, label=r'Small-force approx.')
    ax.plot(forces_indiscrete, lf_14, color='red', label=r'Microscopic probs')
    ax.set_ylabel(r'Average length $\langle L \rangle (f)$')
    ax.set_xlabel('Force $f$')
    ax.legend()
    plt.savefig('out/av_length_by_force3.png', dpi=300)
    plt.show()
