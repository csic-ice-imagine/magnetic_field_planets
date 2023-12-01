# lowes_spec.py has the function that calculates the Lowes spectrum at a given radius, and also has two other functions which
# plot one or several spectra depending on the options given in main.py

import numpy as np
import matplotlib.pyplot as plt

def lowes_spec(NPOL, rc, g, h):
    Rn = np.zeros(NPOL)
    for n in range(0,NPOL):
        for m in range(0,n):
            Rn[n] = Rn[n] + (g[n,m]**2+h[n,m]**2)
        Rn[n] = rc**(-2*n-4)*(n+1)*Rn[n]
    return Rn

def plot_lowes(planet, rc, rc_file, Rn, movie=False):
    plt.clf()
    plt.grid(ls='--', alpha=.5)
    plt.scatter(np.arange(1, len(Rn)), Rn[1:])
    plt.xlabel('Harmonic degree (n) at $r =$' + '%.3f'%rc + '$R_P$')
    plt.ylabel(r'R$_n$ (nT)$^2$')
    plt.yscale('log')
    if movie: 
        plt.savefig(planet + "_movie/" + planet + '_lowes_spectrum_r_' + '%03d'%rc_file + '.png')
        print(planet + "_movie/" + planet + '_lowes_spectrum_r_' + '%03d'%rc_file + '.png')
    else: 
        plt.savefig(planet + "/" + planet + '_lowes_spectrum_r_' + rc_file + '.png')
        print(planet + "/" + planet + '_lowes_spectrum_r_' + rc_file + '.png')

def plot_multiple_lowes(planet, lowes_radii, Rn):
    plt.clf()
    colors = plt.cm.jet(np.linspace(0.10,0.90, len(lowes_radii)))
    plt.grid(ls='--', alpha=.5)
    for r in range(0, len(lowes_radii)):
        plt.plot(np.arange(1, np.shape(Rn)[1]), Rn[r,1:], marker='o', color=colors[r], label=r'r = ' + f'{lowes_radii[r]:.2f}' + r'R$_P$')
    plt.xlabel('Harmonic degree (n)')
    plt.ylabel(r'R$_n$ (nT)$^2$')
    plt.legend(loc='best')
    plt.yscale('log')
    plt.savefig(planet + "/" + planet + '_lowes_spectrum.png')
    print(planet + "/" + planet + '_lowes_spectrum.png')