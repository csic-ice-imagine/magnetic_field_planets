# lowes_spec.py has the function that calculates the Lowes spectrum at a given radius, and also has two other functions which
# plot one or several spectra depending on the options given in main.py

import numpy as np
import matplotlib.pyplot as plt

def lowes_spec(NPOL,    # Number of spherical harmonics used for the magnetic reconstruction
               rc,      # Value of the radius for the Lowes spectrum to be calculated
               g,       # Values for the g constants in the Schmidt polynomial
               h):      # Values for the h constants in the Schmidt polynomial
    Rn = np.zeros(NPOL)
    for n in range(0,NPOL):
        for m in range(0,n):
            Rn[n] = Rn[n] + (g[n,m]**2+h[n,m]**2)
        Rn[n] = rc**(-2*n-4)*(n+1)*Rn[n]
    return Rn

def plot_lowes(planet,             # String used for saving plots in the right folder name
               rc,                 # Float used for the x-label plotting
               rc_file,            # String used name for radius the file name (in case of movie it is the frame)
               Rn,                 # Precalculated Lewis spectrum
               movie=False,        # Switch to save the plot into the correct movie folder
               year=None,          # Year for the year Earth movie plots (used only on plot title)
               years=False):       # Switch to save the plot into the correct movie folder for the Earth yearly movie
    plt.clf()
    plt.grid(ls='--', alpha=.5)
    plt.scatter(np.arange(1, len(Rn)), Rn[1:])
    plt.xlabel('Harmonic degree (n) at $r =$' + '%.3f'%rc + '$R_P$')
    plt.ylabel(r'R$_n$ (nT)$^2$')
    plt.yscale('log')
    if movie:
        if years:
            plt.title(str(year))
            rc_file_true = '%.3f'%rc
            rc_file_true = rc_file_true.replace(".","_")
            plt.savefig(planet + "_movie_years_r_" + rc_file_true + "/" + planet + '_lowes_spectrum_r_' + '%03d'%rc_file + '.png')
            print(planet + "_movie_years_r_" + rc_file_true + "/" + planet + '_lowes_spectrum_r_' + '%03d'%rc_file + '.png')
        else:
            plt.savefig(planet + "_movie/" + planet + '_lowes_spectrum_r_' + '%03d'%rc_file + '.png')
            print(planet + "_movie/" + planet + '_lowes_spectrum_r_' + '%03d'%rc_file + '.png')
    else: 
        plt.savefig(planet + "/" + planet + '_lowes_spectrum_r_' + rc_file + '.png')
        print(planet + "/" + planet + '_lowes_spectrum_r_' + rc_file + '.png')

def plot_multiple_lowes(planet,        # String used for saving plots in the right folder name
                        lowes_radii,   # Set of radii for which to plot the Lowes spectra on the same plot 
                        Rn):           # All different precaluculated and joint in a single array
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