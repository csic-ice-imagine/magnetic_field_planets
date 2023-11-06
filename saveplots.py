import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import cm

cmap1 = cm.get_cmap('RdBu_r', 255)
cmap2 = cm.get_cmap('inferno', 255)

def planeproj(planet, rc, rc_file, phi, theta, Nr, Ntheta, Nphi, potential, fieldr, fieldtheta, fieldphi, fieldmod):
    Phi, Theta = np.meshgrid(360 * (1 - phi / 2 / np.pi), - 180 * theta / np.pi + 90)
    names = [r'Potential (Gauss · 1 $R_P$) at $r =$' + str(rc) + '$R_P$', '$B_r$ (Gauss) at $r =$' + str(rc) + '$R_P$',
             '$B_θ$ (Gauss) at $r =$' + str(rc) + '$R_P$', '$B_{\phi}$ (Gauss) at $r =$' + str(rc) + '$R_P$',
             '$|B|$ (Gauss) at $r =$' + str(rc) + '$R_P$']
    files = [planet + '_r_' + rc_file + '_potential.png', planet + '_r_' + rc_file + '_fieldr.png', planet + '_r_' + rc_file + '_fieldtheta.png', planet + '_r_' + rc_file + '_fieldphi.png',
             planet + '_r_' + rc_file + '_fieldmod.png']
    magntidues = [potential, fieldr, fieldtheta, fieldphi, fieldmod]
    printvalues = np.zeros([Nr, Ntheta, Nphi])
    for number in range(0, 5):
        plt.clf()
        plt.xticks([0, 60, 120, 180, 240, 300, 360])
        plt.yticks([-90, -60, -30, 0, 30, 60, 90])
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        printvalues = magntidues[number]
        limit = max(np.absolute(np.max(printvalues)), np.absolute(np.min(printvalues)))
        if number==4: 
            realmap = cmap2
            vmin, vmax = np.min(printvalues), np.max(printvalues)
        else:
            realmap = cmap1
            vmin, vmax = -limit, limit
        plt.contourf(Phi, Theta, printvalues[0, :, :], cmap=realmap, vmin=vmin, vmax=vmax, levels=30)
        plt.gca().invert_xaxis()
        cbar = plt.colorbar(orientation="horizontal", pad=.15, shrink=0.5)
        cbar.set_label(names[number])
        cbar.ax.tick_params(labelsize=11)
        print("Saving " + planet + "/" + files[number])
        try: plt.savefig(planet + "/" + files[number])
        except:
            os.mkdir(planet)
            plt.savefig(planet + "/" + files[number])

def mollweideproj(planet, rc, rc_file, phi, theta, Nr, Ntheta, Nphi, potential, fieldr, fieldtheta, fieldphi, fieldmod):
    import cartopy.crs as ccrs
    Phi, Theta = np.meshgrid(360 * (1 - phi / 2 / np.pi), - 180 * theta / np.pi + 90)
    names = [r'Potential (Gauss · 1 $R_P$) at $r =$' + str(rc) + '$R_P$', '$B_r$ (Gauss) at $r =$' + str(rc) + '$R_P$',
             '$B_θ$ (Gauss) at $r =$' + str(rc) + '$R_P$', '$B_{\phi}$ (Gauss) at $r =$' + str(rc) + '$R_P$',
             '$|B|$ (Gauss) at $r =$' + str(rc) + '$R_P$']
    files = [planet + '_r_' + rc_file + '_potential_Mollweide', planet + '_r_' + rc_file + '_fieldr_Mollweide', planet + '_r_' + rc_file + '_fieldtheta_Mollweide', planet + '_r_' + rc_file + '_fieldphi_Mollweide', planet + '_r_' + rc_file + '_fieldmod_Mollweide']
    magntidues = [potential, fieldr, fieldtheta, fieldphi, fieldmod]
    printvalues = np.zeros([Nr, Ntheta, Nphi])
    for number in range(0, 5):
        plt.clf()
        ax = plt.axes(projection=ccrs.Mollweide())
        printvalues = magntidues[number]
        limit = max(np.absolute(np.max(printvalues)), np.absolute(np.min(printvalues)))
        if number==4: 
            realmap = cmap2
            vmin, vmax = np.min(printvalues), np.max(printvalues)
        else:
            realmap = cmap1
            vmin, vmax = -limit, limit
        plt.contourf(-Phi, Theta, printvalues[0, :, :], cmap=realmap, levels=35, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())
        if planet=="Earth": ax.coastlines()
        cbar = plt.colorbar(orientation="horizontal", pad=.1, shrink=0.5)
        cbar.set_label(names[number])
        cbar.ax.tick_params(labelsize=11)
        print("Saving " + planet + "/" + files[number])
        try: plt.savefig(planet + "/" + files[number])
        except:
            os.mkdir(planet)
            plt.savefig(planet + "/" + files[number])