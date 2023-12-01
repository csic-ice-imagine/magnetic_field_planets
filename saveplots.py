# saveplots.py contains two functions that either plot the given magnitudes in a plane projections or Mollweide
# projection. Take into account that for the Molleweide projection and the coastline in the plane Earth projection
# you need to have cartopy installed.

import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import cm

plt.rcParams['figure.figsize'] = (10, 7)
plt.rcParams['font.size'] = 16
plt.rcParams['lines.linewidth'] = 1
plt.rcParams["figure.autolayout"] = True
cmap1 = cm.get_cmap('RdBu_r', 255)
cmap2 = cm.get_cmap('inferno', 255)

def planeproj(planet, rc, rc_file, phi, theta, potential, fieldr, fieldtheta, fieldphi, fieldmod, ccrs_library):
    Phi, Theta = np.meshgrid(360 * (1 - phi / 2 / np.pi), - 180 * theta / np.pi + 90)
    names = [r'Potential (Gauss · 1 $R_P$) at $r =$' + str(rc) + '$R_P$', 
             '$B_r$ (Gauss) at $r =$' + str(rc) + '$R_P$',
             '$B_θ$ (Gauss) at $r =$' + str(rc) + '$R_P$', 
             '$B_{\phi}$ (Gauss) at $r =$' + str(rc) + '$R_P$',
             '$|B|$ (Gauss) at $r =$' + str(rc) + '$R_P$']
    files = [planet + '_r_' + rc_file + '_potential.png', 
             planet + '_r_' + rc_file + '_fieldr.png', 
             planet + '_r_' + rc_file + '_fieldtheta.png', 
             planet + '_r_' + rc_file + '_fieldphi.png',
             planet + '_r_' + rc_file + '_fieldmod.png']
    magnitudes = [potential, fieldr, fieldtheta, fieldphi, fieldmod]
    for index, magnitude in enumerate(magnitudes):
        plt.clf()
        if ccrs_library:
            import cartopy.crs as ccrs
            ax = plt.axes(projection=ccrs.PlateCarree())
            if planet=="Earth": ax.coastlines()
        else:
            ax = plt.axes()
        ax.set_xticks([-180,-120,-60,0,60,120,180,240,300,360])
        ax.set_yticks([-90,-60,-30,0,30,60,90])
        ax.set_ylabel("Latitude")
        ax.set_xlabel("Longitude")
        limit = max(np.absolute(np.max(magnitude)), np.absolute(np.min(magnitude)))
        if index==4: 
            realmap = cmap2
            vmin, vmax = np.min(magnitude), np.max(magnitude)
        else:
            realmap = cmap1
            vmin, vmax = -limit, limit
        if planet == "Earth":
            plt.contourf(Phi, Theta, np.flip(magnitude[0, :, :], 1), cmap=realmap, vmin=vmin, vmax=vmax, levels=30)
        else:
            plt.contourf(Phi, Theta, magnitude[0, :, :], cmap=realmap, vmin=vmin, vmax=vmax, levels=30)
        cbar = plt.colorbar(orientation="horizontal", pad=.15, shrink=0.5)
        cbar.set_label(names[index])
        cbar.ax.tick_params(labelsize=11)
        print(planet + "/" + files[index])
        try: plt.savefig(planet + "/" + files[index])
        except:
            os.mkdir(planet)
            plt.savefig(planet + "/" + files[index])

def mollweideproj(planet, rc, rc_file, phi, theta, potential, fieldr, fieldtheta, fieldphi, fieldmod):
    import cartopy.crs as ccrs
    Phi, Theta = np.meshgrid(360 * (1 - phi / 2 / np.pi), - 180 * theta / np.pi + 90)
    names = [r'Potential (Gauss · 1 $R_P$) at $r =$' + str(rc) + '$R_P$', 
             '$B_r$ (Gauss) at $r =$' + str(rc) + '$R_P$',
             '$B_θ$ (Gauss) at $r =$' + str(rc) + '$R_P$', 
             '$B_{\phi}$ (Gauss) at $r =$' + str(rc) + '$R_P$',
             '$|B|$ (Gauss) at $r =$' + str(rc) + '$R_P$']
    files = [planet + '_r_' + rc_file + '_potential_Mollweide.png', 
             planet + '_r_' + rc_file + '_fieldr_Mollweide.png', 
             planet + '_r_' + rc_file + '_fieldtheta_Mollweide.png', 
             planet + '_r_' + rc_file + '_fieldphi_Mollweide.png', 
             planet + '_r_' + rc_file + '_fieldmod_Mollweide.png']
    magnitudes = [potential, fieldr, fieldtheta, fieldphi, fieldmod]
    for index, magnitude in enumerate(magnitudes):
        plt.clf()
        ax = plt.axes(projection=ccrs.Mollweide())
        limit = max(np.absolute(np.max(magnitude)), np.absolute(np.min(magnitude)))
        if index==4: 
            realmap = cmap2
            vmin, vmax = np.min(magnitude), np.max(magnitude)
        else:
            realmap = cmap1
            vmin, vmax = -limit, limit
        plt.contourf(-Phi, Theta, magnitude[0, :, :], cmap=realmap, levels=35, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())
        if planet=="Earth": ax.coastlines()
        cbar = plt.colorbar(orientation="horizontal", pad=.1, shrink=0.5)
        cbar.set_label(names[index])
        cbar.ax.tick_params(labelsize=11)
        print(planet + "/" + files[index])
        try: plt.savefig(planet + "/" + files[index])
        except:
            os.mkdir(planet)
            plt.savefig(planet + "/" + files[index])


def planeproj_1_magnitude(planet, rc, frame, phi, theta, magnitude, index, ccrs_library, year=None, years=False):

    rc_file = '%.3f'%rc
    rc_file = rc_file.replace(".","_")

    Phi, Theta = np.meshgrid(360 * (1 - phi / 2 / np.pi), - 180 * theta / np.pi + 90)
    names = [r'Potential (Gauss · 1 $R_P$) at $r =$' + '%.3f'%rc + '$R_P$', 
             '$B_r$ (Gauss) at $r =$' + '%.3f'%rc + '$R_P$',
             '$B_θ$ (Gauss) at $r =$' + '%.3f'%rc + '$R_P$', 
             '$B_{\phi}$ (Gauss) at $r =$' + '%.3f'%rc + '$R_P$',
             '$|B|$ (Gauss) at $r =$' + '%.3f'%rc + '$R_P$']
    files = [planet + '_potential_r_' + '%03d'%frame + '.png', 
             planet + '_fieldr_r_' + '%03d'%frame + '.png', 
             planet + '_fieldtheta_r_' + '%03d'%frame + '.png', 
             planet + '_fieldphi_r_' + '%03d'%frame + '.png',
             planet + '_fieldmod_r_' + '%03d'%frame + '.png']
    
    plt.clf()
    if ccrs_library:
        import cartopy.crs as ccrs
        ax = plt.axes(projection=ccrs.PlateCarree())
        if planet=="Earth": ax.coastlines()
    else:
        ax = plt.axes()
    ax.set_xticks([-180,-120,-60,0,60,120,180,240,300,360])
    ax.set_yticks([-90,-60,-30,0,30,60,90])
    ax.set_ylabel("Latitude")
    ax.set_xlabel("Longitude")
    limit = max(np.absolute(np.max(magnitude)), np.absolute(np.min(magnitude)))
    if index==4: 
        realmap = cmap2
        vmin, vmax = np.min(magnitude), np.max(magnitude)
    else:
        realmap = cmap1
        vmin, vmax = -limit, limit
    if planet == "Earth":
        plt.contourf(Phi, Theta, np.flip(magnitude), cmap=realmap, vmin=vmin, vmax=vmax, levels=30)
    else:
        plt.contourf(Phi, Theta, magnitude, cmap=realmap, vmin=vmin, vmax=vmax, levels=30)
    cbar = plt.colorbar(orientation="horizontal", pad=.15, shrink=0.5)
    cbar.set_label(names[index])
    cbar.ax.tick_params(labelsize=11)
    if years:
        plt.title(str(year))
        print("Earth_movie_years_r_" + rc_file + "/" + files[index])
    else:
        print(planet + "_movie/" + files[index])
    try: 
        if years:
            plt.savefig("Earth_movie_years_r_" + rc_file + "/" + files[index])
        else:
            plt.savefig(planet + "_movie/" + files[index])
    except:
        if years:
            os.mkdir("Earth_movie_years_r_" + rc_file + "/")
            plt.savefig("Earth_movie_years_r_" + rc_file + "/" + files[index])
        else:
            os.mkdir(planet + "_movie/")
            plt.savefig(planet + "_movie/" + files[index])


def mollweideproj_1_magnitude(planet, rc, frame, phi, theta, magnitude, index, year=None, years=False):
    import cartopy.crs as ccrs

    rc_file = '%.3f'%rc
    rc_file = rc_file.replace(".","_")
    
    Phi, Theta = np.meshgrid(360 * (1 - phi / 2 / np.pi), - 180 * theta / np.pi + 90)
    names = [r'Potential (Gauss · 1 $R_P$) at $r =$' + '%.3f'%rc + '$R_P$', 
             '$B_r$ (Gauss) at $r =$' + '%.3f'%rc + '$R_P$',
             '$B_θ$ (Gauss) at $r =$' + '%.3f'%rc + '$R_P$', 
             '$B_{\phi}$ (Gauss) at $r =$' + '%.3f'%rc + '$R_P$',
             '$|B|$ (Gauss) at $r =$' + '%.3f'%rc + '$R_P$']
    files = [planet + '_potential_Mollweide_r_' + '%03d'%frame + '.png', 
             planet + '_fieldr_Mollweide_r_' + '%03d'%frame + '.png', 
             planet + '_fieldtheta_Mollweide_r_' + '%03d'%frame + '.png', 
             planet + '_fieldphi_Mollweide_r_' + '%03d'%frame + '.png', 
             planet + '_fieldmod_Mollweide_r_' + '%03d'%frame + '.png']
    
    plt.clf()
    ax = plt.axes(projection=ccrs.Mollweide())
    limit = max(np.absolute(np.max(magnitude)), np.absolute(np.min(magnitude)))
    if index==4: 
        realmap = cmap2
        vmin, vmax = np.min(magnitude), np.max(magnitude)
    else:
        realmap = cmap1
        vmin, vmax = -limit, limit
    plt.contourf(-Phi, Theta, magnitude, cmap=realmap, levels=35, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())
    if planet=="Earth": ax.coastlines()
    cbar = plt.colorbar(orientation="horizontal", pad=.1, shrink=0.5)
    cbar.set_label(names[index])
    cbar.ax.tick_params(labelsize=11)
    if years:
        plt.title(str(year))
        print("Earth_movie_years_r_" + rc_file + "/" + files[index])
    else:
        print(planet + "_movie/" + files[index])
    try: 
        if years:
            plt.savefig("Earth_movie_years_r_" + rc_file + "/" + files[index])
        else:
            plt.savefig(planet + "_movie/" + files[index])
    except:
        if years:
            os.mkdir("Earth_movie_years_r_" + rc_file + "/")
            plt.savefig("Earth_movie_years_r_" + rc_file + "/" + files[index])
        else:
            os.mkdir(planet + "_movie/")
            plt.savefig(planet + "_movie/" + files[index])

        