# Version by Daniele with more plots at the same time, and with generic fields

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import csv

# Spherical grid numbers
Ntheta = 100
Nphi = Ntheta*2
Nr = 3
a = 1  # Should be 6371.2/72492 (Earth/Jupiter), but renormalize to 1, since r/a is what matters
rc = 1  # Radius considered in the map plot (CMB = 0.455/0.85 for earth/jupiter).
dr = 0.0001
const = 1.  # To fo from nanotesla to gauss
planet = "Jupiter"  # Can be Earth, Jupiter
year = 2020  # Can be 1900, 1905, 1910, ..., to 2020. Only in the case of the Earth

printall = True  # Saves all 2D plots

plot2d = False  # Shows a given plot (for testing)
Plot = 0
# if 0 you will plot the divergence
# if 1 the curlmod
# if 2 the curvaturemod
# if 3 the curlr
# if 4 the curltheta
# if 5 the curlphi
# if 6 the curvaturer
# if 7 the curvaturetheta
# if 8 the curvaturephi

# Definition of the spherical grid
phi = np.linspace(1 * np.pi / Nphi, 2 * np.pi *(1 - 0.5/Nphi), num=Nphi)
theta = np.linspace(.5 * np.pi / Ntheta, np.pi * (1 - .5 / Ntheta), num=Ntheta)
radius = np.linspace(rc * (1. - dr), rc * (1. + dr), num=Nr)
potential = np.zeros([Nr, Ntheta, Nphi])
# Initialize all components of the magnetic field (shperical and modulus)
fieldr = np.zeros([Nr, Ntheta, Nphi])
fieldtheta = np.zeros([Nr, Ntheta, Nphi])
fieldphi = np.zeros([Nr, Ntheta, Nphi])
fieldmod = np.zeros([Nr, Ntheta, Nphi])

# This part reads the g's and h's of the spherical expansion of the magnetic field, depending on the planet
if planet == "Earth":
    NPOL = 14
    g = np.zeros([NPOL, NPOL])
    h = np.zeros([NPOL, NPOL])
    index = int((year - 1900) / 5 + 3)
    file = open("igrf13coeffs.txt", "r")
    lines = file.readlines()[4:]
    for n in range(0, len(lines)):
        file_list = [i for i in lines[n].split()]
        if file_list[0] == 'g':
            g[int(file_list[1]), int(file_list[2])] = float(file_list[index])
        else:
            h[int(file_list[1]), int(file_list[2])] = float(file_list[index])
elif planet == "Jupiter":
    NPOL = 11
    g = np.zeros([NPOL, NPOL])
    h = np.zeros([NPOL, NPOL])
    file = open("grl57087-sup-0005-2018GL077312-ds01.txt", "r")
    lines = file.readlines()[1:]
    for n in range(0, len(lines)):
        file_list = [i for i in lines[n].split()]
        if file_list[3] == 'g':
            g[int(file_list[4]), int(file_list[5])] = float(file_list[1])
        else:
            h[int(file_list[4]), int(file_list[5])] = float(file_list[1])
else:
    NPOL = 0

#NPOL=10
#g[:,:] = 0.
#h[:,:] = 0.
#g[1,1] = 1.
#h[1,1] = 0.

# This part defines all K's and S's used in the recursive formulas for
# the Schmidt quasi-normalized associated Legendre polynomials,
# which are orthogonal, but not normal.
K = np.zeros([NPOL, NPOL])
S = np.zeros([NPOL, NPOL])
S[0, 0] = 1
delta = 0

for n in range(1, NPOL):
    for m in range(0, n + 1):
        K[n, m] = ((n - 1) ** 2 - m ** 2) / (2 * n - 1) / (2 * n - 3)
        if n == 1:
            K[n, m] = 0
        if m == 0:
            S[n, 0] = (2 * n - 1) * S[n - 1, 0] / n
        else:
            if m == 1:
                delta = 1
            else:
                delta = 0
            S[n, m] = S[n, m - 1] * np.sqrt((n - m + 1) * (delta + 1) / (n + m))

# This part defines the Gaussian-normalized and
# the Schmidt quasi-normalized associated Legendre polynomials
Pgauss = np.zeros([NPOL, NPOL, Ntheta])
derivPgauss = np.zeros([NPOL, NPOL, Ntheta])
P = np.zeros([NPOL, NPOL, Ntheta])
derivP = np.zeros([NPOL, NPOL, Ntheta])

# DV: I have simplified here a lot (one could put also K and S here somehow)
# (not sure if it would be faster to have theta inside, instead of outside)

for ntheta in range(0, Ntheta):
    Pgauss[0, 0, ntheta] = 1
    Pgauss[1, 0, ntheta] = np.cos(theta[ntheta])
    derivPgauss[0, 0, ntheta] = 0
    derivPgauss[1, 0, ntheta] = - np.sin(theta[ntheta])
    for n in range(1, NPOL):
        Pgauss[n, n, ntheta] = np.sin(theta[ntheta]) * Pgauss[n - 1, n - 1, ntheta]
        derivPgauss[n, n, ntheta] = np.sin(theta[ntheta]) * derivPgauss[n - 1, n - 1, ntheta] + np.cos(
            theta[ntheta]) * Pgauss[n - 1, n - 1, ntheta]
        for m in range(0, n):
            Pgauss[n, m, ntheta] = np.cos(theta[ntheta]) * Pgauss[n - 1, m, ntheta] - K[n, m] * Pgauss[
                n - 2, m, ntheta]
            derivPgauss[n, m, ntheta] = np.cos(theta[ntheta]) * derivPgauss[n - 1, m, ntheta] - np.sin(
                theta[ntheta]) * Pgauss[n - 1, m, ntheta] - K[n, m] * derivPgauss[n - 2, m, ntheta]
    P[:, :, ntheta] = S[:, :] * Pgauss[:, :, ntheta]
    derivP[:, :, ntheta] = S[:, :] * derivPgauss[:, :, ntheta]


# This function calculates the potential, fieldr, fieldtheta and fieldphi for a given r, theta and phi
def function(radius0, ntheta, phi0):
    poten, fr, ftheta, fphi = 0, 0, 0, 0
    for n in range(1, NPOL):
        suma, sumatheta, sumaphi = 0, 0, 0
        for m in range(0, n + 1):
            suma += P[n, m, ntheta] * (g[n, m] * np.cos(m * phi0) + h[n, m] * np.sin(m * phi0))
            sumatheta += derivP[n, m, ntheta] * (g[n, m] * np.cos(m * phi0) + h[n, m] * np.sin(m * phi0))
            sumaphi += m * P[n, m, ntheta] * (- g[n, m] * np.sin(m * phi0) + h[n, m] * np.cos(m * phi0))
        poten += suma * (1 / radius0) ** (n + 1)
        fr += suma * (n + 1) * (1 / radius0) ** (n + 2)
        ftheta += - sumatheta * (1 / radius0) ** (n + 2)
        fphi += - sumaphi * (1 / radius0) ** (n + 2) / np.sin(theta[ntheta])
    return poten / const, fr / const, ftheta / const, fphi / const


# Loops for all r, theta and phi and obtaining all the corresponding potential and fields
for i in range(0, Nr):
    print('Nr =', i)
    for j in range(0, Ntheta):
        for k in range(0, Nphi):
            potential[i, j, k], fieldr[i, j, k], fieldtheta[i, j, k], fieldphi[i, j, k] = function(radius[i], j, phi[k])

# Defines all the direction vector magnetic field
fieldmod = np.sqrt(fieldr ** 2 + fieldtheta ** 2 + fieldphi ** 2)
directionr = fieldr / fieldmod
directionphi = fieldphi / fieldmod
directiontheta = fieldtheta / fieldmod

# Initializes all possible derivatives
derivrfieldr, derivrdirectionr = np.zeros([Ntheta, Nphi]), np.zeros([Ntheta, Nphi])
derivrfieldtheta, derivrdirectiontheta = np.zeros([Ntheta, Nphi]), np.zeros([Ntheta, Nphi])
derivrfieldphi, derivrdirectionphi = np.zeros([Ntheta, Nphi]), np.zeros([Ntheta, Nphi])

derivthetafieldr, derivthetadirectionr = np.zeros([Ntheta, Nphi]), np.zeros([Ntheta, Nphi])
derivthetafieldtheta, derivthetadirectiontheta = np.zeros([Ntheta, Nphi]), np.zeros([Ntheta, Nphi])
derivthetafieldphi, derivthetadirectionphi = np.zeros([Ntheta, Nphi]), np.zeros([Ntheta, Nphi])

derivphifieldr, derivphidirectionr = np.zeros([Ntheta, Nphi]), np.zeros([Ntheta, Nphi])
derivphifieldtheta, derivphidirectiontheta = np.zeros([Ntheta, Nphi]), np.zeros([Ntheta, Nphi])
derivphifieldphi, derivphidirectionphi = np.zeros([Ntheta, Nphi]), np.zeros([Ntheta, Nphi])

divd2 = np.zeros([Ntheta,Nphi])
divergence = np.zeros([Ntheta, Nphi])
volume = np.zeros([Ntheta, Nphi])

# Calculates all the field and direction derivatives.
h_r = radius[2] - radius[1]
h_theta = rc * (theta[2] - theta[1])
for j in range(1, Ntheta - 1):
    h_phi = rc * np.sin(theta[j]) * (phi[2] - phi[1])

    derivphifieldr[j, 0] = .5 * (fieldr[i, j, 1] - fieldr[i, j, Nphi - 1]) / h_phi
    derivphidirectionr[j, 0] = .5 * (directionr[i, j, 1] - directionr[i, j, Nphi - 1]) / h_phi
    derivphifieldtheta[j, 0] = .5 * (fieldtheta[i, j, 1] - fieldtheta[i, j, Nphi - 1]) / h_phi
    derivphidirectiontheta[j, 0] = .5 * (directiontheta[i, j, 1] - directiontheta[i, j, Nphi - 1]) / h_phi
    derivphifieldphi[j, 0] = .5 * (fieldphi[i, j, 1] - fieldphi[i, j, Nphi - 1]) / h_phi
    derivphidirectionphi[j, 0] = .5 * (directionphi[i, j, 1] - directionphi[i, j, Nphi - 1]) / h_phi

    derivphifieldr[j, Nphi-1] = .5 * (fieldr[i, j, 0] - fieldr[i, j, Nphi - 2]) / h_phi
    derivphidirectionr[j, Nphi-1] = .5 * (directionr[i, j, 0] - directionr[i, j, Nphi - 2]) / h_phi
    derivphifieldtheta[j, Nphi-1] = .5 * (fieldtheta[i, j, 0] - fieldtheta[i, j, Nphi - 2]) / h_phi
    derivphidirectiontheta[j, Nphi-1] = .5 * (directiontheta[i, j, 0] - directiontheta[i, j, Nphi - 2]) / h_phi
    derivphifieldphi[j, Nphi-1] = .5 * (fieldphi[i, j, 0] - fieldphi[i, j, Nphi - 2]) / h_phi
    derivphidirectionphi[j, Nphi-1] = .5 * (directionphi[i, j, 0] - directionphi[i, j, Nphi - 2]) / h_phi

    for k in range(1, Nphi-1):
        i = 1
        volume[j,k] = h_r*h_phi*h_theta
        divd2[j,k] = h_r**(-2) + h_theta**(-2) + h_phi**(-2)

        derivrfieldr[j, k] = .5 * (fieldr[i + 1, j, k] - fieldr[i - 1, j, k]) / h_r
        derivrdirectionr[j, k] = .5 * (directionr[i + 1, j, k] - directionr[i - 1, j, k]) / h_r
        derivrfieldtheta[j, k] = .5 * (fieldtheta[i + 1, j, k] - fieldtheta[i - 1, j, k]) / h_r
        derivrdirectiontheta[j, k] = .5 * (directiontheta[i + 1, j, k] - directiontheta[i - 1, j, k]) / h_r
        derivrfieldphi[j, k] = .5 * (fieldphi[i + 1, j, k] - fieldphi[i - 1, j, k]) / h_r
        derivrdirectionphi[j, k] = .5 * (directionphi[i + 1, j, k] - directionphi[i - 1, j, k]) / h_r

        derivthetafieldr[j, k] = .5 * (fieldr[i, j + 1, k] - fieldr[i, j - 1, k]) / h_theta
        derivthetadirectionr[j, k] = .5 * (directionr[i, j + 1, k] - directionr[i, j - 1, k]) / h_theta
        derivthetafieldtheta[j, k] = .5 * (fieldtheta[i, j + 1, k] - fieldtheta[i, j - 1, k]) / h_theta
        derivthetadirectiontheta[j, k] = .5 * (directiontheta[i, j + 1, k] - directiontheta[i, j - 1, k]) / h_theta
        derivthetafieldphi[j, k] = .5 * (fieldphi[i, j + 1, k] - fieldphi[i, j - 1, k]) / h_theta
        derivthetadirectionphi[j, k] = .5 * (directionphi[i, j + 1, k] - directionphi[i, j - 1, k]) / h_theta

        derivphifieldr[j, k] = .5 * (fieldr[i, j, k + 1] - fieldr[i, j, k - 1]) / h_phi
        derivphidirectionr[j, k] = .5 * (directionr[i, j, k + 1] - directionr[i, j, k - 1]) / h_phi
        derivphifieldtheta[j, k] = .5 * (fieldtheta[i, j, k + 1] - fieldtheta[i, j, k - 1]) / h_phi
        derivphidirectiontheta[j, k] = .5 * (directiontheta[i, j, k + 1] - directiontheta[i, j, k - 1]) / h_phi
        derivphifieldphi[j, k] = .5 * (fieldphi[i, j, k + 1] - fieldphi[i, j, k - 1]) / h_phi
        derivphidirectionphi[j, k] = .5 * (directionphi[i, j, k + 1] - directionphi[i, j, k - 1]) / h_phi
        
        divergence[j, k] = 0.5/rc**2 * ( (rc + dr)**2*fieldr[i+1,j,1] - (rc - dr)**2*fieldr[i-1,j,1] ) / h_r + \
                       0.5/(rc*np.sin(theta[j])) * ( fieldtheta[i,j+1,k]*np.sin(theta[j+1]) - fieldtheta[i,j-1,k]*np.sin(theta[j-1]) )/ h_theta + \
                       derivphifieldphi[j, k] / rc / np.sin(theta[j])

curlr = np.zeros([Ntheta, Nphi])
curltheta = np.zeros([Ntheta, Nphi])
curlphi = np.zeros([Ntheta, Nphi])
curlmod = np.zeros([Ntheta, Nphi])
curvaturer = np.zeros([Ntheta, Nphi])
curvaturetheta = np.zeros([Ntheta, Nphi])
curvaturephi = np.zeros([Ntheta, Nphi])
curvaturemod = np.zeros([Ntheta, Nphi])

for j in range(1, Ntheta - 1):
    i = 1

#    divergence[j, :] = 2. * fieldr[i, j, :] / rc + derivrfieldr[j, :] + \
#                       fieldtheta[i, j, :] * np.cos(theta[j]) / rc / np.sin(theta[j]) + \
#                       derivthetafieldtheta[j, :] / rc + derivphifieldphi[j, :] / rc / np.sin(theta[j])

    curlr[j, :] = fieldphi[i, j, :] / rc / np.tan(theta[j]) + \
                  derivthetafieldphi[j, :] / rc - derivphifieldtheta[j, :] / rc / np.sin(theta[j])
    curltheta[j, :] = derivphifieldr[j, :] / rc / np.sin(theta[j]) - fieldphi[i, j, :] / rc - derivrfieldphi[j, :]
    curlphi[j, :] = fieldtheta[i, j, :] / rc + derivrfieldtheta[j, :] - derivthetafieldr[j, :] / rc

    curvaturer[j, :] = directionr[i, j, :] * derivrdirectionr[j, :] + \
                       directiontheta[i, j, :] * derivthetadirectionr[j, :] / a + \
                       directionphi[i, j, :] * derivthetadirectionr[j, :] / a / np.sin(theta[j]) - \
                       (directiontheta[i, j, :] ** 2 + directionphi[i, j, :] ** 2) / rc

    curvaturetheta[j, :] = directionr[i, j, :] * derivrdirectiontheta[j, :] + \
                           directiontheta[i, j, :] * derivthetadirectiontheta[j, :] / a + \
                           directionphi[i, j, :] * derivthetadirectiontheta[j, :] / a / np.sin(theta[j]) + \
                           directiontheta[i, j, :] * directionr[i, j, :] / rc - \
                           directionphi[i, j, :] ** 2 / rc / np.tan(theta[j])

    curvaturephi[j, :] = directionr[i, j, :] * derivrdirectionphi[j, :] + \
                         directiontheta[i, j, :] * derivthetadirectionphi[j, :] / a + \
                         directionphi[i, j, :] * derivthetadirectionphi[j, :] / a / np.sin(theta[j]) + \
                         directionphi[i, j, :] * directionr[i, j, :] / rc + \
                         directionphi[i, j, :] * directiontheta[i, j, :] / rc / np.tan(theta[j])

curlmod = np.sqrt(curlr ** 2 + curltheta ** 2 + curlphi ** 2)
curvaturemod = np.sqrt(curvaturer ** 2 + curvaturetheta ** 2 + curvaturephi ** 2)

field_L2 = np.sum(volume*fieldmod**2*divd2)/np.sum(volume)
div_L2 = np.sum(volume*divergence**2)/np.sum(volume)
curl_L2 = np.sum(volume*curlmod**2)/np.sum(volume)
print("B/dr, Div and Curl L2:",field_L2,div_L2,curl_L2,np.sum(volume))

# Plot parameters
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['font.size'] = 10
plt.rcParams['lines.linewidth'] = 1
plt.rcParams["figure.autolayout"] = True
cmap1 = cm.get_cmap('RdBu_r', 255)
cmap2 = cm.get_cmap('inferno', 255)
Phi, Theta = np.meshgrid(360 * (1 - phi / 2 / np.pi), - 180 * theta / np.pi + 90)

# Prints the plot you chose at the beginning of the code
if plot2d:
    plt.xticks([0, 60, 120, 180, 240, 300, 360])
    plt.yticks([-90, -60, -30, 0, 30, 60, 90])
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    if Plot == 0:
        plt.contourf(Phi, Theta, divergence, cmap=cmap1, levels=30)
    elif Plot == 1:
        plt.contourf(Phi, Theta, curlmod, cmap=cmap1, levels=30)
    elif Plot == 2:
        plt.contourf(Phi, Theta, curvaturemod, cmap=cmap1, levels=30)
    elif Plot == 3:
        plt.contourf(Phi, Theta, curlr, cmap=cmap1, levels=30)
    elif Plot == 4:
        plt.contourf(Phi, Theta, curltheta, cmap=cmap2, levels=30)
    elif Plot == 5:
        plt.contourf(Phi, Theta, curlphi, cmap=cmap2, levels=30)
    elif Plot == 6:
        plt.contourf(Phi, Theta, curvaturer, cmap=cmap2, levels=30)
    elif Plot == 7:
        plt.contourf(Phi, Theta, curvaturetheta, cmap=cmap2, levels=30)
    elif Plot == 8:
        plt.contourf(Phi, Theta, curvaturephi, cmap=cmap2, levels=30)
    cbar = plt.colorbar(orientation="horizontal", pad=.1, shrink=0.5)
    cbar.set_label('Whatever you chose')
    plt.gca().invert_xaxis()
    plt.show()

# Saves all plots in a png format
if printall:

    names = [r'$\nabla \cdot B h$ (Gauss)', r'$|\nabla x B|h$ (Gauss)', r'$(\nabla x B)_r h$ (Gauss)', r'$(\nabla x B)_\theta h$ (Gauss)', r'$(\nabla x B)_\phi h$ (Gauss)', r'$| \kappa |$', 'B', 'Br', 'Btheta', 'Bphi']
    files = [planet+'_div.png',planet+'_curl_mod.png',planet+'_curlr.png',planet+'_curlth.png',planet+'_curlphi.png',planet+'_k_mod.png',planet+'_b_mod.png',planet+'_br.png',planet+'_btheta.png',planet+'_bphi.png']
    magnitudes = [divergence*h_theta, curlmod*h_theta, curlr*h_theta, curltheta*h_theta, curlphi*h_theta, curvaturemod, fieldmod[1,:,:], fieldr[1,:,:], fieldtheta[1,:,:], fieldphi[1,:,:]]

    maxf = [ 100.,10., 1., 1., 1.,100.,100., 100., 100., 100.]
    minf = [-100.,0.,-1.,-1.,-1.,  0.,  0.,-100.,-100.,-100.]
    printvalues = np.zeros([Ntheta, Nphi])

    for number in range(0,len(magnitudes)):
        plt.clf()
        plt.xticks([0, 60, 120, 180, 240, 300, 360])
        plt.yticks([-90, -60, -30, 0, 30, 60, 90])
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")

        printvalues = magnitudes[number]
        printvalues[np.where(printvalues > maxf[number])] = maxf[number]
        printvalues[np.where(printvalues < minf[number])] = minf[number]
        if number == 1 or number == 2:
            map = cmap2
        else:
            map = cmap1
        plt.contourf(Phi, Theta, printvalues, cmap=map, levels=30)
        plt.gca().invert_xaxis()

        cbar = plt.colorbar(orientation="horizontal", pad=.1, shrink=0.5)
        cbar.set_label(names[number])
        plt.savefig(files[number])
