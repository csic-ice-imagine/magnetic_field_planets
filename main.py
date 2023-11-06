import numpy as np
import matplotlib.pyplot as plt
import reader
import schmidt
import saveoutput
import saveplots
import lowes_spec
import magnitudes

# Number of points in the latitudinal (North-South) direction, which automatically sets the longitudinal (East-West) direction. 
Ntheta = 50
Nphi = 2*Ntheta
# Number of shells to be calculated
Nr = 1 

a = 1  # Should be 6371.2/72492 (Earth/Jupiter), but renormalize to 1, since r/a is what matters
rc = 1.00  # Radius considered in the map plot (CMB = 0.455/0.85 for Earth/Jupiter), and name of the given files
rc_file = str(rc)
rc_file = rc_file.replace(".","_") # String used for naming the output files

# You can choose either Earth, Jupiter, Jupiter_2021, Saturn, Neptune or Uranus, if you put anything else you will have 0 to everything. 
# If you choose Earth, you also need to choose a year, which can only be: 1900, 1905, 1910, ..., to 2020.
planet, year = "Earth", 2020

# Saves a csv file with potencial, Br, Btheta, Bphi and Bmod. Use with only 1 Nr, for spherical plots
# Saves a vtu file (Paraview for 3D visualization) with Bx, By, Bz. Use more than 1 Nr, as cells are used
filecvs, filevtu = False, False
 
# Saves all 2D/Mollweide plots in png format (in case of the Earth it plots the coastline behind)
planeproj, mollweideproj = True, True

# Calculates many magnitudes
plot_magnitudes = False

# Definition of the spherical grid matrices
#phi    = np.linspace(1 * np.pi / Nphi, 2 * np.pi *(1 + 1/Nphi), num=Nphi)
phi    = np.linspace(0, 2*np.pi, num=Nphi)
theta  = np.linspace(.1 * np.pi / Ntheta, np.pi * (1 - .1 / Ntheta), num=Ntheta)
if magnitudes:
    Nr = 3     # To calulate derivatives and laplacians in a given radius we need at least two other shells of points
    dr = 0.001 # Radial distance between the three shells (in units of the chosen shell set to 1)
    radius = np.linspace(rc * (1. - dr), rc * (1. + dr), num=Nr)
else: radius = np.linspace(rc, rc + 3 * rc, num=Nr)

# Saves the Lowes spectrum for the given radius
lowes = True
# Saves the Lowes spectrum for a number of radii
multiple_lowes_r, lowes_radii = False, np.array([1.45,1.30,1.15,1.00,0.85,0.70,0.55])

if planet=="Earth":
    NPOL=14
    const = 1e5  # To go from nanotesla to gauss (usually the plots are using gauss) if necessary
elif planet=="Jupiter":
    NPOL=11
    const = 1e5
elif planet=="Jupiter_2021":
    NPOL=13
    const = 1e5
elif planet=="Saturn":
    NPOL=7
    const = 1e5
elif planet=="Uranus":
    NPOL=4
    const=1      # Uranus and Neptune coefficients are already in Gauss 
elif planet=="Neptune":
    NPOL=4
    const=1


# Initialize all components of the magnetic field (spherical, cartesian and modulus)
potential = np.zeros([Nr, Ntheta, Nphi])
fieldr = np.zeros([Nr, Ntheta, Nphi])
fieldtheta = np.zeros([Nr, Ntheta, Nphi])
fieldphi = np.zeros([Nr, Ntheta, Nphi])
fieldmod = np.zeros([Nr, Ntheta, Nphi])
fieldx = np.zeros([Nr, Ntheta, Nphi])
fieldy = np.zeros([Nr, Ntheta, Nphi])
fieldz = np.zeros([Nr, Ntheta, Nphi])

# This part defines the K and S matrices with dimension NPOL x NPOL, depending on the 
# planet and the year
g, h = reader.reader(planet, year, NPOL)

# This part defines the K and S matrices with dimension NPOL x NPOL
K, S = schmidt.KandS(NPOL)

# This part defines the Gaussian-normalized and
# the Schmidt quasi-normalized associated Legendre polynomials
P, derivP = schmidt.Schmidtcoefficients(NPOL, Ntheta, theta, K, S)

# Loops for all r, theta and phi and obtaining all the corresponding potential and fields
for j in range(0, Ntheta):
    for k in range(0, Nphi):
        potential[:, j, k], fieldr[:, j, k], fieldtheta[:, j, k], fieldphi[:, j, k] = schmidt.potentialfunction(radius[:], j, phi[k], theta, NPOL, P, derivP, const, g, h)
        fieldmod[:, j, k] = np.sqrt(fieldr[:, j, k] ** 2 + fieldtheta[:, j, k] ** 2 + fieldphi[:, j, k] ** 2)
        fieldx[:, j, k] = np.cos(phi[k]) * np.sin(theta[j]) * fieldr[:, j, k] + np.cos(phi[k]) * np.cos(theta[j]) * \
                            fieldtheta[:, j, k] - np.sin(phi[j]) * fieldphi[:, j, k]
        fieldy[:, j, k] = np.sin(phi[k]) * np.sin(theta[j]) * fieldr[:, j, k] + np.sin(phi[k]) * np.cos(theta[j]) * \
                            fieldtheta[:, j, k] + np.cos(phi[j]) * fieldphi[:, j, k]
        fieldz[:, j, k] = np.cos(theta[j]) * fieldr[:, j, k] - np.sin(theta[j]) * fieldtheta[:, j, k]

# Prints a cvs files with all calculated spherical quantities
if filecvs: saveoutput.savecsv(Nr, Ntheta, Nphi, radius, phi, theta, potential, fieldr, fieldtheta, fieldphi, fieldmod, planet, year)

# Prints a vtu file with the cartesian magnetic field, for Paraview 3D visualization
if filevtu: saveoutput.savevtu(Nr, Ntheta, Nphi, radius, phi, theta, fieldx, fieldy, fieldz, planet, year)

plt.rcParams['figure.figsize'] = (10, 7)
plt.rcParams['font.size'] = 16
plt.rcParams['lines.linewidth'] = 1
plt.rcParams["figure.autolayout"] = True

# Saves plane projection of the magnitudes
if planeproj: saveplots.planeproj(planet, rc, rc_file, phi, theta, Nr, Ntheta, Nphi, potential, fieldr, fieldtheta, fieldphi, fieldmod)

# Saves Mollweide projection of the magnitudes
if mollweideproj: saveplots.mollweideproj(planet, rc, rc_file, phi, theta, Nr, Ntheta, Nphi, potential, fieldr, fieldtheta, fieldphi, fieldmod)

if plot_magnitudes: magnitudes.printMagnitudes(planet, Ntheta, Nphi, radius, rc, rc_file, a, dr, phi, theta, fieldr, fieldtheta, fieldphi, fieldmod)

# Obtain the Lowes spectrum for a the given plotted radius and plot it
if lowes:
    Rn = lowes_spec.lowes_spec(NPOL, rc, g, h)
    lowes_spec.plot_lowes(planet, rc, rc_file, Rn)

# Obtian many Lowes spectrum ad different radii and plot them
if multiple_lowes_r:
    Rn = np.zeros([len(lowes_radii), NPOL])
    for r in range(0,len(lowes_radii)):
        Rn[r,:] = lowes_spec.lowes_spec(NPOL, lowes_radii[r], g, h)
    lowes_spec.plot_multiple_lowes(planet, lowes_radii, Rn)
