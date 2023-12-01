import numpy as np
import matplotlib.pyplot as plt
import reader, schmidt, saveoutput, saveplots, lowes_spec, magnitudes

# Plot resolution
Ntheta = 50      # Latitudinal points (North-South direction)
Nphi = 2*Ntheta  # Longitudial points (East-West direction)
Nr = 1           # Radial points (change only to generate 3D output)

a = 1  # Should be 6371.2/72492 (Earth/Jupiter), but renormalize to 1, since r/a is what matters. It is used only for the calculation of some magnitudes
rc = 1.00 # Radius considered in the map plot (CMB = 0.455/0.85 for Earth/Jupiter), and name of the given files

# Planet (or satellite) to choose. Raw data is located in folder data/
planet, year = "Jupiter_2021", 2020
# You can choose either Earth, Jupiter, Jupiter_2021, Saturn, Neptune, Uranus, Mercury
# and Ganymede. Anything else will make the code stop.  If you choose Earth, you also need
# to choose a year, which can only be: 1900, 1905, 1910, ..., to 2020.
 
# Switches to save projections in plane and Mollweide projections. Coastlines are included in Earth plots.
planeproj, mollweideproj = True, True
# If you have successfully installed the ccrs library you can put the Earth coastline in the Earth plane projections also, using the boolean ccrs_library
ccrs_library = True

# ATTENTION: To plot using the Mollweide projection you need the ccrs library. The combination mollweideproj=True, ccrs_library=False will crash if you have 
# not installed this library

# Switch to plot the curl, divergence and curvature of the magnetic field
plot_magnitudes = True

# Switch to save the Lowes spectrum for the given radius
lowes = True
# Switch to save the Lowes spectrum for a number of radii
multiple_lowes_r, lowes_radii = True, np.array([1.45,1.30,1.15,1.00,0.85,0.70,0.55])

# Saves a csv file with potencial, Br, Btheta, Bphi and Bmod. Use with only 1 Nr, for spherical plots
# Saves a vtu file (Paraview for 3D visualization) with Bx, By, Bz. Use more than 1 Nr, as cells are used
filecvs, filevtu = False, False


rc_file = '%.2f'%rc
rc_file = rc_file.replace(".","_") # String used for naming the output files

# Definition of the spherical grid matrices
#phi    = np.linspace(1 * np.pi / Nphi, 2 * np.pi *(1 + 1/Nphi), num=Nphi)
phi    = np.linspace(0, 2*np.pi, num=Nphi)
theta  = np.linspace(.1 * np.pi / Ntheta, np.pi * (1 - .1 / Ntheta), num=Ntheta)
if plot_magnitudes:
    Nr = 3     # To calulate derivatives and laplacians in a given radius we need at least two other shells of points
    dr = 0.001 # Radial distance between the three shells (in units of the chosen shell set to 1)
    radius = np.linspace(rc * (1. - dr), rc * (1. + dr), num=Nr)
else: radius = np.linspace(rc, rc + 3 * rc, num=Nr)


# Depending on the planet you choose, the data will have different multipole definition and will be normalized in nT or G
if planet=="Earth":
    NPOL,NPOL_EXT=14,0
    const=1e5  # To go from nanotesla to gauss (usually the plots are using gauss) if necessary
elif planet=="Jupiter":
    NPOL,NPOL_EXT=11,0
    const=1e5
elif planet=="Jupiter_2021":
    NPOL,NPOL_EXT=31,2
    const=1e5
elif planet=="Saturn":
    NPOL,NPOL_EXT=7,0
    const=1e5
elif planet=="Uranus":
    NPOL,NPOL_EXT=4,0
    const=1      # Uranus and Neptune coefficients are already in Gauss 
elif planet=="Neptune":
    NPOL,NPOL_EXT=4,0
    const=1
elif planet=="Mercury":
    NPOL,NPOL_EXT=4,2
    const=1e5
elif planet=="Ganymede":
    NPOL,NPOL_EXT=3,0
    const=1e5
else:
    NPOL,NPOL_EXT=0,0

# This part defines the K and S matrices with dimension NPOL x NPOL, depending on the 
# planet and the year
g, h, G, H = reader.reader(planet, year, NPOL, NPOL_EXT)

# This part defines the K and S matrices with dimension NPOL x NPOL
K, S = schmidt.KandS(NPOL)

# This part defines the Gaussian-normalized and
# the Schmidt quasi-normalized associated Legendre polynomials
P, derivP = schmidt.Schmidtcoefficients(NPOL, Ntheta, theta, K, S)

# Initialize all components of the magnetic field (spherical, cartesian and modulus)
potential, potential_EXT = np.zeros([Nr, Ntheta, Nphi]), np.zeros([Nr, Ntheta, Nphi])
fieldr, fieldr_EXT = np.zeros([Nr, Ntheta, Nphi]), np.zeros([Nr, Ntheta, Nphi])
fieldtheta, fieldtheta_EXT = np.zeros([Nr, Ntheta, Nphi]), np.zeros([Nr, Ntheta, Nphi])
fieldphi, fieldphi_EXT = np.zeros([Nr, Ntheta, Nphi]), np.zeros([Nr, Ntheta, Nphi])
fieldmod = np.zeros([Nr, Ntheta, Nphi])
fieldx = np.zeros([Nr, Ntheta, Nphi])
fieldy = np.zeros([Nr, Ntheta, Nphi])
fieldz = np.zeros([Nr, Ntheta, Nphi])

# Loops for all r, theta and phi and obtaining all the corresponding potential and fields
for j in range(0, Ntheta):
    for k in range(0, Nphi):
        potential[:, j, k], fieldr[:, j, k], fieldtheta[:, j, k], fieldphi[:, j, k] = schmidt.potentialfunction(radius[:], j, phi[k], theta, NPOL, P, derivP, const, g, h)
        if NPOL_EXT != 0:
            potential_EXT[:, j, k], fieldr_EXT[:, j, k], fieldtheta_EXT[:, j, k], fieldphi_EXT[:, j, k] = schmidt.potentialfunctionexternal(radius[:], j, phi[k], theta, NPOL_EXT, P, derivP, const, G, H)
            potential[:, j, k] += potential_EXT[:, j, k]
            fieldr[:, j, k] += fieldr_EXT[:, j, k]
            fieldtheta[:, j, k] += fieldtheta_EXT[:, j, k]
            fieldphi[:, j, k] += fieldphi_EXT[:, j, k]
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

# Plot parameters
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 16
plt.rcParams['lines.linewidth'] = 1
plt.rcParams["figure.autolayout"] = True

print("---------------------------------------------")
print("Saving plots for " + planet + " r = " + str(rc))

# Saves plane projection of the magnitudes
if planeproj: saveplots.planeproj(planet, rc, rc_file, phi, theta, potential, fieldr, fieldtheta, fieldphi, fieldmod, ccrs_library)

# Saves Mollweide projection of the magnitudes
if mollweideproj: saveplots.mollweideproj(planet, rc, rc_file, phi, theta, potential, fieldr, fieldtheta, fieldphi, fieldmod)

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

print("---------------------------------------------")
