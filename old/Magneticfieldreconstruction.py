import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import csv

# Spherical grid numbers
Nphi = 100
Ntheta = 50
Nr = 1
a = 1  # Should be 6371.2/72492 (Earth/Jupiter), but renormalize to 1, since r/a is what matters
rc = .85  # Radius considered in the map plot (CMB = 0.455/0.85 for earth/jupiter).
# Remember to change the name of the plots to include the distance!!!
const = 1e5  # To fo from nanotesla to gauss (usually the plots are using gauss)

planet = "Jupiter"  # Can be Earth, Jupiter
year = 2020  # Can be 1900, 1905, 1910, ..., to 2020. Only in the case of the Earth

filecvs = False  # Saves a cvs file with potencial, Br, Btheta, Bphi and Bmod. Use with only 1 Nr, for spherical plots
filevtu = False  # Saves a vtu file (Paraview) with Bx, By, Bz. Use more than 1 Nr, as cells are used

printall = False  # Saves all 2D plots

plot2d = True  # Shows a given plot (for testing)
Plot = 1
# if 0 you will plot the potential
# if 1 the radial field
# if 2 the theta field
# if 3 the phi field
# if 4 the field modulus

# Definition of the spherical grid
phi = np.linspace(0, 2 * np.pi, num=Nphi)
theta = np.linspace(.1 * np.pi / Ntheta, np.pi * (1 - .1 / Ntheta), num=Ntheta)
radius = np.linspace(rc, rc + 9 * rc, num=Nr)
potential = np.zeros([Nr, Ntheta, Nphi])
# Initialize all components of the magnetic field (shperical, cartesian and modulus)
fieldr = np.zeros([Nr, Ntheta, Nphi])
fieldtheta = np.zeros([Nr, Ntheta, Nphi])
fieldphi = np.zeros([Nr, Ntheta, Nphi])
fieldmod = np.zeros([Nr, Ntheta, Nphi])
fieldx = np.zeros([Nr, Ntheta, Nphi])
fieldy = np.zeros([Nr, Ntheta, Nphi])
fieldz = np.zeros([Nr, Ntheta, Nphi])

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

# This part defines all K's and S's used in the recursive formulas for
# the Schmidt quasi-normalized associated Legendre polynomials
K = np.zeros([NPOL, NPOL])
S = np.zeros([NPOL, NPOL])
S[0, 0] = 1
delta = 0
for n in range(1, NPOL):
    for m in range(0, n + 1):
        K[n, m] = ((n - 1) ** 2 - m ** 2) / (2 * n - 1) / (2 * n - 3)
        if n == 1:
            K[n, m] = 0
        if m == 1:
            delta = 1
        else:
            delta = 0
        if m == 0:
            S[n, 0] = (2 * n - 1) * S[n - 1, 0] / n
        else:
            S[n, m] = S[n, m - 1] * np.sqrt((n - m + 1) * (delta + 1) / (n + m))

# This part defines the Gaussian-normalized and
# the Schmidt quasi-normalized associated Legendre polynomials
Pgauss = np.zeros([NPOL, NPOL, Ntheta])
derivPgauss = np.zeros([NPOL, NPOL, Ntheta])
P = np.zeros([NPOL, NPOL, Ntheta])
derivP = np.zeros([NPOL, NPOL, Ntheta])
for ntheta in range(0, Ntheta):
    Pgauss[0, 0, ntheta] = 1
    derivPgauss[0, 0, ntheta] = 0

for ntheta in range(0, Ntheta):
    for n in range(1, NPOL):
        for m in range(0, n + 1):
            for i in range(1, n + 1):
                for j in range(0, i + 1):
                    if i == j:
                        Pgauss[i, i, ntheta] = np.sin(theta[ntheta]) * Pgauss[i - 1, i - 1, ntheta]
                        derivPgauss[i, i, ntheta] = np.sin(theta[ntheta]) * derivPgauss[i - 1, i - 1, ntheta] + np.cos(
                            theta[ntheta]) * Pgauss[i - 1, i - 1, ntheta]
                    elif i == 1:
                        Pgauss[1, 0, ntheta] = np.cos(theta[ntheta])
                        derivPgauss[1, 0, ntheta] = - np.sin(theta[ntheta])
                    else:
                        Pgauss[i, j, ntheta] = np.cos(theta[ntheta]) * Pgauss[i - 1, j, ntheta] - K[i, j] * Pgauss[
                            i - 2, j, ntheta]
                        derivPgauss[i, j, ntheta] = np.cos(theta[ntheta]) * derivPgauss[i - 1, j, ntheta] - np.sin(
                            theta[ntheta]) * Pgauss[i - 1, j, ntheta] - K[i, j] * derivPgauss[i - 2, j, ntheta]
            P[n, m, ntheta] = S[n, m] * Pgauss[n, m, ntheta]
            derivP[n, m, ntheta] = S[n, m] * derivPgauss[n, m, ntheta]


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
            fieldmod[i, j, k] = np.sqrt(fieldr[i, j, k] ** 2 + fieldtheta[i, j, k] ** 2 + fieldphi[i, j, k] ** 2)
            fieldx[i, j, k] = np.cos(phi[k]) * np.sin(theta[j]) * fieldr[i, j, k] + np.cos(phi[k]) * np.cos(theta[j]) * \
                              fieldtheta[i, j, k] - np.sin(phi[j]) * fieldphi[i, j, k]
            fieldy[i, j, k] = np.sin(phi[k]) * np.sin(theta[j]) * fieldr[i, j, k] + np.sin(phi[k]) * np.cos(theta[j]) * \
                              fieldtheta[i, j, k] + np.cos(phi[j]) * fieldphi[i, j, k]
            fieldz[i, j, k] = np.cos(theta[j]) * fieldr[i, j, k] - np.sin(theta[j]) * fieldtheta[i, j, k]

# Prints a cvs files with all calculated spherical quantities
if filecvs:
    if planet == "Earth":
        names = ['Earthfieldpot_' + str(year) + '.csv', 'Earthfieldr_' + str(year) + '.csv',
                 'Earthfieldtheta_' + str(year) + '.csv', 'Earthfieldphi_' + str(year) + '.csv',
                 'Earthfieldmod_' + str(year) + '.csv']
    elif planet == "Jupiter":
        names = ['Jupiterfieldpot.csv', 'Jupiterfieldr.csv', 'Jupiterfieldtheta.csv', 'Jupiterfieldphi.csv',
                 'Jupiterfieldmod.csv']
    else:
        names = ['error.csv', 'error.csv', 'error.csv', 'error.csv', 'error.csv']

    heads = [["X", "Y", "Z", "Pot"], ["X", "Y", "Z", "Br"], ["X", "Y", "Z", "Btheta"], ["X", "Y", "Z", "Bphi"],
             ["X", "Y", "Z", "Bmod"]]
    magntidues = [potential, fieldr, fieldtheta, fieldphi, fieldmod]
    printvalues = np.zeros([Nr, Ntheta, Nphi])
    for number in range(0, 5):
        outputfile = names[number]
        head = heads[number]
        printvalues = magntidues[number]
        with open(outputfile, 'w', newline='') as file:
            writer = csv.writer(file, delimiter=',')
            writer.writerow(head)
            for i in range(0, Nr):
                for j in range(0, Ntheta):
                    for k in range(0, Nphi):
                        values = [radius[i] * np.cos(phi[k]) * np.sin(theta[j]),
                                  radius[i] * np.sin(phi[k]) * np.sin(theta[j]), radius[i] * np.cos(theta[j]),
                                  printvalues[i, j, k]]
                        writer.writerow(values)

# Prints a vtu file with the cartesian magnetic field, for Paraview 3D visualization
if filevtu:
    cellx = np.zeros([Nr, Ntheta, Nphi])
    celly = np.zeros([Nr, Ntheta, Nphi])
    cellz = np.zeros([Nr, Ntheta, Nphi])

    for i in range(0, Nr - 1):
        for j in range(0, Ntheta - 1):
            for k in range(0, Nphi - 1):
                cellx[i, j, k] = 0.125 * (
                            fieldx[i, j, k] + fieldx[i, j + 1, k] + fieldx[i, j + 1, k + 1] + fieldx[i, j, k + 1] +
                            fieldx[i + 1, j, k] + fieldx[i + 1, j + 1, k] + fieldx[i + 1, j + 1, k + 1] + fieldx[
                                i + 1, j, k + 1])
                celly[i, j, k] = 0.125 * (
                            fieldy[i, j, k] + fieldy[i, j + 1, k] + fieldy[i, j + 1, k + 1] + fieldy[i, j, k + 1] +
                            fieldy[i + 1, j, k] + fieldy[i + 1, j + 1, k] + fieldy[i + 1, j + 1, k + 1] + fieldy[
                                i + 1, j, k + 1])
                cellz[i, j, k] = 0.125 * (
                            fieldz[i, j, k] + fieldz[i, j + 1, k] + fieldz[i, j + 1, k + 1] + fieldz[i, j, k + 1] +
                            fieldz[i + 1, j, k] + fieldz[i + 1, j + 1, k] + fieldz[i + 1, j + 1, k + 1] + fieldz[
                                i + 1, j, k + 1])

    N_Points = Nr * Ntheta * Nphi
    N_Cells = (Nr - 1) * (Ntheta - 1) * (Nphi - 1)
    id_nodes = np.zeros([Nr + 1, Ntheta, Nphi], dtype=int)

    if planet == "Earth":
        name = 'Earthfieldvec.vtu'
    elif planet == "Jupiter":
        name = 'Jupiterfieldvec.vtu'
    else:
        name = 'error.vtu'
    with open(name, 'w') as file:
        file.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        file.write('  <UnstructuredGrid>\n')
        file.write('    <Piece NumberOfPoints="{}" NumberOfCells="{}">\n'.format(N_Points, N_Cells))
        file.write('      <Points>\n')
        file.write('        <DataArray type="Float32" format="ascii" NumberOfComponents="3">\n')
        identif = 0
        for i in range(0, Nr):
            for j in range(0, Ntheta):
                for k in range(0, Nphi):
                    id_nodes[i, j, k] = identif
                    file.write('        {}   {}   {}\n'.format(radius[i] * np.cos(phi[k]) * np.sin(theta[j]),
                                                               radius[i] * np.sin(phi[k]) * np.sin(theta[j]),
                                                               radius[i] * np.cos(theta[j])))
                    identif += 1
        file.write('        </DataArray>\n')
        file.write('      </Points>\n')
        file.write('      <Cells>\n')
        file.write('        <DataArray type="Int32" Name="connectivity" format="ascii">\n')
        for i in range(0, Nr - 1):
            for j in range(0, Ntheta - 1):
                for k in range(0, Nphi - 1):
                    file.write('        {}    {}    {}    {}    {}    {}    {}    {}\n'.format(id_nodes[i, j, k],
                                                                                               id_nodes[i, j + 1, k],
                                                                                               id_nodes[
                                                                                                   i, j + 1, k + 1],
                                                                                               id_nodes[i, j, k + 1],
                                                                                               id_nodes[i + 1, j, k],
                                                                                               id_nodes[
                                                                                                   i + 1, j + 1, k],
                                                                                               id_nodes[
                                                                                                   i + 1, j + 1, k + 1],
                                                                                               id_nodes[
                                                                                                   i + 1, j, k + 1]))
        file.write('        </DataArray>\n')
        file.write('        <DataArray type="Int32" Name="offsets" format="ascii">\n')
        for i in range(1, N_Cells + 1):
            file.write('      {}\n'.format(i * 8))
        file.write('        </DataArray>\n')
        file.write('        <DataArray type="UInt8" Name="types" format="ascii">\n')
        for i in range(1, N_Cells + 1):
            file.write('      {}\n'.format(12))
        file.write('        </DataArray>\n')
        file.write('      </Cells>\n')
        file.write('      <PointData>\n')
        file.write('      </PointData>\n')
        file.write('      <CellData>\n')
        file.write('        <DataArray Name="B" NumberOfComponents="3" type="Float32" format="ascii">')
        for i in range(0, Nr - 1):
            for j in range(0, Ntheta - 1):
                for k in range(0, Nphi - 1):
                    file.write('        {}    {}    {}\n'.format(cellx[i, j, k], celly[i, j, k], cellz[i, j, k]))
        file.write('        </DataArray>\n')
        file.write('      </CellData>\n')
        file.write('    </Piece>\n')
        file.write('  </UnstructuredGrid>\n')
        file.write('</VTKFile>\n')

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
        plt.contourf(Phi, Theta, potential[0, :, :], cmap=cmap1, levels=30)
    elif Plot == 1:
        plt.contourf(Phi, Theta, fieldr[0, :, :], cmap=cmap1, levels=30)
    elif Plot == 2:
        plt.contourf(Phi, Theta, fieldtheta[0, :, :], cmap=cmap1, levels=30)
    elif Plot == 3:
        plt.contourf(Phi, Theta, fieldphi[0, :, :], cmap=cmap1, levels=30)
    elif Plot == 4:
        plt.contourf(Phi, Theta, fieldmod[0, :, :], cmap=cmap2, levels=30)
    cbar = plt.colorbar(orientation="horizontal", pad=.1, shrink=0.5)
    cbar.set_label('Whatever you chose')
    plt.gca().invert_xaxis()
    plt.show()

# Saves all plots in a png format
if printall:
    if planet == "Earth":
        names = ['Potential (Gauss · 1 $R_T$) at $r =$ $R_T$', '$B_r$ (Gauss) at $r =$ $R_T$',
                 '$B_θ$ (Gauss) at $r =$ $R_T$', '$B_φ$ (Gauss) at $r =$ $R_T$', '$|B|$ (Gauss) at $r =$ $R_T$']
        files = ['Earthpotential.png', 'Earthfieldr.png', 'Earthfieldtheta.png', 'Earthfieldphi.png',
                 'Earthfieldmod.png']
    elif planet == "Jupiter":
        names = ['Potential (Gauss · 1 $R_J$) at $r = 0.85$ $R_J$', '$B_r$ (Gauss) at $r = 0.85$ $R_J$',
                 '$B_θ$ (Gauss) at $r = 0.85$ $R_J$', '$B_φ$ (Gauss) at $r = 0.85$ $R_J$',
                 '$|B|$ (Gauss) at $r = 0.85$ $R_J$']
        files = ['Jupiterpotential.png', 'Jupiterfieldr.png', 'Jupiterfieldtheta.png', 'Jupiterfieldphi.png',
                 'Jupiterfieldmod.png']
    else:
        names = ['Potential?', '$B_r$ ?', '$B_θ$ ?', '$B_φ$ ?', '$|B|$ ?']
        files = ['errorpot.png', 'errorr.png', 'errortheta.png', 'errorphi.png', 'errormod.png']
    magntidues = [potential, fieldr, fieldtheta, fieldphi, fieldmod]
    printvalues = np.zeros([Nr, Ntheta, Nphi])
    for number in range(0, 5):
        plt.clf()
        plt.xticks([0, 60, 120, 180, 240, 300, 360])
        plt.yticks([-90, -60, -30, 0, 30, 60, 90])
        plt.xlabel("Longitude")
        plt.ylabel("Latitude")
        printvalues = magntidues[number]
        if number == 4:
            realmap = cmap2
        else:
            realmap = cmap1
        plt.contourf(Phi, Theta, printvalues[0, :, :], cmap=realmap, levels=30)
        plt.gca().invert_xaxis()
        cbar = plt.colorbar(orientation="horizontal", pad=.1, shrink=0.5)
        cbar.set_label(names[number])
        plt.savefig(files[number])
