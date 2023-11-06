import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import cm

# Plot parameters
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['font.size'] = 10
plt.rcParams['lines.linewidth'] = 1
plt.rcParams["figure.autolayout"] = True
cmap1 = cm.get_cmap('RdBu_r', 255)
cmap2 = cm.get_cmap('inferno', 255)


def printMagnitudes(planet, Ntheta, Nphi, radius, rc, rc_file, a, dr, phi, theta, fieldr, fieldtheta, fieldphi, fieldmod):

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
            divd2[j,k] = h_r**(2) + h_theta**(2) + h_phi**(2)

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

        divergence[j, :] = 2. * fieldr[i, j, :] / rc + derivrfieldr[j, :] + \
                        fieldtheta[i, j, :] * np.cos(theta[j]) / rc / np.sin(theta[j]) + \
                        derivthetafieldtheta[j, :] / rc + derivphifieldphi[j, :] / rc / np.sin(theta[j])

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

    field_L2 = np.sum(volume*fieldmod**2)
    div_L2 = np.sum(volume*divergence**2*divd2)
    curl_L2 = np.sum(volume*curlmod**2*divd2)
    print("B**2, (Div*h)**2 and (Curl*h)**2, integrated:",field_L2,div_L2,curl_L2,np.sum(volume),"dr,Nth,Np:",dr,Ntheta,Nphi)


    names = [r'$\nabla \cdot B h$ (Gauss)', r'$|\nabla x B|h$ (Gauss)', r'$(\nabla x B)_r h$ (Gauss)', r'$(\nabla x B)_\theta h$ (Gauss)', r'$(\nabla x B)_\phi h$ (Gauss)', r'$| \kappa |$', 'B', 'Br', 'Btheta', 'Bphi']
    files = [planet+'_div.png',planet+'_curl_mod.png',planet+'_curlr.png',planet+'_curlth.png',planet+'_curlphi.png',planet+'_k_mod.png',planet+'_b_mod.png',planet+'_br.png',planet+'_btheta.png',planet+'_bphi.png']
    magnitudes = [divergence*np.sqrt(divd2), curlmod*np.sqrt(divd2), curlr*np.sqrt(divd2), curltheta*np.sqrt(divd2), curlphi*np.sqrt(divd2), curvaturemod, fieldmod[1,:,:], fieldr[1,:,:], fieldtheta[1,:,:], fieldphi[1,:,:]]

    maxf = [ 100.,10., 1., 1., 1.,100.,100., 100., 100., 100.]
    minf = [-100.,0.,-1.,-1.,-1.,  0.,  0.,-100.,-100.,-100.]
    printvalues = np.zeros([Ntheta, Nphi])
    Phi, Theta = np.meshgrid(360 * (1 - phi / 2 / np.pi), - 180 * theta / np.pi + 90)

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