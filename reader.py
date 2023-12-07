# reader.py reads the files in data/ depending on which planet is chosen. 
# Still to be implemented the G and H constants fro Jupiter 2021 for external sources 
# Tables from Jupiter adn Earth are directly obtained from free sources, the others
# have been constructed manually from data on "old" Voyager papers or a reduced table 
# from Saturn (it has few non-zero multipoles).

import numpy as np

def reader(planet, year, NPOL, NPOL_EXT):
    g = np.zeros([NPOL, NPOL])
    h = np.zeros([NPOL, NPOL])
    G = np.zeros([NPOL_EXT, NPOL_EXT])
    H = np.zeros([NPOL_EXT, NPOL_EXT])
    if planet == "Earth":
        # Initialitze the set of g and h constants to be filled by table values.
        # Only half of them is filled (as n<m) and maybe there is a better way to implement this.
        # The variable index is the column corresponding to that year
        index = int((year - 1900) / 5 + 3)
        file = open("data/igrf13coeffs.txt", "r")
        # Remove the introductory lines, careful not to skip the first multipole.
        lines = file.readlines()[4:]
        for n in range(0, len(lines)):
            file_list = [i for i in lines[n].split()]
            if file_list[0] == 'g':
                g[int(file_list[1]), int(file_list[2])] = float(file_list[index])
            else:
                h[int(file_list[1]), int(file_list[2])] = float(file_list[index])
    elif planet == "Jupiter":
        file = open("data/grl57087-sup-0005-2018GL077312-ds01.txt", "r")
        lines = file.readlines()[1:]
        for n in range(0, len(lines)):
            file_list = [i for i in lines[n].split()]
            if file_list[3] == 'g':
                g[int(file_list[4]), int(file_list[5])] = float(file_list[1])
            else:
                h[int(file_list[4]), int(file_list[5])] = float(file_list[1])
    elif planet == "Jupiter_2021":
        file = open("data/2021JE007055-sup-0002-Table+SI-S01.txt", "r")
        lines = file.readlines()[1:]
        for n in range(0, len(lines)):
            file_list = [i for i in lines[n].split()]
            error = float(file_list[2].replace('(','').replace(')',''))
            if n < 195 or n > 959:
                if file_list[3] == 'g':
                    g[int(file_list[4]), int(file_list[5])] = float(file_list[1])
                elif file_list[3] == 'h':
                    h[int(file_list[4]), int(file_list[5])] = float(file_list[1])
                elif file_list[3] == 'G':
                    G[int(file_list[4]), int(file_list[5])] = float(file_list[1])
                elif file_list[3] == 'H':
                    H[int(file_list[4]), int(file_list[5])] = float(file_list[1])
    elif planet == "Saturn":
        file = open("data/saturn_models_dougherty18.txt", "r")
        lines = file.readlines()[2:]
        for n in range(0, len(lines)):
            file_list = [i for i in lines[n].split()]
            if file_list[0] == 'g':
                g[int(file_list[1]), int(file_list[2])] = float(file_list[4])
            else:
                h[int(file_list[1]), int(file_list[2])] = float(file_list[4])
    elif planet == "Uranus":
        file = open("data/uranus_model_connerney87.txt", "r")
        lines = file.readlines()[2:]
        for n in range(0, len(lines)):
            file_list = [i for i in lines[n].split()]
            if file_list[0] == 'g':
                g[int(file_list[1]), int(file_list[2])] = float(file_list[3])
            else:
                h[int(file_list[1]), int(file_list[2])] = float(file_list[3])
    elif planet == "Neptune":
        file = open("data/neptune_models_selesnick92.txt", "r")
        lines = file.readlines()[2:]
        for n in range(0, len(lines)):
            file_list = [i for i in lines[n].split()]
            if file_list[0] == 'g':
                g[int(file_list[1]), int(file_list[2])] = float(file_list[7])
            else:
                h[int(file_list[1]), int(file_list[2])] = float(file_list[7])
    elif planet == "Mercury":
        file = open("data/mercury_model_toepfet21.txt", "r")
        lines = file.readlines()[2:]
        for n in range(0, len(lines)):
            file_list = [i for i in lines[n].split()]
            if file_list[0] == 'g':
                g[int(file_list[1]), int(file_list[2])] = float(file_list[3])
            elif file_list[0] == 'h':
                h[int(file_list[1]), int(file_list[2])] = float(file_list[3])
            elif file_list[0] == 'G':
                G[int(file_list[1]), int(file_list[2])] = float(file_list[3])
            elif file_list[0] == 'H':
                H[int(file_list[1]), int(file_list[2])] = float(file_list[3])
            else:
                continue
    elif planet == "Ganymede":
        file = open("data/ganymede_models_weber22.txt", "r")
        lines = file.readlines()[2:]
        for n in range(0, len(lines)):
            file_list = [i for i in lines[n].split()]
            if file_list[0] == 'g':
                g[int(file_list[1]), int(file_list[2])] = float(file_list[4])
            else:
                h[int(file_list[1]), int(file_list[2])] = float(file_list[4])
    
    return g,h,G,H
