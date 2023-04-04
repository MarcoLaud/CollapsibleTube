# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 13:45:32 2022

@author: marco
"""

from numpy import asarray, cos, genfromtxt, loadtxt, mean, pi, sin, std, zeros
from matplotlib.pyplot import figure, legend, plot, title, xlabel, ylabel
from pathlib import Path
from os.path import basename

E = 1e6
nu = 0.49

folders = ["elle", "dtilde", "gamma"]
letters = ["l", "d", "g"]
analyses = [asarray([11, 12, 13, 14, 15, 16, 17, 18]),
            asarray([30, 35, 40, 45, 50, 55, 60]),
            asarray([6, 7, 8, 9])]

# folders = ["dtilde", "gamma"]
# letters = ["d", "g"]
# analyses = [asarray([30, 35, 45, 50, 55, 60]),
#             asarray([6, 7, 8, 9])]

# folder = "elle"
# letter = "l"
# analysis = asarray([11, 12, 13, 14, 15, 16, 17, 18])

# folder = "dtilde"
# letter = "d"
# analysis = asarray([30, 35, 45, 50, 55, 60])

# folder = "gamma"
# letter = "g"
# analysis = asarray([6, 7, 8, 9])


figure()
title('Tube law - Non-dimensional units')
xlabel('Normalised area [-]')
ylabel('Normalised pressure [-]')
# ylim(-0.25, 0)

for ff, folder in enumerate(folders):
    letter = letters[ff]
    analysis = analyses[ff]
    for ps, case in enumerate(analysis):
    
        pathlist = Path('./{}/{}={}/'.format(folder, letter, str(case))).glob('rtz*.csv')
        path_list = [str(element) for element in pathlist]
        sort_path_list = sorted(path_list,
                                key=lambda i: int(
                                        basename(i).split('_')[1].split('.')[0]))
        
        
        pp = genfromtxt("./{}/pressure.csv".format(folder), delimiter=",", skip_header=1)
        pp = pp[:, 1] * (-1)
        areas = zeros(len(sort_path_list))
        
        for jj, elem in enumerate(sort_path_list):
            flag = False
            area = 0
            data = genfromtxt(elem, delimiter=',', skip_header=1)
            rr = data[:, 1]
            th = data[:, 2]
        
            th, rr = zip(*sorted(zip(th, rr)))
            rr = asarray(rr)
            th = asarray(th)
        
            if th[0] > th[1]:
                rr = rr[::-1]
                th = th[::-1]
        
            if max([th[ii+1]-th[ii] for ii in range(len(th)-1)]) > 1:
                flag = True
                RR = [rr[ii] for ii in range(len(rr)) if rr[ii]*cos(th[ii]) > 0 and
                      rr[ii]*sin(th[ii]) > 0]
                TH = [th[ii] for ii in range(len(rr)) if rr[ii]*cos(th[ii]) > 0 and
                      rr[ii]*sin(th[ii]) > 0]
                rr, th = zip(*sorted(zip(RR, TH)))
        
            for ii in range(len(rr)-1):
                if flag:
                    area += 2.0 * mean([rr[ii], rr[ii+1]])**2 * abs(th[ii] - th[ii+1])
                else:
                    area += 0.5 * mean([rr[ii], rr[ii+1]])**2 * abs(th[ii] - th[ii+1])
        
            print('Deformed area - {} : {} m^2'.format(case, area))
            areas[jj] = area

        if letter == "l":
            aa = (areas)/(pi*0.003**2)*(case/10.)
            bb = (pp*(1-0.49**2)/E)*((case/10.)**(-2))*6/0.12
            plot(aa, bb,
                 label="{}={}".format(letter, case/10),
                 linestyle="--")
    
        elif letter == "d":
            aa = (areas)/(pi*0.003**2)*(1.1)
            bb = (pp*(1-0.49**2)/E)*(1.1**(-2))*(case/5.)/0.12
            plot(aa, bb,
                 label="{}={}".format(letter, case/10),
                 linestyle="--")
    
        elif letter == "g":
            aa = (areas)/(pi*0.003**2)*(1.1)
            bb = (pp*(1-0.49**2)/E)*(1.1**(-2))*6/(case/50.)
            plot(aa, bb,
                 label="{}={}".format(letter, case/10),
                 linestyle="--")
legend()

avg_areas = []
avg_pp = []

for ii, folder in enumerate(folders):
    data_pr = loadtxt("./{}/critical_pr_{}.txt".format(folder, folder))
    if folder == "elle":
        analysis = analyses[ii] / 10
        area_ND = data_pr[:, -1] / (pi*0.003**2) * analysis
        pr_ND = (data_pr[:, 1]*(1-0.49**2)/E)*(analysis)**(-2)*6/0.12
    if folder == "dtilde":
        analysis = analyses[ii] / 5
        area_ND = (data_pr[:, -1])/(pi*0.003**2)*(1.1)
        pr_ND = (data_pr[:, 1]*(1-0.49**2)/E)*(1.1**(-2))*(analysis)/0.12
    if folder == "gamma":
        analysis = analyses[ii] / 50
        area_ND = (data_pr[1:, -1])/(pi*0.003**2)*(1.1)
        pr_ND = (data_pr[1:, 1]*(1-0.49**2)/E)*(1.1**(-2))*6/(analysis)

    avg_areas.append(mean(area_ND))
    avg_pp.append(mean(pr_ND))

area_buck = mean(avg_areas)
area_buck_std = std(avg_areas)
pp_buck = mean(avg_pp) * -1
pp_buck_std = std(avg_pp)


print("Normalised buckling critical area: {} +- {}".format(round(area_buck, 3),
                                                           round(area_buck_std, 3)))
print("Normalised buckling critical pressure: {} +- {}".format(round(pp_buck, 3),
                                                           round(pp_buck_std, 3)))

plot(area_buck, pp_buck, linestyle="", marker="o", color="black")

avg_areas_con = []
avg_pp_con = []

for ii, folder in enumerate(folders):
    data_pr = loadtxt("./{}/contact_critical_pr_{}.txt".format(folder, folder))
    if folder == "elle":
        analysis = analyses[ii] / 10
        area_ND = data_pr[:, -1] / (pi*0.003**2) * analysis
        pr_ND = (data_pr[:, 1]*(1-0.49**2)/E)*(analysis)**(-2)*6/0.12
    if folder == "dtilde":
        analysis = analyses[ii] / 5
        area_ND = (data_pr[:, -1])/(pi*0.003**2)*(1.1)
        pr_ND = (data_pr[:, 1]*(1-0.49**2)/E)*(1.1**(-2))*(analysis)/0.12
    if folder == "gamma":
        analysis = analyses[ii] / 50
        area_ND = (data_pr[1:, -1])/(pi*0.003**2)*(1.1)
        pr_ND = (data_pr[1:, 1]*(1-0.49**2)/E)*(1.1**(-2))*6/(analysis)

    avg_areas_con.append(mean(area_ND))
    avg_pp_con.append(mean(pr_ND))

area_con = mean(avg_areas_con)
area_con_std = std(avg_areas_con)
pp_con = mean(avg_pp_con)
pp_con_std = std(avg_pp_con)


print("Normalised contact critical area: {} +- {}".format(round(area_con, 3),
                                                          round(area_con_std, 3)))
print("Normalised contact critical pressure: {} +- {}".format(round(pp_con, 3),
                                                           round(pp_con_std,3)))

plot(area_con, pp_con, linestyle="", marker="o", color="black")





