# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 13:15:32 2022

@author: marco
"""

from numpy import asarray, cos, genfromtxt, mean, pi, ones, sin, std, zeros
from matplotlib.pyplot import errorbar, figure, legend, plot, title, xlabel,\
    ylabel
from pathlib import Path
from os.path import basename


import matplotlib.pyplot as plt
#tube_number = 3    # questo funziona abbastanza col tubo 3
#ts = 0.1
#critical_time = 20
#extr_pr = 4600

tube_number = 8
#ts = 0.1
#critical_time = 250
#extr_pr = 4000

data_p_exp = genfromtxt('ptm.csv', delimiter=',', skip_header=1)
err_p_exp = genfromtxt('ptm_pm.csv', delimiter=',', skip_header=1)
data_A_exp = genfromtxt('area.csv', delimiter=',', skip_header=1)

data_p_exp = data_p_exp[tube_number-1, :] * (-100)  # P<0 convention, mbar->Pa
err_p_exp = err_p_exp[tube_number-1, :] * 100       # mbar -> Pa
data_A_exp = data_A_exp[tube_number-1, :] * 1e-6    # mm^2 -> m^2

diff_area = max(abs(data_A_exp[1] - data_A_exp[0]),
                abs(data_A_exp[2] - data_A_exp[1]),
                abs(data_A_exp[2] - data_A_exp[0]))

#A_err = ones(len(err_p_exp)) * std(data_A_exp[:3])
A_err = ones(len(err_p_exp)) * diff_area

geom_param = genfromtxt('summary.csv', delimiter=',', skip_header=1)
llzero = round(geom_param[tube_number-1, 1]/geom_param[tube_number-1, 2], 2)

pathlist = Path('./tube{}'.format(tube_number)).glob('rtz*.csv')
path_list = [str(element) for element in pathlist]
sort_path_list = sorted(path_list,
                        key=lambda i: int(
                                basename(i).split('_')[1].split('.')[0]))


fig = figure(figsize=(8, 5))
#title('Tube law - Pressure ramp comparison')
xlabel('Normalised area [-]')
ylabel('Pressure [Pa]')

ax = fig.add_subplot(111)
ax.spines["bottom"].set_visible(False)
ax.spines["right"].set_visible(False)
#ax.spines["left"].set_visible(False)

ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
ax.yaxis.tick_left()
ax.yaxis.set_label_position('left')

#ax.set_xlim()

#pp = [((-ii*ts)/critical_time) * extr_pr for ii in range(len(sort_path_list))]
pp = genfromtxt("./tube{}/pressure.csv".format(tube_number),
                delimiter=",", skip_header=1)
pp = pp[:, 1] * (-1)
areas = zeros(len(sort_path_list))
JJ = []

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

    print('Deformed area - Tube {} : {} m^2'.format(tube_number, area))
    areas[jj] = area

l, caps, c = errorbar((data_A_exp * 1e6)/(pi*3.15**2), data_p_exp,
                      xerr=(A_err*1e6)/(pi*3.15**2), yerr=err_p_exp,
                      elinewidth=2, markeredgewidth=2,
                      linestyle='', marker='o', color="black",
                      label='Gregory et al. [14]')

plot((areas * 1e6)/(pi*3.15**2), pp, label='Simulation data', linestyle="-.",
     color="grey", linewidth=2)

#plt.scatter(areas[JJ]*1e6/(pi*3.15**2), pp[JJ], color="red")
legend()

#fig.savefig("tube8.svg", format="svg", dpi=1200)
