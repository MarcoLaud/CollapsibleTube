# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 13:15:32 2022

@author: marco
"""

from numpy import asarray, cos, genfromtxt, mean, pi, sin, std, zeros
from matplotlib.pyplot import errorbar, figure, legend, plot, title, xlabel,\
    ylabel
from pathlib import Path
from os.path import basename


# tube_number = 1

# data_p_exp = genfromtxt('ptm.csv', delimiter=',', skip_header=1)
# err_p_exp = genfromtxt('ptm_pm.csv', delimiter=',', skip_header=1)
# data_A_exp = genfromtxt('area.csv', delimiter=',', skip_header=1)

# data_p_exp = data_p_exp[tube_number-1, :] * (-100)  # P<0 convention, mbar->Pa
# err_p_exp = err_p_exp[tube_number-1, :] * 100       # mbar -> Pa
# data_A_exp = data_A_exp[tube_number-1, :] * 1e-6    # mm^2 -> m^2

# A_err = ones(len(err_p_exp)) * std(data_A_exp[:3])


analysis = ["400Pas", "800Pas", "1600Pas"]
plot_style = ["--", ":", "-."]


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

# errorbar((data_A_exp * 1e6)/(pi*3.15**2), data_p_exp,
#          xerr=(A_err*1e6)/(pi*3.15**2), yerr=err_p_exp,
#          linestyle='', marker='o', label='Exp data', color="black")


for ps, case in enumerate(analysis):

    pathlist = Path('./{}/'.format(case)).glob('rtz*.csv')
    path_list = [str(element) for element in pathlist]
    sort_path_list = sorted(path_list,
                            key=lambda i: int(
                                    basename(i).split('_')[1].split('.')[0]))


    pp = genfromtxt("./pressure_{}.csv".format(case), delimiter=",", skip_header=1)
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

        print('Deformed area - {} : {} m^2'.format(case, area))
        areas[jj] = area
        # cross_check = areas

    # for jj in range(len(cross_check)-1):
    #     if abs(cross_check[jj] - cross_check[jj+1]) > 6e-6:
    #         ind = [kk for kk in range(jj, jj+8)]
    #         areas = delete(areas, ind)
    #         pp = delete(pp, ind)
    #         print("Cross-check revealed wrong computations! Points removed.")

    ll = case.split("P")[0] + " Pa/s"
    plot((areas)/(areas[0]), pp, label=ll, linestyle=plot_style[ps],
         linewidth=2)


legend()

# fig.savefig("ramp.svg", format="svg", dpi=1200)
