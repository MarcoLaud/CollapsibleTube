# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 13:45:32 2022

@author: marco
"""

from numpy import asarray, cos, diag, exp, genfromtxt, linspace, loadtxt, mean,\
                  sin, sqrt, tanh, zeros
from matplotlib.pyplot import errorbar, figure, fill_between, legend, plot,\
                              text, title, xlabel, xlim, ylabel, ylim
from scipy.optimize import curve_fit
from pathlib import Path
from os.path import basename


def order_par(xx, A, c, p_cr, beta):
    return A*tanh(c*(p_cr/xx - 1)**beta)


def linear_function(xx, m, q):
    return m*xx + q


def zerofive(xx, yy):
    return len(xx)*[mean(yy)]


def tanh_gen(xx, a, b, c):
    return a*tanh(b*xx) + c

def sigmoid(xx, a, b, c):
    return a/(1+exp(-b*xx)) + c


# case = 11
# p_crit = 2450
# p_min = 1000

# case = 12
# p_crit = 3100
# p_min = 1500

case = 13
p_crit = 3400
p_min = 1900

#case = 14
#p_crit = 3500
#p_min = 2000

#case = 15
#p_crit = 3550
#p_min = 1000

#case = 16
#p_crit = 3550
#p_min = 800

# case = 17
# p_crit = 3550
# p_min = 1400

# case = 18
# p_crit = 3550
# p_min = 1500


figure()
title("Order parameter")
xlabel("Pressure [Pa]")
ylabel("Area - A_buckling [m^2]")


pathlist = Path('./l={}/'.format(str(case))).glob('rtz*.csv')
path_list = [str(element) for element in pathlist]
sort_path_list = sorted(path_list,
                        key=lambda i: int(
                                basename(i).split('_')[1].split('.')[0]))


pp = genfromtxt("./pressure.csv", delimiter=",", skip_header=1)
pp = pp[:, 1] * (-1)
PP = pp * (-1)

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

    areas[jj] = area

press = []
idx = []

for kk, elem in enumerate(PP):
    if elem > p_min and elem < p_crit:
        idx.append(kk)
        press.append(elem)

aa = areas[idx]
aa = aa - min(aa)

par, cov = curve_fit(order_par, press, aa, p0=[max(aa), 1, p_crit, 0.5])

plot(press, aa, label="simulation data")
plot(press, order_par(press, par[0], par[1], par[2], par[3]), linestyle="--",
     label="a*tanh(b*(p_cr/p - 1)^c")
legend()
#
print(par)
print(sqrt(diag(cov)))

# Fitting critical exponents
data_exp = loadtxt("critical_exp.txt")
ps = data_exp[:, 0] / 10.
crit_exp = data_exp[:, 1]
errs = data_exp[:, 2]

par_lin, cov_lin = curve_fit(linear_function, ps, crit_exp, p0=[0, 0.5])

figure()
title("Critical exponents dependence on pre-stretch ratio")
xlabel("Pre-stretch ratio [-]")
ylabel("Critical exponent [-]")
ylim(0, 1)
errorbar(ps, crit_exp, yerr=errs, linestyle="", marker="o")
plot(ps, zerofive(ps, crit_exp), linestyle="--",
     label="Mean value = {}".format(round(mean(crit_exp),2)))
legend()


# Fitting critical pressures
data_pr = loadtxt("critical_pr.txt")
crit_pr = data_pr[:, 1]
errs_pr = data_pr[:, 2]
par_tan, cov_tan = curve_fit(tanh_gen, ps, crit_pr, p0=[max(crit_pr), 1, -1e6])
XX = linspace(min(ps), max(ps), 50)

figure()
title("Pre-stretch - Critical pressures")
xlabel("Pre-stretch ratio [-]")
ylabel("Critical pressure [Pa]")
errorbar(ps, crit_pr, linestyle="", marker="o")
plot(XX, tanh_gen(XX, par_tan[0], par_tan[1], par_tan[2]), linestyle="--",
     label="a*tanh(b*l)+c")
legend()

xx = linspace(1, 2, 50)
ymin = min(tanh_gen(xx, par_tan[0], par_tan[1], par_tan[2]))
ymax = max(tanh_gen(xx, par_tan[0], par_tan[1], par_tan[2]))

figure()
title("Pre-stretch - Phase diagram")
xlabel("Pre-stretch ratio [-]")
ylabel("Pressure [Pa]")
xlim(1, max(xx))
ylim(ymin, ymax + 500)
plot(xx, tanh_gen(xx, par_tan[0], par_tan[1], par_tan[2]), color="black",
     label="Critical pressure")
fill_between(xx, ymin, tanh_gen(xx, par_tan[0], par_tan[1], par_tan[2]))
fill_between(xx, tanh_gen(xx, par_tan[0], par_tan[1], par_tan[2]), ymax + 500,
             alpha=0.5)
text(1.05, 3500, "Buckling state")
text(1.3, 3000, "Non-buckling state")
legend()
