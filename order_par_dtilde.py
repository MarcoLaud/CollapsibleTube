# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 13:45:32 2022

@author: marco
"""

from numpy import asarray, cos, diag, genfromtxt, linspace, loadtxt, mean,\
                  pi, sin, sqrt, tanh, zeros
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


def amico_fun(xx, a, b, c):
    return a*(xx)**(b) + c


#def vonmises(E,nu,gamma,d):
#    return (E/(4*(1-nu)))*gamma**3 +\
#            ((pi**2)/(pi**2 + 2*d**2))*\
#            (((E*(7-nu))/(12*(1-nu**2)))*gamma**3+((E*gamma)/(3)))

def vonmises(E,nu,gamma,d):
    return ((E * 3)/(12 * (1-nu**2)))*(2*gamma)**3+\
    (1 / (1 + ((2/pi)*2*d)**2))*\
    (((E/12) * ((7-nu)/(1-nu**2)) * (2*gamma)**3) +
     ((E/3) * (2*gamma)))


def zarandi(a, d, b, c):
    return 98.06 * (a*d**b + c)

# case = 30
# p_crit = 2450
# p_min = 1200

# case = 35
# p_crit = 1900
# p_min = 1000

# case = 40
# p_crit = 1500
# p_min = 600

# case = 45
# p_crit = 1300
# p_min = 600

# case = 50
# p_crit = 1100
# p_min = 600

# case = 55
# p_crit = 1000
# p_min = 500

case = 60
p_crit = 900
p_min = 300


figure()
title("Order parameter")
xlabel("Pressure [Pa]")
ylabel("Area - A_buckling [m^2]")


pathlist = Path('./d={}/'.format(str(case))).glob('rtz*.csv')
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

print(par)
print(sqrt(diag(cov)))

# Linear fit of the critical exponents
data_exp = loadtxt("critical_exp.txt")
ps = data_exp[:, 0] / 10.
crit_exp = data_exp[:, 1]
errs = data_exp[:, 2]

par_lin, cov_lin = curve_fit(linear_function, ps, crit_exp, p0=[0, 0.5])

figure()
title("Critical exponents dependence on length-diameter ratio")
xlabel("Length-diameter ratio [-]")
ylabel("Critical exponent [-]")
ylim(0, 1)
errorbar(ps, crit_exp, yerr=errs, linestyle="", marker="o")
plot(ps, zerofive(ps, crit_exp), linestyle="--",
     label="Mean value = {}".format(round(mean(crit_exp),2)))
legend()

# Fit of the critical pressures
data_pr = loadtxt("critical_pr.txt")
crit_pr = data_pr[:, 1]
errs_pr = data_pr[:, 2]
par_am, cov_am = curve_fit(amico_fun, ps, crit_pr, p0=[1e3, -2, 3000])
XX = linspace(min(ps), max(ps), 50)

figure()
title("Length-diameter ratio - Critical pressures")
xlabel("Length-diameter ratio [-]")
ylabel("Critical pressure [Pa]")
errorbar(ps, crit_pr, linestyle="", marker="o")
plot(XX, amico_fun(XX, par_am[0], par_am[1], par_am[2]), linestyle="--",
     label="a*d^(b) + c")
legend()

xx = linspace(1, 6, 50)
ymin = min(amico_fun(xx, par_am[0], par_am[1], par_am[2]))
ymax = max(amico_fun(xx, par_am[0], par_am[1], par_am[2]))

figure()
title("Length-diameter ratio - Phase diagram")
xlabel("Length-diameter ratio [-]")
ylabel("Pressure [Pa]")
xlim(min(xx), max(xx))
ylim(ymin, ymax + 500)
plot(xx, amico_fun(xx, par_am[0], par_am[1], par_am[2]), color="black",
     label="Critical pressure")
fill_between(xx, ymin, amico_fun(xx, par_am[0], par_am[1], par_am[2]))
fill_between(xx, amico_fun(xx, par_am[0], par_am[1], par_am[2]), ymax + 500,
             alpha=0.5)
text(2.5, 7500, "Buckling state")
text(1.1, 2700, "Non-buckling state")
legend()

# Comparison with von Mises
gamma = 0.06    # h/d, for von mises we need to to l/r = 2 x l/d
EE = 1e6  # Pa
NU = 0.49

figure()
title("Length-diameter ratio - Phase diagram")
xlabel("Length-diameter ratio [-]")
ylabel("Critical pressure [Pa]")
errorbar(ps, crit_pr, linestyle="", marker="o", label="Simulation data")
plot(XX, amico_fun(XX, par_am[0], par_am[1], par_am[2]), linestyle="--",
     label="a*d^(b) + c")
plot(XX, vonmises(EE, NU, gamma, XX), label="von Mises")
plot(XX, zarandi(150, XX, -3.3, 3.230), label="zarandi")
legend()
