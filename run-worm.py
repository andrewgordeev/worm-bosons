#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" main.cpp compiled to /main before running this script"""

import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

""" Project 1, part 2 """

# """ Parameters """ 
# num_events = 100000
# num_burn = 10000
# N = 2
# mu = 1.4
# beta = 12
# t = 1
# epsilon_list = [0.03,0.02,0.01,0.008,0.005,0.002,0.001]
# energy_list = np.array([])
# number_list = np.array([])

# """ Running worm code and saving results"""
# for epsilon in epsilon_list:
#     process = subprocess.run(['./main', str(num_events), str(num_burn), str(N), str(epsilon), str(mu), str(beta), str(t)],  
#                capture_output=True, text=True)
#     lines = process.stdout.splitlines()
#     energy_list = np.append(energy_list, float(lines[-3][8:lines[-3].find('±')]))
#     number_list = np.append(number_list, float(lines[-2][8:lines[-2].find('±')]))
    
# """ Polynomial fit """
# energy_fit = np.polyfit(epsilon_list, energy_list, deg=3)
# number_fit = np.polyfit(epsilon_list, number_list, deg=3)

# print("<n> =", number_fit[-1])
# print("<e> =", energy_fit[-1])

""" Project 2, part 1 """

# """ Parameters """ 
# num_events = 100000
# num_burn = 10000
# N = 4
# beta = 100
# t = 1
# epsilon = 0.01
# mu_range = np.linspace(0, 0.5, 51)
# number_list = np.array([])
# number_err_list = np.array([])

# """ Running worm code and saving results"""
# for mu in mu_range:
#     print(mu)
#     process = subprocess.run(['./main', str(num_events), str(num_burn), str(N), str(epsilon), str(mu), str(beta), str(t)],  
#                 capture_output=True, text=True)
#     lines = process.stdout.splitlines()
#     number_list = np.append(number_list, float(lines[-2][8:lines[-2].find('±')]))
#     number_err_list = np.append(number_err_list, float(lines[-2][lines[-2].find('±')+1:]))
    
# """ Plotting results """
# plt.errorbar(mu_range, number_list, yerr=number_err_list, fmt = 'kx-')
# plt.xlabel(r'$\mu$')
# plt.ylabel(r'$<n>$')

""" Project 2, part 2 """

# """ Parameters """ 
# num_events = 100000
# num_burn = 10000
# N = 20
# mu = 1.0
# t = 1
# epsilon = 0.01
# beta_range = np.linspace(0.75, 1.1, 8)
# chi_list = np.array([])
# chi_err_list = np.array([])

# """ Running worm code and saving results"""
# for beta in beta_range:
#     print(beta)
#     process = subprocess.run(['./main', str(num_events), str(num_burn), str(N), str(epsilon), str(mu), str(beta), str(t)],  
#                 capture_output=True, text=True)
#     lines = process.stdout.splitlines()
#     chi_list = np.append(chi_list, float(lines[-1][15:lines[-1].find('±')]))
#     chi_err_list = np.append(chi_err_list, float(lines[-1][lines[-1].find('±')+1:]))
    
# """ Plotting results """
# plt.errorbar(beta_range, chi_list, yerr=chi_err_list, fmt = 'kx-')
# plt.xlabel(r'$\beta$')
# plt.ylabel(r'$\chi_\omega$')

""" Project 2, part 3 """

""" Parameters """ 
num_events = 100000
num_burn = 10000
N_list = [16, 20, 24]
mu = 1.0
t = 1
epsilon = 0.01
beta_range = np.linspace(0.88, 0.95, 8)
number_list = [np.array([]),np.array([]),np.array([])]
number_err_list = [np.array([]),np.array([]),np.array([])]
chi_list = [np.array([]),np.array([]),np.array([])]
chi_err_list = [np.array([]),np.array([]),np.array([])]

""" Running worm code and saving results"""
for i in range(len(N_list)):
    for beta in beta_range:
        print(N_list[i], beta)
        process = subprocess.run(['./main', str(num_events), str(num_burn), str(N_list[i]), str(epsilon), str(mu), str(beta), str(t)],  
                                  capture_output=True, text=True)
        lines = process.stdout.splitlines()
        number_list[i] = np.append(number_list[i], float(lines[-2][8:lines[-2].find('±')]))
        number_err_list[i] = np.append(number_err_list[i], float(lines[-2][lines[-2].find('±')+1:]))
        chi_list[i] = np.append(chi_list[i], float(lines[-1][15:lines[-1].find('±')]))
        chi_err_list[i] = np.append(chi_err_list[i], float(lines[-1][lines[-1].find('±')+1:]))
    
""" Plotting results """
for i in range(len(N_list)):
    plt.errorbar(beta_range, chi_list[i], yerr=chi_err_list[i], fmt = 'x-', label="N = " + str(N_list[i]))
plt.xlabel(r'$\beta$')
plt.ylabel(r'$N \chi_\omega$')
plt.legend()

""" Calculating intersection - only for 3 entries in N_list """
# Linear interpolating splines for each susceptibility
susc1 = UnivariateSpline(beta_range, chi_list[0], k=1, s=0)
susc2 = UnivariateSpline(beta_range, chi_list[1], k=1, s=0)
susc3 = UnivariateSpline(beta_range, chi_list[2], k=1, s=0)

# Finding 3 intersections
beta_space = np.linspace(0.88, 0.95, 1000000)
beta1 = beta_space[np.argwhere(np.diff(np.sign(N_list[1]*susc2(beta_space) - N_list[0]*susc1(beta_space)))).flatten()]
beta2 = beta_space[np.argwhere(np.diff(np.sign(N_list[2]*susc3(beta_space) - N_list[0]*susc1(beta_space)))).flatten()]
beta3 = beta_space[np.argwhere(np.diff(np.sign(N_list[1]*susc2(beta_space) - N_list[2]*susc3(beta_space)))).flatten()]
betas = np.array([beta1, beta2, beta3])

beta_crit = np.mean(betas)
beta_crit_err = np.std(betas)

# Calculating particle density
nums1 = UnivariateSpline(beta_range, number_list[0], k=1, s=0)
nums2 = UnivariateSpline(beta_range, number_list[1], k=1, s=0)
nums3 = UnivariateSpline(beta_range, number_list[2], k=1, s=0)

rho = np.zeros(3)
rho[0] = nums1(beta_crit)/N_list[0]
rho[1] = nums2(beta_crit)/N_list[1]
rho[2] = nums3(beta_crit)/N_list[2]

rho_avg = np.mean(rho)
rho_err = np.std(rho)