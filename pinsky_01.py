__author__ = 'kuo xiao, derived from Pinsky et al 1994 and Gabbiani prsolve.m'
# here uses the euler method for estimate the solution for equations

from pinsky_functions import *

# set params
e_na = 120
e_k = -15
e_ca = 140
e_leak = 0
g_na_max = 30
g_kdr_max = 15
g_ca_max = 10
g_kc_max = 15
g_kahp_max = 0.8
g_leak_max = 0.1
c_m = 3
p_cable = 0.5

# set more params
i_soma = 1
i_dend = 0
g_cable = 2.1
dt = 0.05

# generate state dict
states = {'t': np.linspace(1, 100, 2000), 'vm_soma': [-4.6], 'vm_dend': [-4.5], 'h': [0.999], 'n': [0.001],
          's': [0.009], 'c': [0.007], 'q': [0.01], 'ca': [0.2]}

for i in range(len(states['t'])):
    # equations for cal the current
    i_leak_s = g_leak_max * (states['vm_soma'][i] - e_leak)
    i_leak_d = g_leak_max * (states['vm_dend'][i] - e_leak)
    m_inf = a_m(states['vm_soma'][i]) / (a_m(states['vm_soma'][i]) + b_m(states['vm_soma'][i]))
    i_na = g_na_max * (m_inf ** 2) * states['h'][i] * (states['vm_soma'][i] - e_na)
    i_kdr = g_kdr_max * states['n'][i] * (states['vm_soma'][i] - e_k)
    i_ca = g_ca_max * (states['s'][i] ** 2) * (states['vm_dend'][i] - e_ca)
    i_kc = g_kc_max * states['c'][i] * states['ca'][i] * (states['vm_dend'][i] - e_k)
    i_kahp = g_kahp_max * states['q'][i] * (states['vm_dend'][i] - e_k)

    # cal the derivatives
    d_vm_soma = (-i_leak_s - i_na - i_kdr + (g_cable / p_cable) *
                 (states['vm_dend'][i] - states['vm_soma'][i]) + i_soma / p_cable) / c_m
    d_vm_dend = (-i_leak_d - i_ca - i_kahp - i_kc +
                 (g_cable / (1-p_cable)) * (states['vm_soma'][i] - states['vm_dend'][i]) +
                 i_dend / (1-p_cable)) / c_m
    d_h = a_h(states['vm_soma'][i]) * (1 - states['h'][i]) - b_h(states['vm_soma'][i]) * states['h'][i]
    d_n = a_n(states['vm_soma'][i]) * (1 - states['n'][i]) - b_n(states['vm_soma'][i]) * states['n'][i]
    d_s = a_s(states['vm_dend'][i]) * (1 - states['s'][i]) - b_s(states['vm_dend'][i]) * states['s'][i]
    d_c = a_c(states['vm_dend'][i]) * (1 - states['c'][i]) - b_c(states['vm_dend'][i]) * states['c'][i]
    d_q = a_q(states['ca'][i]) * (1 - states['q'][i]) - b_q(states['ca'][i]) * states['q'][i]
    d_ca = -0.13 * i_ca - 0.075 * states['ca'][i]

    # update the states
    states['vm_soma'].append(states['vm_soma'][-1] + d_vm_soma * dt)
    states['vm_dend'].append(states['vm_dend'][-1] + d_vm_dend * dt)
    states['h'].append(states['h'][-1] + d_h * dt)
    states['n'].append(states['n'][-1] + d_n * dt)
    states['s'].append(states['s'][-1] + d_s * dt)
    states['c'].append(states['c'][-1] + d_c * dt)
    states['q'].append(states['q'][-1] + d_q * dt)
    states['ca'].append(states['ca'][-1] + d_ca * dt)


plt.figure(dpi=200)
plt.subplot(8, 1, 1)
plt.plot(states['t'], states['vm_soma'][1:])
plt.subplot(8, 1, 2)
plt.plot(states['t'], states['vm_dend'][1:])
plt.subplot(8, 1, 3)
plt.plot(states['t'], states['h'][1:])
plt.subplot(8, 1, 4)
plt.plot(states['t'], states['n'][1:])
plt.subplot(8, 1, 5)
plt.plot(states['t'], states['s'][1:])
plt.subplot(8, 1, 6)
plt.plot(states['t'], states['c'][1:])
plt.subplot(8, 1, 7)
plt.plot(states['t'], states['q'][1:])
plt.subplot(8, 1, 8)
plt.plot(states['t'], states['ca'][1:])
plt.show()
