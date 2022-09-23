__author__ = 'kuo xiao, derived from Pinsky et al 1994 and Gabbiani prsolve.m'
# here uses the runge-kutta 4 method for estimate the solution for equations

from pinsky_functions import *

# set params
e_na = 120.
e_k = -15.
e_ca = 140.
e_leak = 0.
g_na_max = 30.
g_kdr_max = 15.
g_ca_max = 10.
g_kc_max = 15.
g_kahp_max = 0.8
g_leak_max = 0.1
c_m = 3.
p_cable = 0.5

# set more params
i_soma = 5.
i_dend = 0.
g_cable = 2.1
dt = 0.05

# generate state dict
states = {'t': np.linspace(1, 1000, 20000), 'vm_soma': [-4.6], 'vm_dend': [-4.5], 'h': [0.999], 'n': [0.001],
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

    # cal the next step
    y_k1 = np.array([states['vm_soma'][i], states['vm_dend'][i], states['h'][i], states['n'][i], states['s'][i], states['c'][i],
           states['q'][i], states['ca'][i]])
    d_y_k1 = get_diff(i_leak_s, i_leak_d, i_na, i_kdr, i_ca, i_kc, i_kahp, y_k1[0], y_k1[1], y_k1[2], y_k1[3], y_k1[4],
                     y_k1[5], y_k1[6], y_k1[7], g_cable, p_cable,i_soma, i_dend, c_m)
    y_k2 = y_k1 + (d_y_k1 * dt / 2)
    d_y_k2 = get_diff(i_leak_s, i_leak_d, i_na, i_kdr, i_ca, i_kc, i_kahp, y_k2[0], y_k2[1], y_k2[2], y_k2[3], y_k2[4],
                      y_k2[5], y_k2[6], y_k2[7], g_cable, p_cable,i_soma, i_dend, c_m)
    y_k3 = y_k1 + (d_y_k2 * dt / 2)
    d_y_k3 = get_diff(i_leak_s, i_leak_d, i_na, i_kdr, i_ca, i_kc, i_kahp, y_k3[0], y_k3[1], y_k3[2], y_k3[3], y_k3[4],
                      y_k3[5], y_k3[6], y_k3[7], g_cable, p_cable,i_soma, i_dend, c_m)
    y_k4 = y_k1 + (d_y_k3 * dt)
    d_y_k4 = get_diff(i_leak_s, i_leak_d, i_na, i_kdr, i_ca, i_kc, i_kahp, y_k4[0], y_k4[1], y_k4[2], y_k4[3], y_k4[4],
                      y_k4[5], y_k4[6], y_k4[7], g_cable, p_cable,i_soma, i_dend, c_m)

    y_final = y_k1 + 1/6 * (d_y_k1 + 2*d_y_k2 + 3*d_y_k3 + d_y_k4) * dt

    # update the states
    states['vm_soma'].append(y_final[0])
    states['vm_dend'].append(y_final[1])
    states['h'].append(y_final[2])
    states['n'].append(y_final[3])
    states['s'].append(y_final[4])
    states['c'].append(y_final[5])
    states['q'].append(y_final[6])
    states['ca'].append(y_final[7])


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
