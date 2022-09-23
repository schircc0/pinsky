import numpy as np
import matplotlib.pyplot as plt


def a_h(v):
    return 0.128 * np.exp((17.-v) / 18.)


def b_h(v):
    return 4. / (1.+np.exp((40.-v) / 5.))


def a_m(v):
    return 0.32 * (13.1-v) / (np.exp((13.1-v) / 4.) - 1.)


def b_m(v):
    return 0.28 * (v-40.1) / (np.exp((v-40.1) / 5.) - 1.)


def a_n(v):
    return 0.016 * (35.1-v) / (np.exp((35.1-v) / 5.) - 1.)


def b_n(v):
    return 0.25 * np.exp(0.5 - 0.025 * v)


def a_s(v):
    return 1.6 / (1. + np.exp(-0.072 * (v-65.)))


def b_s(v):
    return (0.02 * (v-51.1)) / (np.exp((v-51.1) / 5.) - 1.)


def a_c(v):
    if v <= 50:
        return np.exp((v-10.)/11. - (v-6.5)/27.) / 18.975
    else:
        return 2. * np.exp((6.5-v) / 27.)


def b_c(v):
    if v <= 50:
        return 2. * np.exp((6.5-v) / 27.) - np.exp((v-10.)/11. - (v-6.5)/27.) / 18.975
    else:
        return 0


def a_q(v):
    return min((0.00002 * v), 0.01)


def b_q(v):
    return 0.001


def get_diff(i_leak_s, i_leak_d, i_na, i_kdr, i_ca, i_kc, i_kahp, vm_soma, vm_dend, h, n, s, c, q, ca, g_cable, p_cable,
             i_soma, i_dend, c_m):
    d_vm_soma = (-i_leak_s - i_na - i_kdr + (g_cable / p_cable) *
                 (vm_dend - vm_soma) + i_soma / p_cable) / c_m
    d_vm_dend = (-i_leak_d - i_ca - i_kahp - i_kc +
                 (g_cable / (1 - p_cable)) * (vm_soma - vm_dend) +
                 i_dend / (1 - p_cable)) / c_m
    d_h = a_h(vm_soma) * (1 - h) - b_h(vm_soma) * h
    d_n = a_n(vm_soma) * (1 - n) - b_n(vm_soma) * n
    d_s = a_s(vm_dend) * (1 - s) - b_s(vm_dend) * s
    d_c = a_c(vm_dend) * (1 - c) - b_c(vm_dend) * c
    d_q = a_q(ca) * (1 - q) - b_q(ca) * q
    d_ca = -0.13 * i_ca - 0.075 * ca
    return np.array([d_vm_soma, d_vm_dend, d_h, d_n, d_s, d_c, d_q, d_ca])
