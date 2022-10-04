from pyomo.environ import *


def a_h(v):
    return 0.128 * exp((17.-v) / 18.)


def b_h(v):
    return 4. / (1.+exp((40.-v) / 5.))


def a_m(v):
    return 0.32 * (13.1-v) / (exp((13.1-v) / 4.) - 1.)


def b_m(v):
    return 0.28 * (v-40.1) / (exp((v-40.1) / 5.) - 1.)


def a_n(v):
    return 0.016 * (35.1-v) / (exp((35.1-v) / 5.) - 1.)


def b_n(v):
    return 0.25 * exp(0.5 - 0.025 * v)


def a_s(v):
    return 1.6 / (1. + exp(-0.072 * (v-65.)))


def b_s(v):
    return 0.02 * (v-51.1) / (exp((v-51.1) / 5.) - 1.)


# def a_c(v):
#     if v <= 50:
#         return exp((v-10.)/11. - (v-6.5)/27.) / 18.975
#     else:
#         return 2. * exp((6.5-v) / 27.)


# def b_c(v):
#     if v <= 50:
#         return 2. * exp((6.5-v) / 27.) - exp((v-10.)/11. - (v-6.5)/27.) / 18.975
#     else:
#         return 0

def a_c(v):

    return exp((v-10.)/11. - (v-6.5)/27.) / 18.975


def b_c(v):

    return 2. * exp((6.5-v) / 27.) - exp((v-10.)/11. - (v-6.5)/27.) / 18.975



def a_q(v):
    return 0.00002 * v


def b_q(v):
    return 0.001


