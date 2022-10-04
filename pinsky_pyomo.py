from pyomo.environ import *
from pyomo.dae import *
from pinsky_functions_pyomo import *
import numpy as np

i_soma = 1
i_dend = 0
g_cable = 2.1
c_m = 3
p_cable = 0.5

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


m = AbstractModel()
m.t = ContinuousSet()
m.MEAS_t = Set(within=m.t)
m.vm_soma_meas = Param(m.MEAS_t)
m.vm_dend_meas = Param(m.MEAS_t)
m.h_meas = Param(m.MEAS_t)
m.n_meas = Param(m.MEAS_t)
m.s_meas = Param(m.MEAS_t)
m.c_meas = Param(m.MEAS_t)
m.q_meas = Param(m.MEAS_t)
m.ca_meas = Param(m.MEAS_t)

m.vm_soma = Var(m.t, within=PositiveReals, initialize=-4.6)
m.vm_dend = Var(m.t, within=PositiveReals, initialize=-4.5)
m.h = Var(m.t, within=PositiveReals, initialize=0.999)
m.n = Var(m.t, within=PositiveReals, initialize=0.001)
m.s = Var(m.t, within=PositiveReals, initialize=0.009)
m.c = Var(m.t, within=PositiveReals, initialize=0.007)
m.q = Var(m.t, within=PositiveReals, initialize=0.01)
m.ca = Var(m.t, within=PositiveReals, initialize=0.2)

m.i_leak_s = Var(m.t, within=PositiveReals)
m.i_leak_d = Var(m.t, within=PositiveReals)
m.m_inf = Var(m.t, within=PositiveReals)
m.i_na = Var(m.t, within=PositiveReals)
m.i_kdr = Var(m.t, within=PositiveReals)
m.i_ca = Var(m.t, within=PositiveReals)
m.i_kc = Var(m.t, within=PositiveReals)
m.i_kahp = Var(m.t, within=PositiveReals)

m.vm_soma_d = DerivativeVar(m.vm_soma, wrt=m.t)
m.vm_dend_d = DerivativeVar(m.vm_dend, wrt=m.t)
m.h_d = DerivativeVar(m.h, wrt=m.t)
m.n_d = DerivativeVar(m.n, wrt=m.t)
m.s_d = DerivativeVar(m.s, wrt=m.t)
m.c_d = DerivativeVar(m.c, wrt=m.t)
m.q_d = DerivativeVar(m.q, wrt=m.t)
m.ca_d = DerivativeVar(m.ca, wrt=m.t)


def vm_soma_d_rule(m, i):
    if i == 0:
        return Constraint.Skip
    return m.vm_soma_d[i] == (-m.i_leak_s[i] - m.i_na[i] - m.i_kdr[i] + (g_cable / p_cable) *
                 (m.vm_dend[i] - m.vm_soma[i]) + i_soma / p_cable) / c_m


def vm_dend_d_rule(m, i):
    if i == 0:
        return Constraint.Skip
    return m.vm_dend_d[i] == (-m.i_leak_d[i] - m.i_ca[i] - m.i_kahp[i] - m.i_kc[i] +
                 (g_cable / (1-p_cable)) * (m.vm_soma[i] - m.vm_dend[i]) +
                 i_dend / (1-p_cable)) / c_m


def h_d_rule(m, i):
    if i == 0:
        return Constraint.Skip
    return m.h_d[i] == a_h(m.vm_soma[i]) * (1 - m.h[i]) - b_h(m.vm_soma[i]) * m.h[i]


def n_d_rule(m, i):
    if i == 0:
        return Constraint.Skip
    return m.n_d[i] == a_n(m.vm_soma[i]) * (1 - m.n[i]) - b_n(m.vm_soma[i]) * m.n[i]


def s_d_rule(m, i):
    if i == 0:
        return Constraint.Skip
    return m.s_d[i] == a_s(m.vm_dend[i]) * (1 - m.s[i]) - b_s(m.vm_dend[i]) * m.s[i]


def c_d_rule(m, i):
    if i == 0:
        return Constraint.Skip
    return m.c_d[i] == a_c(m.vm_dend[i]) * (1 - m.c[i]) - b_c(m.vm_dend[i]) * m.c[i]


def q_d_rule(m, i):
    if i == 0:
        return Constraint.Skip
    return m.q_d[i] == a_q(m.ca[i]) * (1 - m.q[i]) - b_q(m.ca[i]) * m.q[i]


def ca_d_rule(m, i):
    if i == 0:
        return Constraint.Skip
    return m.ca_d[i] == -0.13 * m.i_ca[i] - 0.075 * m.ca[i]


def i_leak_s_rule(m, i):
    if i == 0:
        return Constraint.Skip
    return m.i_leak_s[i] == g_leak_max * (m.vm_soma[i] - e_leak)


def i_leak_d_rule(m, i):
    if i == 0:
        return Constraint.Skip
    return m.i_leak_d[i] == g_leak_max * (m.vm_dend[i] - e_leak)


def m_inf_rule(m, i):
    if i == 0:
        return Constraint.Skip
    return m.m_inf[i] == a_m(m.vm_soma[i]) / (a_m(m.vm_soma[i]) + b_m(m.vm_soma[i]))


def i_na_rule(m, i):
    if i == 0:
        return Constraint.Skip
    return m.i_na[i] == g_na_max * (m.m_inf[i] ** 2) * m.h[i] * (m.vm_soma[i] - e_na)


def i_kdr_rule(m, i):
    if i == 0:
        return Constraint.Skip
    return m.i_kdr[i] == g_kdr_max * m.n[i] * (m.vm_soma[i] - e_k)


def i_ca_rule(m, i):
    if i == 0:
        return Constraint.Skip
    return m.i_ca[i] == g_ca_max * (m.s[i] ** 2) * (m.vm_dend[i] - e_ca)


def i_kc_rule(m, i):
    if i == 0:
        return Constraint.Skip
    return m.i_kc[i] == g_kc_max * m.c[i] * (m.ca[i] / 250) * (m.vm_dend[i] - e_k)


def i_kahp_rule(m, i):
    if i == 0:
        return Constraint.Skip
    return m.i_kahp[i] == g_kahp_max * m.q[i] * (m.vm_dend[i] - e_k)


def _obj(m):
    return sum((m.vm_soma[i]-m.vm_soma_meas[i])**2+(m.vm_dend[i]-m.vm_dend_meas[i])**2 for i in m.MEAS_t)


m.vm_soma_d_const = Constraint(m.t, rule=vm_soma_d_rule)
m.vm_dend_d_const = Constraint(m.t, rule=vm_dend_d_rule)
m.h_d_const = Constraint(m.t, rule=h_d_rule)
m.n_d_const = Constraint(m.t, rule=n_d_rule)
m.s_d_const = Constraint(m.t, rule=s_d_rule)
m.c_d_const = Constraint(m.t, rule=c_d_rule)
m.q_d_const = Constraint(m.t, rule=q_d_rule)
m.ca_d_const = Constraint(m.t, rule=ca_d_rule)
m.i_leak_s_const = Constraint(m.t, rule=i_leak_s_rule)
m.i_leak_d_const = Constraint(m.t, rule=i_leak_d_rule)
m.m_inf_const = Constraint(m.t, rule=m_inf_rule)
m.i_na_const = Constraint(m.t, rule=i_na_rule)
m.i_kdr_const = Constraint(m.t, rule=i_kdr_rule)
m.i_ca_const = Constraint(m.t, rule=i_ca_rule)
m.i_kc_const = Constraint(m.t, rule=i_kc_rule)
m.i_kahp_const = Constraint(m.t, rule=i_kahp_rule)
m.obj = Objective(rule=_obj)


m.pprint()

instance = m.create_instance('config/exp.dat')
instance.t.pprint()

discretizer = TransformationFactory('dae.collocation')
discretizer.apply_to(instance, nfe=30)#,ncp=3)

solver = SolverFactory('ipopt')

results = solver.solve(instance, tee=True)
