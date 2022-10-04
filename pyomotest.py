import pyomo.environ as pyo

model = pyo.ConcreteModel()
model.x = pyo.Var([1, 2], domain=pyo.NonNegativeReals)
model.OBJ = pyo.Objective(expr=2 * model.x[1] + 3 * model.x[2], sense=pyo.minimize)
model.Constraint1 = pyo.Constraint(expr=3 * model.x[1] + 4 * model.x[2] >= 1)
solver = pyo.SolverFactory('glpk')
solver.solve(model)
model.pprint()
print(model.OBJ())
print(model.solutions.solutions)
for i in model.x:
    print(model.x[i].value)
for v in model.component_data_objects(pyo.Var):
    print(str(v), v.value)
print('1')

# model2 = pyo.AbstractModel()
# model2.m = pyo.Param(within=pyo.NonNegativeIntegers)
# model2.n = pyo.Param(within=pyo.NonNegativeIntegers)
#
# model2.I = pyo.RangeSet(1, model2.m)
# model2.J = pyo.RangeSet(1, model2.n)
#
# model2.a = pyo.Param(model2.I, model2.J)
# model2.b = pyo.Param(model2.I)
# model2.c = pyo.Param(model2.J)
#
# model2.x = pyo.Var(model2.J, domain=pyo.NonNegativeReals)
#
#
# def obj_expression(m):
#     return pyo.summation(m.c, m.x)
#
#
# model2.obj = pyo.Objective(rule=obj_expression)
#
#
# def ax_constraint_rule(m,i):
#     return sum(m.a[i,j] * m.x[j] for j in m.J) >= m.b[i]
#
#
# model2.AxbConstraint = pyo.Constraint(model2.I, rule=ax_constraint_rule)