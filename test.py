from dolfin import *
import numpy as np 

LENGTH = 1.0
WIDTH = 0.2
HEIGHT = 0.1

LAMA_MU = 1.0
LAMA_LAMBDA = 1.25

def epsilon(u):
    engineering_strain = 0.5 * (nabla_grad(u) + nabla_grad(u).T)
    return engineering_strain

def sigma(u):
    cauchy_stress = LAMA_LAMBDA * tr(epsilon(u))*Identity(2) + 2 * LAMA_MU * epsilon(u)
    return cauchy_stress

mesh = UnitSquareMesh(10,10)
V = VectorFunctionSpace(mesh, "Lagrange", 1)
# plot(mesh)

# Boundary Conditions
def clamped_boundary(x, on_boundary):
    return on_boundary and x[0] < DOLFIN_EPS

db = DirichletBC(
    V,
    Constant((0.0, 0.0)),
    clamped_boundary,
)

u_trial = TrialFunction(V)
v_test = TestFunction(V)

forcing = Constant((1, 1))

lhs = inner(sigma(u_trial), nabla_grad(v_test)) * dx
rhs = dot(forcing, v_test) * dx

u_solution = Function(V)

problem = LinearVariationalProblem(lhs, rhs, u_solution, bcs=db)

solver = LinearVariationalSolver(problem)
solver.solve()

print("pause..................................")
print(u_solution.vector()[:])


