from dolfin import *
import numpy as np 

LENGTH = 1.0
WIDTH = 0.2
HIGHT = 0.1

LAMA_MU = 1.0
LAMA_LAMBDA = 1.25

E = LAMA_MU*(3*LAMA_LAMBDA+2*LAMA_MU)/(LAMA_MU+LAMA_LAMBDA)
nu = LAMA_LAMBDA/2.0/(LAMA_MU+LAMA_LAMBDA)

E1, E2, mu12, G12 = 1.3, 0.77, 0.33, 0.48
# 弹性矩阵
C2D_Aniso = np.array([
    [E/(1-nu**2)   ,   nu*E/(1-nu**2),  0           ],
    [nu*E/(1-nu**2),   E/(1-nu**2)   ,  0           ],
    [0             ,   0             ,  E/2.0/(1+nu)]
])

C2D_Iso = np.array([
    [E1/(1-mu12**2), mu12*E2/(1-mu12**2), 0],
    [mu12*E2/(1-mu12**2), E2/(1-mu12**2), 0],
    [0, 0, G12]
])

# 定义本构
def epsilon(u):
    engineering_strain = 0.5 * (nabla_grad(u) + nabla_grad(u).T)
    return engineering_strain

def sigma(u):
    cauchy_stress = LAMA_LAMBDA * tr(epsilon(u))*Identity(2) + 2 * LAMA_MU * epsilon(u)
    return cauchy_stress

def sigma_tensor(u):
    # 计算应变张量
    epsilon_ij = epsilon(u)
    ep = np.array([epsilon_ij[0,0], epsilon_ij[1,1], epsilon_ij[0,1]])
    # 使用弹性系数计算应力张量
    sigma_ij = np.dot(C2D_Iso, ep)

    return as_tensor([[sigma_ij[0], sigma_ij[2]],
                      [sigma_ij[2], sigma_ij[1]]])
        

# 创建网格并显示
mesh = UnitSquareMesh(20, 20) 
V = VectorFunctionSpace(mesh, "Lagrange", 1)
plot(mesh)
# 狄利克雷边界
def clamped_boundary(x, on_boundary):
    return on_boundary and x[0] < DOLFIN_EPS

bc = DirichletBC(V, Constant((0.0, 0.0)), clamped_boundary)

# 自然边界
class RightEnd(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and abs(x[0] - 1) < DOLFIN_EPS and (abs(x[1] - 0.5) < 0.1)
right_end_boundary = RightEnd()

class TopEnd(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and abs(x[0] - 0.5) < 0.1 and (abs(x[1] - 1) < DOLFIN_EPS)
top_end_boundary = TopEnd()

boundary_mark = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
boundary_mark.set_all(0)
right_end_boundary.mark(boundary_mark, 1)
top_end_boundary.mark(boundary_mark, 2)

# expression = Expression(Constant((0.0, 0.2)), degree=2)
# expression = Constant((0.0, 1.0))
# u_ex = project(expression, V
# plot(u_ex, title="forcing vectors")
u_trial = TrialFunction(V)
v_test = TestFunction(V)
forcing = Constant((50, 0))
forcing2 = Constant((0, 20))
lhs = inner(sigma_tensor(u_trial), nabla_grad(v_test)) * dx
# rhs = dot(forcing, v_test) * dx
rhs = (dot(forcing, v_test) * ds(subdomain_data=boundary_mark, domain=mesh, subdomain_id=1)
       +
       dot(forcing2, v_test) * ds(subdomain_data=boundary_mark, domain=mesh, subdomain_id=2))
# rhs = dot(forcing, v_test) * dx


# plot(rhs)
u_solution = Function(V)

problem = LinearVariationalProblem(lhs, rhs, u_solution, bcs=[bc])
solver = LinearVariationalSolver(problem)
solver.solve()

plot(u_solution,
     title = "my fancy plot")

print("max u: " + str(u_solution.vector().max()))

# print(u_solution.vector()[:])