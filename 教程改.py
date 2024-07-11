import fenics as fe

CANTILEVER_LENGTH = 1.0
CANTILEVER_WIDTH = 0.2

N_POINTS_LENGTH = 10
N_POINTS_WIDTH = 3

LAME_MU = 1.0
LAME_LAMBDA = 1.25
DENSITY = 1.0
ACCELERATION_DUE_TO_GRAVITY = 0.016

def main():
    # Mesh and Vector Function Space
    mesh = fe.BoxMesh(
        fe.Point(0.0, 0.0, 0.0),
        fe.Point(CANTILEVER_LENGTH, CANTILEVER_WIDTH, CANTILEVER_WIDTH),
        N_POINTS_LENGTH,
        N_POINTS_WIDTH,
        N_POINTS_WIDTH,
    )
    lagrange_vector_space_first_order = fe.VectorFunctionSpace(
        mesh,
        "Lagrange",
        1,
    )
    
    # Boundary Conditions
    def clamped_boundary(x, on_boundary):
        return on_boundary and x[0] < fe.DOLFIN_EPS

    dirichlet_clamped_boundary = fe.DirichletBC(
        lagrange_vector_space_first_order,
        fe.Constant((0.0, 0.0, 0.0)),
        clamped_boundary,
    )

    # Define strain and stress
    def epsilon(u):
        engineering_strain = 0.5 * (fe.nabla_grad(u) + fe.nabla_grad(u).T)
        return engineering_strain
    
    def sigma(u):
        cauchy_stress = (
            LAME_LAMBDA * fe.tr(epsilon(u)) * fe.Identity(3)
            +
            2 * LAME_MU * epsilon(u)
        )
        return cauchy_stress
    
    # Define weak form
    u_trial = fe.TrialFunction(lagrange_vector_space_first_order)
    v_test = fe.TestFunction(lagrange_vector_space_first_order)
    forcing = fe.Constant((0.0, 0.0, - DENSITY * ACCELERATION_DUE_TO_GRAVITY))
    traction = fe.Constant((0.0, 0.0, 0.0))

    weak_form_lhs = fe.inner(sigma(u_trial), epsilon(v_test)) * fe.dx  # Crucial to use inner and not dot
    weak_form_rhs = (
        fe.dot(forcing, v_test) * fe.dx
        +
        fe.dot(traction, v_test) * fe.ds
    )

    # Compute solution
    u_solution = fe.Function(lagrange_vector_space_first_order)
    fe.solve(
        weak_form_lhs == weak_form_rhs,
        u_solution,
        dirichlet_clamped_boundary,
    )
    u_vector = u_solution.vector()[:]
    u_mesh = u_solution.function_space().mesh()
    import matplotlib.figure as plt
    plt.plot(u_solution)
    u_solution.function_space()
    # u_solution()[:]
    # Compute the von Mises stress
    deviatoric_stress_tensor = (
        sigma(u_solution)
        -
        1/3 * fe.tr(sigma(u_solution)) * fe.Identity(3)
    )
    von_Mises_stress = fe.sqrt(3/2 * fe.inner(deviatoric_stress_tensor, deviatoric_stress_tensor))

    lagrange_scalar_space_first_order = fe.FunctionSpace(
        mesh,
        "Lagrange",
        1,
    )
    von_Mises_stress = fe.project(von_Mises_stress, lagrange_scalar_space_first_order)

    # Write out fields for visualization with Paraview
    u_solution.rename("Displacement Vector", "")
    von_Mises_stress.rename("von Mises stress", "")

    beam_deflection_file = fe.XDMFFile("beam_deflection.xdmf")
    beam_deflection_file.parameters["flush_output"] = True
    beam_deflection_file.parameters["functions_share_mesh"] = True
    beam_deflection_file.write(u_solution, 0.0)
    beam_deflection_file.write(von_Mises_stress, 0.0)
    print(u_solution)

if __name__ == "__main__":
    main()