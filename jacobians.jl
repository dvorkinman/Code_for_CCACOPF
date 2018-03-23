using JuMP

# extracts a block of the jacobian from the JuMP model m
# at the current solution point.
# rows correspond to those indexed by nonlinear_constraints
# columns correspond to those indexed by variables
function extract_jacobian(m, nonlinear_constraints, variables)
    n_linconstr = MathProgBase.numlinconstr(m)
    n_quadconstr = MathProgBase.numquadconstr(m)
    n_allconstr = MathProgBase.numconstr(m)
    n_var = MathProgBase.numvar(m)

    # nonlinear constraints appear after the linear and quadratic constraints
    d = JuMP.NLPEvaluator(m)
    MathProgBase.initialize(d, [:Jac])

    I,J = MathProgBase.jac_structure(d)

    current_solution = [getvalue(Variable(m,i)) for i in 1:n_var]

    V = zeros(length(I))
    MathProgBase.eval_jac_g(d,V,current_solution)

    jac = sparse(I,J,V,n_allconstr,n_var)

    g = zeros(n_allconstr)
    MathProgBase.eval_g(d,g,current_solution)

    rows = [linearindex(ref)+n_linconstr+n_quadconstr for ref in nonlinear_constraints]
    cols = [linearindex(v) for v in variables]
    L,U = JuMP.constraintbounds(m)
    @assert all(L[rows] .== 0)
    @assert all(U[rows] .== 0)

    return g[rows], jac[rows,cols]
end
