## Crews_Econ38001_DecisionRules

function Crews_Econ38001_DecisionRules(w::Float64, E::Array{Float64,2}, K::Array{Float64,2},
    theta::Float64, nu::Float64, beta::Float64, delta::Float64, xi::Float64, thresh::Float64,
    p_LL::Float64, p_HH::Float64)
    #=
    The function takes ten arguments:
    + initial wage guess `w`
    + 2D matrix `E` of epsilons
    + 2D matrix `K` of capital values
    + capital elasticity `theta`
    + labor demand elasticity `nu`
    + discount factor `beta`
    + depreciation `delta`
    + adjustment cost `xi`
    + adjustment threshold `thresh`
    + transition probabilities `p_LL`, `p_HH`

    The function solves for
    + matrix `H` that gives the policy function
    =#

    # grid size
    N = size(E)[2]

    # Compute Pi
    Pi = ((nu/w) * E .* K .^(theta)).^(nu/(1-nu)) .* (E .* K .^(theta)) - w*((nu/w) * E .* K .^(theta)).^(1/(1-nu))

    # Because I'll use Cartesian indices, the following matrix will be helpful
    K3 = zeros(2, N, N)
    for k = 1:N
        K3[:,k,:] = K
    end

    # Initialize =while= loop
    err = 100
    bound = 10^(-6)
    V_1 = Crews_Econ38001_ValueFuncIter(Pi, Pi, E, K3, beta, delta, xi, thresh, p_LL, p_HH, w)[1]

    # =while= loop
    V_n1 = V_1
    while err > bound
        V_n2 = Crews_Econ38001_ValueFuncIter(V_n1, Pi, E, K3, beta, delta, xi, thresh, p_LL, p_HH, w)[1]
        V_n = V_n1
        V_n1 = V_n2
        err = maximum(abs.(V_n1 - V_n))
    end

    # Get decision rule
    H = Crews_Econ38001_ValueFuncIter(V_n1, Pi, E, K3, beta, delta, xi, thresh, p_LL, p_HH, w)[2]
    return H, p_LL, p_HH
end
