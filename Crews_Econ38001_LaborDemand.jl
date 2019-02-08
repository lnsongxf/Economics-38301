## Crews_Econ38001_LaborDemand

function Crews_Econ38001_LaborDemand(w::Float64, E::Array{Float64,2}, K::Array{Float64,2},
    theta::Float64, nu::Float64, beta::Float64, delta::Float64, xi::Float64, p_LL::Float64,
    p_HH::Float64)
    #=
    The function takes ten arguments:
    + initial wage guess w
    + 2D matrix E of epsilons
    + 2D matrix K of capital values
    + capital elasticity theta
    + labor demand elasticity nu
    + discount factor beta
    + depreciation delta
    + adjustment cost xi
    + transition probabilities p_HH, p_LL

    The function solves for
    + labor demand N_d that matches the guessed wage
    =#

    H, p_LL, p_HH = Crews_Econ38001_DecisionRules(w, E, K, theta, nu, beta, delta, xi, p_LL, p_HH)
    G = Crews_Econ38001_StationaryDist(H, p_LL, p_HH)

    # firm labor demand function
    n_d = ((nu/w) * E .* K .^(theta)).^(1/(1-nu))
    N_temp = n_d * G'
    N_d = N_temp[1,1] + N_temp[2,2]
    return N_d
end
