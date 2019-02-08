## Crews_Econ38001_LaborDemand

function Crews_Econ38001_LaborDemand(w::Float64, E::Array{Float64,2}, K::Array{Float64,2},
    theta::Float64, nu::Float64, beta::Float64, delta::Float64, xi::Float64, thresh::Float64,
    p_LL::Float64, p_HH::Float64)
    #=
    The function takes eleven arguments:
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
    + labor demand `N_d` that matches the guessed wage
    =#
    N = size(E)[2]

    H, p_LL, p_HH = Crews_Econ38001_DecisionRules(w, E, K, theta, nu, beta, delta, xi, thresh, p_LL, p_HH)
    G = Crews_Econ38001_StationaryDist(H, p_LL, p_HH)

    # firm labor demand function
    n_d = ((nu/w) * E .* K .^(theta)).^(1/(1-nu))
    XI = zeros(2,N)
    for i = 1:2
        for j = 1:N
            if abs(H[i,j] - (1-delta)*K[i,j]) > thresh * K[i,j]
                XI[i,j] = xi
            end
        end
    end
    n_d = n_d + XI

    N_temp = n_d * G'
    N_d = N_temp[1,1] + N_temp[2,2]
    return N_d
end
