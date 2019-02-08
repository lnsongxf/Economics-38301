## Crews_Econ38001_StationaryDist

function Crews_Econ38001_StationaryDist(H::Array{Float64,2}, p_LL::Float64,
    p_HH::Float64)
    #=
    The function takes three arguments:
    + a decision rule H
    + transition probabilities p_HH, p_LL

    The function solves for
    + a stationary distribution G
    =#

    # grid size
    N = size(E)[2]

    # Construct indicator matrix
    IndicG = zeros(2,N,2,N)
    for i = 1:2
        for j = 1:N
            for k = 1:2
                for l = 1:N
                    if H[k,l] == K[i,j]
                        IndicG[i,j,k,l] = 1
                    end
                end
            end
        end
    end

    # Construct probability matrix
    Prob = zeros(2,N,2,N)
    for i = 1:2
        for j = 1:N
            for k = 1:2
                for l = 1:N
                    if i == 1 && k == 1
                        Prob[i,j,k,l] = p_LL
                    elseif i == 2 && k == 1
                        Prob[i,j,k,l] = (1.0 - p_LL)
                    elseif i == 1 && k == 2
                        Prob[i,j,k,l] = (1.0 - p_HH)
                    else
                        Prob[i,j,k,l] = p_HH
                    end
                end
            end
        end
    end

    # Reshape
    IndicG = reshape(IndicG, 2*N, 2*N)
    Prob = reshape(Prob, 2*N, 2*N)

    # Multiply elementwise to get =T=
    T = IndicG .* Prob

    # Initialize iteration
    g_0 = repeat([1/(2*N)], 2*N, 1)
    err = 100
    bound = 10^(-15)
    g_n = g_0

    # Iterate to get stationary =g=
    while err > bound
        g_n1 = T * g_n
        err = maximum(abs.(g_n1 - g_n))
        g_n = g_n1
    end

    G = reshape(g_n, 2, N)
    return G
end
