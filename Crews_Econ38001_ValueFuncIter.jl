## Crews_Econ38001_ValueFuncIter

function Crews_Econ38001_ValueFuncIter(V_n::Array{Float64,2}, Pi::Array{Float64,2},
    E::Array{Float64,2}, K3::Array{Float64,3}, beta::Float64, delta::Float64,
    xi::Float64, thresh::Float64, p_LL::Float64, p_HH::Float64, w::Float64)
    #=
    The function takes eleven arguments:
    + matrix `V_n` that represents the value function (nth iteration)
    + matrix `Pi` of maximized payoffs
    + 2D matrix `E` of epsilons
    + 3D matrix `K` of capital values
    + discount factor `beta`
    + depreciation `delta`
    + adjustment cost `xi`
    + adjustment threshold `thresh`
    + transition probabilities `p_HH`, `p_LL`
    + wage `w`

    The function solves for
    + matrix `V_n1` that gives the next iteration of the value function
    + matrix `H_n1` that gives the policy function
    =#
    N = size(V_n)[2]

    Bell = zeros(2,N,N)
    for i = 1:2
        for j = 1:N
            for k = 1:N
                # ternary operator
                abs(K3[1,1,k] - (1-delta)*K3[1,1,j]) > thresh * K3[1,1,j] ? indicPay = 1 : indicPay = 0
                #(2-i) =1 if i=1, =0 if i=2
                #(i-1) =0 if i=1, =1 if i=2
                contL = beta*(p_LL*V_n[1,k] + (1-p_LL)*V_n[2,k])
                contH = beta*((1-p_HH)*V_n[1,k] + p_HH*V_n[2,k])
                cont = (2-i)*contL + (i-1)*contH
                Bell[i,j,k] = Pi[i,j] - (K3[1,1,k] - (1-delta)*K3[1,1,j]) - indicPay*xi*w + cont
            end
        end
    end

    V, I = findmax(Bell, dims=3) # [max values, Cartesian index of max]
    V_n1 = reshape(V, 2, N)
    H_n1 = dropdims(K3[I], dims=3)
    return V_n1, H_n1
end
