module Metodos
    using SparseArrays
    export calc_x, dif_finitas

    function calc_x(delta_x, L)
        x = Float64[]

        for i in 0:delta_x:L
            append!(x, i)
        end

        return x
    end

    function gauss_seidel_prof(_A::Array{Float64,2},_b::Array{Float64}, _tol::Float64,_nit::Int64)
        m,n = size(_A)
        x = zeros(Float64, n, 1)
        C = copy(_A)
        d = copy(_b)
    
        for i=1:n
            C[i,i] = 0.0
            d[i,1] /= _A[i,i]
            
            for j=1:n
                C[i,j] /= _A[i,i]
            end
        end
        
        iter = 0
        eps = typemax(Float64)
        xold = 0.0
        while (iter < _nit) && (eps > _tol)
            eps = 0.0
            for i=1:n
                
                xold = x[i,1]
                x[i,1] = d[i,1] - (C[i,:]'*x)[1,1]
                eps = max(abs((x[i,1]-xold)/x[i,1]), eps)
            end
            iter += 1
        end
    
        return x[:,1]
    
    end

    function gauss_seidel(_A::Array{Float64,2}, _b::Array{Float64}, _tol::Float64, _nit::Int64)
        m,n = size(_A)
        if m != n
            return nothing
        end
        x = zeros(Float64, n)
        C = copy(_A)
        d = copy(_b)
    
        for i=1:n
            C[i,i] = 0.0
            d[i,1] /= _A[i,i]
            
            for j=1:n
                C[i,j] /= _A[i,i]
            end
        end
        
        iter = 0
        eps = typemax(Float64)
        xold = 0.0
        while (iter < _nit) && (eps > _tol)
            eps = 0.0
            for i=1:n
                xold = x[i]
                x[i] = d[i] - (C[i,:]'*x)
                eps = max(abs((x[i]-xold)/x[i]), eps)
            end
            iter += 1
        end
    
        return x
    
    end

    function gauss_seidel_3_vetores(A_D::Array{Float64}, A_DP::Array{Float64}, A_E::Array{Float64}, b::Array{Float64}, x::Array{Float64}, _tol=nothing, _nit=nothing)
        if _tol === nothing
            _tol = 0.0001
        end
        if _nit === nothing
            _nit = 100
        end

        N = length(A_DP)

        b ./= A_DP # divide b pela diagonal principal
        A_D ./= A_DP
        A_E ./= A_DP

        iter = 0
        eps = typemax(Float64)
        xold = 0.0
        println(x, "\n")
        while (iter < _nit) && (eps > _tol)
            eps = 0.0
            for i=1:N
                xold = x[i]

                if i == 1
                    x[i] = b[i] - (A_D[i]*x[i+1])
                elseif i == N
                    x[i] = b[i] - (A_E[i]*x[i-1])
                else
                    x[i] = b[i] - (A_E[i]*x[i-1]) - (A_D[i]*x[i+1])
                end
                
                eps = max(abs((x[i]-xold)/x[i]), eps)
            end

            iter += 1
        end

    end

    function gauss_seidel_bloco(bloco::Array{Float64}, b::Array{Float64})
    end

    function dif_finitas(pvc, method::String, _tol=nothing, _nit=nothing)
        if !(method in ["sparse_direto", "direto", "3_vetores", "bloco"])
            return nothing
        end

        N = length(pvc.x)
        b = fill(pvc.h_linha*(pvc.delta_x^2)*pvc.t_inf, N-2)
        b[1] += pvc.t_a # primeira linha do sistema
        b[N-2] += pvc.t_b # ultima linha do sistema
        
        if method == "sparse_direto" || method == "direto"
            if occursin("sparse", method)
                A = spzeros(N-2, N-2)
            else
                A = zeros(Float64, N-2, N-2)
            end

            # primeiro
            A[1, 1] = 2 + pvc.h_linha*pvc.delta_x^2
            A[1, 2] = -1

            # ultimo
            A[N-2, N-3] = -1
            A[N-2, N-2] = 2 + pvc.h_linha*pvc.delta_x^2
            
            for i in 2:N-3
                A[i, i-1] = -1
                A[i, i] = (2 + pvc.h_linha * pvc.delta_x^2)
                A[i, i+1] = -1

                #b[i] = pvc.h_linha*(pvc.delta_x^2)*pvc.t_inf
            end

            #display(A)
            #print("\n")

            #result = A\b

            result = gauss_seidel(A, b, 0.001, 100)
        
        elseif method == "3_vetores"
            A_D = fill(-1.0, N-3)
            insert!(A_D, N-2, 0) # o vetor a direita da diagonal nao chega ate a ultima linha
            A_DP = fill(2 + pvc.h_linha * pvc.delta_x^2, N-2) # Diagonal Principal
            A_E = fill(-1.0, N-3)
            insert!(A_E, 1, 0) # o vetor a esquerda da diagonal nao chega na primeira linha

            result = zeros(Float64, N-2)

            gauss_seidel_3_vetores(A_D, A_DP, A_E, b, result, _tol, _nit)

        elseif method == "bloco"
            bloco = [-1, 2 + pvc.h_linha * pvc.delta_x^2, -1]

            result = gauss_seidel_bloco(bloco, b)
        end

        insert!(result, 1, pvc.t_a)
        insert!(result, N, pvc.t_b)

        return result
    end
end