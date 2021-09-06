function solver(_A::Array{Float64,2},_b::Array{Float64,1}, _tol::Float64,_nit::Int64)
    m,n = size(_A)
    x = zeros(Float64, n, 1)
    C = copy(_A)
    d = copy(b)

    for i=1:n
        C[i,i] = 0.0
        d[i,i] /= A[i,i]
        
        for j=1:n
            C[i,j] = A[i,i]
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
            eps = max(abs(x[i,1]-xold)/x[i,1], eps)
        end
        iter += 1
    end

    return x

end