module PVCClass
    export PVC, pvc_solution

    struct PVC
        h_linha::Float64
        h::Float64
        r::Float64
        k::Float64
        t_inf::Float64
        t_a::Float64
        t_b::Float64
        L::Float64
        delta_x::Float64
        x::Array
    end

    # Solucao analitica
    function pvc_solution(pvc::PVC)
        lambda = sqrt(pvc.h_linha)

        A = ((pvc.t_a - pvc.t_inf)*MathConstants.e^(-lambda*pvc.L) - (pvc.t_b - pvc.t_inf)) / (MathConstants.e^(-lambda*pvc.L) - MathConstants.e^(lambda*pvc.L))
        B = ((pvc.t_b - pvc.t_inf) - (pvc.t_a - pvc.t_inf)*MathConstants.e^(lambda*pvc.L)) / (MathConstants.e^(-lambda*pvc.L) - MathConstants.e^(lambda*pvc.L))

        N = length(pvc.x)
        y = Array{Float64, 1}(undef, N)
        for i in 2:N-1
            y[i] = pvc.t_inf + A * MathConstants.e^(lambda*pvc.x[i]) + B * MathConstants.e^(-lambda*pvc.x[i])
        end
        y[1] = pvc.t_a
        y[N] = pvc.t_b
        return y
    end
end