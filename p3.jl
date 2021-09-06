using Plots
using JSON
include("./metodos.jl")
import .Metodos
include("./pvc.jl")
using .PVCClass

function readJSON(_filename::String)
    data = Dict()
    data = JSON.parsefile(_filename)
    return data
end

function delta_x_plot(pvc, d_out::Dict, delta_x::Float64)
    d_out["delta_x = $delta_x"] = Dict()

    d_out["delta_x = $delta_x"]["Solucao Analítica"] = pvc_solution(pvc)
    #func = plot(pvc.x, d_out["delta_x = $delta_x"]["Solucao Analítica"], marker=:circle, label="S. Analítica")

    d_out["delta_x = $delta_x"]["Diferenças Finitas"] = Metodos.dif_finitas(pvc, "3_vetores")
    #plot!(pvc.x, d_out["delta_x = $delta_x"]["Diferenças Finitas"], marker=:circle, label="Dif. Finitas", xlabel="Posição da Barra", ylabel="Temperatura")

    #filename = string("pvc_deltaX_", delta_x, ".png")
    #savefig(func, filename)
end

function main(_filename = nothing)
    d_out = Dict()
    if _filename === nothing
        # Problema teste
        delta_x = 2.0
        x = Metodos.calc_x(delta_x, 10)
        pvc = PVC(0.05, 1, 0.2, 200, 200, 243, 400, 10, delta_x, x)
        delta_x_plot(pvc, d_out, delta_x)
        #=
        delta_x = 1.0
        x = Metodos.calc_x(delta_x, 10)
        pvc = PVC(0.05, 1, 0.2, 200, 200, 243, 400, 10, delta_x, x)
        delta_x_plot(pvc, d_out, delta_x)

        delta_x = 0.5
        x = Metodos.calc_x(delta_x, 10)
        pvc = PVC(0.05, 1, 0.2, 200, 200, 243, 400, 10, delta_x, x)
        delta_x_plot(pvc, d_out, delta_x)

        delta_x = 0.25
        x = Metodos.calc_x(delta_x, 10)
        pvc = PVC(0.05, 1, 0.2, 200, 200, 243, 400, 10, delta_x, x)
        delta_x_plot(pvc, d_out, delta_x)
        =#
    else
        data = readJSON(_filename)
        pvc = PVC(data["h_linha"], data["h"], data["r"], data["k"], data["t_inf"], data["t_a"], data["t_b"], data["L"], data["delta_x"], data["x"])
        d_out["Solucao Analítica"] = pvc_solution(pvc)
        d_out["Diferenças Finitas"] = Metodos.dif_finitas(pvc)
    end

    fout = open("output.json", "w")
    JSON.print(fout, sort(d_out, rev=true), 2)
    close(fout)
end


if length(ARGS) == 1
    main(ARGS[1])
else
    main()
end
