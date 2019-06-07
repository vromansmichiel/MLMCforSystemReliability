using JLD2,Plots,Distributions,LaTeXStrings,MultilevelEstimators,Plots.PlotMeasures
p = plot()
nsim=3
namefig="system20_k3/system20_k3.png"
runtimeMC = zeros(10,nsim)
runtimeMLMC = zeros(10,nsim)
runtimeMLQMC = zeros(10,nsim)
tolerance = zeros(10)
α = []
β = []
γ = []
for k = 1:nsim
    @load("system20_k3/system20_k3_$(k)_MLQMC.jld2",history)
    for i = 1:length(history.data)
        for lev = 0:length(history.data[end][:nb_of_samples])-1
            global runtimeMLQMC[i,k]=runtimeMLQMC[i,k]+history.data[i][:T][Level(lev)]*history.data[i][:nb_of_samples][lev+1]
        end
    end
    @load("system20_k3/system20_k3_$(k)_MLMC.jld2",history)
    for i = 1:length(history.data)
        global runtimeMC[i,k]=history.data[i][:nb_of_samples][1]*maximum(history.data[i][:T])[2]
        for lev = 0:length(history.data[end][:nb_of_samples])-1
            global runtimeMLMC[i,k]=runtimeMLMC[i,k]+history.data[i][:T][Level(lev)]*history.data[i][:nb_of_samples][lev+1]
        end
        global tolerance[i]=history.data[i][:tol]
    end
    scatter!(tolerance,runtimeMC[:,k],lw=12,label = "",color=:red,markersize=10);
    scatter!(tolerance,runtimeMLMC[:,k],lw=12,label = "",xaxis=:log,yaxis=:log ,color=:blue,markershape=:diamond,markersize=10);
    scatter!(tolerance,runtimeMLQMC[:,k],lw=12,label = "",xaxis=:log,yaxis=:log ,color=:green,markershape=:star5,markersize=10);
    push!(α,history.data[end][:α])
    push!(β,history.data[end][:β])
    push!(γ,history.data[end][:γ])
    end
    println(" α = ",mean(α),"\n β = ",mean(β),"\n γ = ",mean(γ) )
y = log10.(mean(runtimeMC,dims=2))
x = log10.(tolerance)
a,b = hcat(fill!(similar(x), 1), x) \ y
linr = similar(tolerance)
for i = 1:length(tolerance)
    linr[i] = 10^(b.*log10(tolerance[i]).+a)
end
println("orde MC= ",b)
plot!(tolerance,linr,lw=12,xaxis=:log,yaxis=:log ,color=:red,label = "MC",linestyle=:solid);
y = log10.(mean(runtimeMLMC,dims=2))
x = log10.(tolerance)
a,b = hcat(fill!(similar(x), 1), x) \ y
linr = similar(tolerance)
for i = 1:length(tolerance)
    linr[i] = 10^(b.*log10(tolerance[i]).+a)
end
println("orde MLMC= ",b)
plot!(tolerance,linr,lw=12,xaxis=:log,yaxis=:log ,color=:blue,label = "MLMC",linestyle=:dash);
y = log10.(mean(runtimeMLQMC,dims=2))
x = log10.(tolerance)
a,b = hcat(fill!(similar(x), 1), x) \ y
linr = similar(tolerance)
for i = 1:length(tolerance)
    linr[i] = 10^(b.*log10(tolerance[i]).+a)
end
println("orde MLQMC= ",b)
plot!(tolerance,linr,lw=12,xaxis=:log,yaxis=:log ,color=:green,label = "MLQMC",linestyle=:dot,xlabel=L"tolerantie",ylabel=L"runtime(s)",size=(1600,1100),legendfontsize=28,guidefontsize=36,tickfontsize=20,margin=18mm)
#plot!(tolerance[end-3:end],10^(-4).*tolerance[end-3:end].^(-2),lw=4,color=:black,label="")
display(p);
savefig(p,namefig)
