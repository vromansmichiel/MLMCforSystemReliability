using JLD2,Distributions,Random,Systems,Reporter,MultilevelEstimators,MultilevelSystems,DelimitedFiles
import MultilevelSystems.level_levensduur,MultilevelSystems.QMC_ML_levensduur,MultilevelSystems.sample_levensduur
include("Toepassingen.jl");
include("dep.jl")
start = time()
G,D = hydrocentrale();
print("Generated a system. Number of components: ", G.aantal,"\n");
println("system: ",G.circuit)
T = Component[];
Levels,T,varT = levels_for_QMC(G);
println("Levels en Sampling voor QMC voltooid")
println(T);println(varT);
print("Levels selected! Number of levels: ",length(Levels),"\n");
lijst = G.lijst
function depend_levensduur(level::Level,ω::Vector{<:Real},G::System,Levels,lijst::Vector{<:System},D)
     # finest grid
     #lijst = component_list(G);
     sampleLevensduur(lijst,ω);
     dependence(G,D);
     Qf = levensduur(Levels[level]);
     # compute difference when not on coarsest grid
     dQ = Qf
     if level != Level(1)
       #sampleLevensduur(G);
       Qc = levensduur(Levels[level-one(level)]);
       dQ -= Qc
     end
     return dQ, Qf
     end
level_levensduur(level,ω)=depend_levensduur(level+one(level),ω,G,Levels,lijst,D)
QMC_ML_levensduur(level,ω) = depend_levensduur(level+one(level),ω,G,Levels,T,D);
SL_levensduur(level,ω) = QMC_SL_levensduur(level+one(level),ω,G,Levels,T);
function distribution_QMC_levensduur(G::System,lijst)
    distributions = [MultilevelEstimators.Weibull() for i in 1:length(lijst)]
    for i = 1:length(lijst)
    k = lijst[i].shape
    l = lijst[i].scale
     distributions[i] = MultilevelEstimators.Weibull(k,l)
    end
    return distributions
end
elapsed = time() - start;
display(elapsed)
tol = 2e-4
Random.seed!(2000);
distributions = distribution_QMC_levensduur(G,T);
estimator = MultilevelEstimators.Estimator(ML(), QMC(), QMC_ML_levensduur, distributions, name="hydro3_MLQMC", max_index_set_param=length(Levels)-1,nb_of_warm_up_samples = 1000);#,robustify_bias_estimate = false);
h2 = run(estimator, tol);
# save samples
S = MultilevelEstimators.samples_diff(estimator)
name = estimator[:name]
sf_dir = string(name[1:end-5], "_samples")
isdir(sf_dir) || mkdir(sf_dir)
for idx in CartesianIndices(S)
    idx_dir = joinpath(sf_dir, join(idx.I,"_"))
    isdir(idx_dir) || mkdir(idx_dir)
    for k in keys(S[idx])
        writedlm(joinpath(idx_dir, string("samples_level_", k[1]-1, ".txt")), S[idx][k])
    end
end
println("finished with estimating");
estimator = MultilevelEstimators.Estimator(ML(), MC(), level_levensduur, distributions, name="hydro3_MLMC", max_index_set_param=length(Levels)-1,nb_of_warm_up_samples = 1000);#,robustify_bias_estimate = false);
h1 = run(estimator, tol);
# save samples
S = MultilevelEstimators.samples_diff(estimator)
name = estimator[:name]
sf_dir = string(name[1:end-5], "_samples")
isdir(sf_dir) || mkdir(sf_dir)
for idx in CartesianIndices(S)
    idx_dir = joinpath(sf_dir, join(idx.I,"_"))
    isdir(idx_dir) || mkdir(idx_dir)
    for k in keys(S[idx])
        writedlm(joinpath(idx_dir, string("samples_level_", k[1]-1, ".txt")), S[idx][k])
    end
end
println("finished with estimating");
systeem = Dict(:systeem => G,:levels => Levels, :lijst => T,:invloed => varT,:setup_time => elapsed)
@save "hydro_system.jld2" systeem
