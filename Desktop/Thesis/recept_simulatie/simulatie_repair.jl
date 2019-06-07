using JLD2,Distributions,Random,Systems,MultilevelEstimators,MultilevelSystems
import MultilevelSystems.level_levensduur,MultilevelSystems.sample_levensduur,MultilevelSystems.QMC_ML_levensduur
include("Statische_systemen/repairProcess.jl")
include("Statische_systemen/ToyExamples.jl");
start = time()
#G = creeerToyExample(3,ncomp=30);
#G = generateSystem("G",25,0.6,0.3,1.)
print("Generated a system. Number of components: ", G.aantal,"\n");
println("system: ",G.circuit)
T = Component[];
Levels,T,varT = levels_for_QMC(G);
println("Levels en Sampling voor QMC voltooid")
println(T);println(varT);
print("Levels selected! Number of levels: ",length(Levels),"\n");
lijst = G.lijst
function repair_levensduur(level::Level,ω::Vector{<:Real},G::System,Levels,lijst::Vector{<:System})
     # finest grid
     #lijst = component_list(G);
     sampleLevensduur(lijst,ω);
         Qf = levensduur_repair(Levels[level]);
         dQ = Qf
     if level != Level(1)
       Qf,dQ = levensduur_repair(Levels[level],Levels[level-one(level)]);
     end
     return dQ, Qf
 end
level_levensduur(level,ω)=repair_levensduur(level+one(level),ω,G,Levels,lijst)
QMC_ML_levensduur(level,ω)=repair_QMC_levensduur(level+one(level),ω,G,Levels,lijst)
function repair_QMC_levensduur(level::Level,ω::Vector{<:Real},G::System,Levels,lijst)
    # finest grid
    sampleLevensduur(lijst,ω);
    Qf = levensduur_repair(Levels[level]);
    # compute difference when not on coarsest grid
    dQ = Qf
    if level != Level(1)
      Qc = levensduur_repair(Levels[level-one(level)]);
      dQ -= Qc
    end
    return dQ, Qf
end
function distribution_QMC_levensduur(G::System,lijst)
    distributions = [MultilevelEstimators.Weibull() for i in 1:length(lijst)]
    for i = 1:length(lijst)
    k = lijst[i].shape
    l = lijst[i].scale
     distributions[i] = MultilevelEstimators.Weibull(k,l)
    end
    return distributions
end
tol=1e-3
elapsed = time() - start;
display(elapsed)

    Random.seed!(2000);
    distributions = distribution_QMC_levensduur(G,T);
    estimator = MultilevelEstimators.Estimator(ML(), MC(), level_levensduur, distributions, name="system25_Brepairable_MLMC", max_index_set_param=length(Levels)-1,nb_of_warm_up_samples = 100);#,robustify_bias_estimate = false);
    h1 = run(estimator, tol);
    #println("finished with estimating");
    #estimator = MultilevelEstimators.Estimator(ML(), QMC(), QMC_ML_levensduur, distributions, name="system_repairable_MLQMC", max_index_set_param=length(Levels)-1,nb_of_warm_up_samples = 100);#,robustify_bias_estimate = false);
    #h2 = run(estimator, tol);
    println("finished with estimating");

systeem = Dict(:systeem => G,:levels => Levels, :lijst => T,:invloed => varT,:setup_time => elapsed)
@save "system25_system.jld2" systeem
