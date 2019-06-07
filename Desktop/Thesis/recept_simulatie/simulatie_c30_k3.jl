using JLD2,Distributions,Random,Systems,MultilevelEstimators,MultilevelSystems
import MultilevelSystems.level_levensduur,MultilevelSystems.sample_levensduur
include("Statische_systemen/ToyExamples.jl");
start = time()
G = creeerToyExample(3,ncomp=30);
print("Generated a system. Number of components: ", G.aantal,"\n");
println("system: ",G.circuit)
T = Component[];
Levels,T,varT = levels_for_QMC(G);
println("Levels en Sampling voor QMC voltooid")
println(T);println(varT);
print("Levels selected! Number of levels: ",length(Levels),"\n");
lijst = G.lijst
level_levensduur(level,ω)=MultilevelSystems.lvl_levensduur(level+one(level),ω,G,Levels,lijst)
QMC_ML_levensduur(level,ω) = MultilevelSystems.QMC_ML_levensduur(level+one(level),ω,G,Levels,T)
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
tol=1e-3
elapsed = time() - start;
display(elapsed)
for i = 1:3
    Random.seed!(2000+i);
    distributions = distribution_QMC_levensduur(G,T);
    estimator = MultilevelEstimators.Estimator(ML(), MC(), level_levensduur, distributions, name="system30_k3_$(i)_MLMC", max_index_set_param=length(Levels)-1,nb_of_warm_up_samples = 100);#,robustify_bias_estimate = false);
    h1 = run(estimator, tol);
    println("finished with estimating");
    estimator = MultilevelEstimators.Estimator(ML(), QMC(), QMC_ML_levensduur, distributions, name="system30_k3_$(i)_MLQMC", max_index_set_param=length(Levels)-1,nb_of_warm_up_samples = 100);#,robustify_bias_estimate = false);
    h2 = run(estimator, tol);
    println("finished with estimating");
end
systeem = Dict(:systeem => G,:levels => Levels, :lijst => T,:invloed => varT,:setup_time => elapsed)
@save "system30_k3_system.jld2" systeem
