using JLD2,Distributions,Random,Systems,Reporter,MultilevelEstimators,MultilevelSystems
import MultilevelSystems.level_levensduur,MultilevelSystems.sample_levensduur
include("Toepassingen.jl");
start = time()
G = brandalarm(3.);
print("Generated a system. Number of components: ", G.aantal,"\n");
println("system: ",G.circuit)
T = Component[];
Levels,T,varT = levels_for_QMC(G);
println("Levels en Sampling voor QMC voltooid")
println(T);println(varT);
print("Levels selected! Number of levels: ",length(Levels),"\n");
lijst = G.lijst
level_levensduur(level,ω)=MultilevelSystems.lvl_levensduur(level+one(level),ω,G,Levels,lijst)
QMC_ML_levensduur(level,ω) = MultilevelSystems.QMC_ML_levensduur(level+one(level),ω,G,Levels,T);
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
tol=5e-4
elapsed = time() - start;
display(elapsed)
Random.seed!(2000);
distributions = distribution_QMC_levensduur(G,lijst);
estimator = MultilevelEstimators.Estimator(ML(), QMC(), QMC_ML_levensduur, distributions, name="brandalarm2_MLQMC", max_index_set_param=length(Levels)-1,nb_of_warm_up_samples = 100);#,robustify_bias_estimate = false);
h2 = run(estimator, tol);
println("finished with estimating");
estimator = MultilevelEstimators.Estimator(ML(), MC(), level_levensduur, distributions, name="brandalarm2_MLMC", max_index_set_param=length(Levels)-1,nb_of_warm_up_samples = 100);#,robustify_bias_estimate = false);
#h1 = run(estimator, tol);
println("finished with estimating");
systeem = Dict(:systeem => G,:levels => Levels, :lijst => T,:invloed => varT,:setup_time => elapsed)
@save "brandalarm_system.jld2" systeem
