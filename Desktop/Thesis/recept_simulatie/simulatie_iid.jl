using JLD2,Distributions,Random,Systems,MultilevelEstimators,MultilevelSystems
import MultilevelSystems.level_levensduur,MultilevelSystems.sample_levensduur,MultilevelSystems.sampling_for_QMC,MultilevelSystems.levels_for_QMC
@load "systeem_iid.jld2"
start = time()
T = Component[];
Levels,T,varT = levels_for_QMC(systeem);
lijst = systeem.lijst
level_levensduur(level,ω)=MultilevelSystems.lvl_levensduur(level+one(level),ω,systeem,Levels,lijst)
QMC_ML_levensduur(level,ω)=MultilevelSystems.lvl_levensduur(level+one(level),ω,systeem,Levels,T)
function distribution_QMC_levensduur(G::System,lijst)
    distributions = [MultilevelEstimators.Uniform{Float64}(0.,1.) for i in 1:length(lijst)]
    for i = 1:length(lijst)
    k = lijst[i].shape
    l = lijst[i].scale
    distributions[i] = MultilevelEstimators.Uniform(k,l)
    end
    return distributions
    end
tol=5e-5
elapsed = time() - start;
display(elapsed)
for i = 1:1
Random.seed!(2000+i);
distributions = distribution_QMC_levensduur(G,T);
estimator = MultilevelEstimators.Estimator(ML(), QMC(), QMC_ML_levensduur, distributions, name="systemiid", max_index_set_param=length(Levels)-1,nb_of_warm_up_samples = 100,save_samples=true);#,robustify_bias_estimate = false);
h2 = run(estimator, tol);
println("finished with estimating");
#estimator = MultilevelEstimators.Estimator(ML(), MC(), level_levensduur, distributions, name="systemiid", max_index_set_param=length(Levels)-1,nb_of_warm_up_samples = 100,save_samples=true);#,robustify_bias_estimate = false);
h1 = run(estimator, tol);
println("finished with estimating");
end
systeem = Dict(:systeem => G,:levels => Levels, :lijst => T,:invloed => varT,:setup_time => elapsed)
@save "system20_uniform_system.jld2" systeem
