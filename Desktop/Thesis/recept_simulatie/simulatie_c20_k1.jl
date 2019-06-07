using JLD2,Systems,MultilevelEstimators,MultilevelSystems
import MultilevelSystems.level_levensduur,MultilevelSystems.sample_levensduur

## Voorbereiding
start = time()
# willekeurig systeem
G = generateSystem("G",20,0.4,0.5,1.);
println("Generated a system. Number of components: ", G.aantal);
println("system: ",G.circuit)
T = Component[];
# bepalen levelstructuur en ordening
Levels,T,varT = levels_for_QMC(G);
println("Levels en Sampling voor QMC voltooid")
println(T);println(varT);
println("Aantal levels: ",length(Levels));
lijst = G.lijst
# simulatiefunctie
level_levensduur(level,W)=MultilevelSystems.lvl_levensduur(level
						+one(level),W,G,Levels,lijst)
levelQ_levensduur(level,W) = MultilevelSystems.lvl_levensduur(level
						+one(level),W,G,Levels,T)
# verdelingen van de componenten
function distribution_levensduur(G::System,lijst)
    distributions = [MultilevelEstimators.Weibull()
    							for i in 1:length(lijst)]
    for i = 1:length(lijst)
    k = lijst[i].shape
    l = lijst[i].scale
     distributions[i] = MultilevelEstimators.Weibull(k,l)
    end
    return distributions
end
# bepaling duur van de set-up
elapsed = time() - start;
display(elapsed)
# nauwkeurigheid
tol=1e-3
id = 1

## MLMC berekening
distributions = distribution_levensduur(G,lijst);
distributionsQ = distribution_levensduur(G,T);
estimator = MultilevelEstimators.Estimator(ML(), QMC(),
   levelQ_levensduur, distributionsQ,
   name="levensduur_$(id)_MLQMC", max_index_set_param=
   length(Levels)-1, nb_of_warm_up_samples = 100);
h2 = run(estimator, tol);

## MLQMC berekening
estimator = MultilevelEstimators.Estimator(ML(), MC(),
   level_levensduur, distributions,
   name="levensduur_$(id)_MLMC", max_index_set_param=
   length(Levels)-1, nb_of_warm_up_samples = 100);
h2 = run(estimator, tol);

## opslaan van het systeem en extra informatie
systeem = Dict(:systeem => G,:levels => Levels, :lijst => T,
                 :invloed => varT,:setup_time => elapsed)
@save "system_$(id).jld2" systeem
