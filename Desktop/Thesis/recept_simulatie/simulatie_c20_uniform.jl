using JLD2,Distributions,Random,Systems,MultilevelEstimators,MultilevelSystems
import MultilevelSystems.level_levensduur,MultilevelSystems.sample_levensduur,MultilevelSystems.sampling_for_QMC,MultilevelSystems.levels_for_QMC
include("Statische_systemen/ToyExamples.jl");
start = time()
function sampling_for_QMC(G::System);
    T = [];
    ETj = [];
    sample_size = 50;
    list = G.lijst;
    ix = [i for i = 1:length(list)];
    A = reshape(Float64[],length(G.lijst),0);
    B = reshape(Float64[],length(G.lijst),0);
    cut,Q = vindCutset(G,eigenschap="systeem");
    for i = 1:sample_size
        ω = [rand(Distributions.Uniform(G.lijst[i].shape,G.lijst[i].scale)) for i = 1:length(G.lijst)]
        A = hcat(A,ω);
        B = hcat(B,[rand(Distributions.Uniform(G.lijst[i].shape,G.lijst[i].scale)) for i = 1:length(G.lijst)]);
        sampleLevensduur(list,ω);
        push!(T,levensduur(cut));
    end
    println("Levensduur is bepaald")
    for k = 1:length(list)
        Tk = []
        for j = 1:sample_size
        ω = A[:,j]
        ω[k] = B[k,j];
        sampleLevensduur(list,ω);
        fab = levensduur(cut);
        push!(Tk,(T[j]-fab)^2);
        end
        push!(ETj,mean(Tk));
        println("variantie component $k is bepaald")
    end
    lijst = sortslices(hcat(ETj,ix),dims=1,rev=true)
    return list[lijst[:,2]];
 end
function levels_for_QMC(G::System);
    T = [];
    ETj = [];
    sample_size = 50;
    list = G.lijst;
    index = [i for i = 1:length(list)];
    A = reshape(Float64[],length(G.lijst),0);
    B = reshape(Float64[],length(G.lijst),0);
    cut,Q = vindCutset(G,minimal = true, eigenschap="systeem");
    L = Int(max(floor(log2(Q)-2),3));
    nr = Int(ceil(Q/(2^L)));
    pilot = reshape(Float64[],Q,0);
    for i = 1:sample_size
        ω = [rand(Distributions.Uniform(G.lijst[i].shape,G.lijst[i].scale)) for i = 1:length(G.lijst)]
        A = hcat(A,ω);
        B = hcat(B,[rand(Distributions.Uniform(G.lijst[i].shape,G.lijst[i].scale)) for i = 1:length(G.lijst)]);
        sampleLevensduur(list,ω);
        pilot=hcat(pilot,levensduur(cut,extend=true));
    end
     T = minimum(pilot,dims=1);
    println("Levensduur is bepaald")
    l = mean(pilot,dims=2); # level 0
    lst = sort(l,dims=1)[nr];
    ix = findall(l.<=lst);
    level_prev = cut[ix];
    ix = findall(l.>lst);
    ctrial = cut[ix];
    T_prev,Tvar = MC_levensduur(G,level_prev,100);
    levels = [level_prev];
    for l = 1:L # level 1 tot L
      l==L ? nr = size(cut,1)-size(level_prev,1) : nr = Int(ceil(Q/2^(L-l)))-Int(ceil(Q/2^(L-(l-1))))
      trial = T_prev .-min.(T_prev,pilot[ix,:])
      delta = mean(trial,dims=2);
      delta = hcat(delta,1:length(delta));
      delta = sortslices(delta,dims=1,by = x->x[1],rev = true);
      ix=convert(Array{Int,1},delta[1:nr,2])
      level = ctrial[ix];
      level = vcat(level_prev,level);
      T_prev,Tvar = MC_levensduur(G,level_prev,100);
      ix = delta[nr+1:end,2];
      ix=convert(Array{Int,1},ix)
      ctrial = ctrial[ix];
      push!(levels,level);
      level_prev = level;
    end
    println("Levels selected")
    for k = 1:length(list)
        Tk = []
        for j = 1:sample_size
        ω = A[:,j]
        ω[k] = B[k,j];
        sampleLevensduur(list,ω);
        fab = levensduur(cut);
        push!(Tk,(T[j]-fab)^2);
        end
        push!(ETj,mean(Tk));
        println("variantie component $k is bepaald")
    end
    lijst = sortslices(hcat(ETj,index),dims=1,rev=true)
    return levels,list[lijst[:,2]],lijst[:,1];
end
    ncomp = 20
    sc = rand(MersenneTwister(20),2.:10.,ncomp)
    C = [System("C[$i]",0.,sc[i],dist=Distributions.Uniform) for i = 1:ncomp];
    G1 = System("G1",Parallel(),System[C[1],C[2]])
    G2 = System("G2",Serie(),System[C[3],G1])
    G3 = System("G3",Serie(),System[C[4],C[5]])
    G4 = System("G4",Serie(),System[C[6],C[7],C[8]])
    G5 = System("G5",Parallel(),System[G3,C[9],G4])
    G6 = System("G6",Serie(),System[G5,C[10]])
    G7 = System("G7",Parallel(),System[G2,G6])
    G8 = System("G8",Parallel(),System[C[11],C[12]])
    G9 = System("G9",Parallel(),System[C[13],C[14]])
    G10 = System("G10",Parallel(),System[C[15],C[16],C[17]])
    G11 = System("G11",Parallel(),System[C[18],C[19]])
    G = System("G",Serie(),System[G7,G8,C[20],G9,G10,G11])

print("Generated a system. Number of components: ", G.aantal,"\n");
println("system: ",G.circuit)
T = Component[];
Levels,T,varT = levels_for_QMC(G);
println("Levels en Sampling voor QMC voltooid")
println(T);println(varT);
print("Levels selected! Number of levels: ",length(Levels),"\n");
lijst = G.lijst
level_levensduur(level,ω)=MultilevelSystems.lvl_levensduur(level+one(level),ω,G,Levels,lijst)
function QMC_ML_levensduur(level,ω)
level = level+one(level)
lijst = T
    # finest grid
    sampleLevensduur(lijst,ω);
    Qf = levensduur(Levels[level]);
    # compute difference when not on coarsest grid
    dQ = Qf
    if level != Level(1)
      Qc = levensduur(Levels[level-one(level)]);
      dQ -= Qc
    end
    return dQ, Qf
end
SL_levensduur(level,ω) = QMC_SL_levensduur(level+one(level),ω,G,Levels,T);
#Levels = levels;#pilot_run(G,levels)
#print("Levels removed!  Number of levels: ",length(Levels),"\n");
function distribution_QMC_levensduur(G::System,lijst)
    distributions = [MultilevelEstimators.Uniform{Float64}(0.,1.) for i in 1:length(lijst)]
    for i = 1:length(lijst)
    k = lijst[i].shape
    l = lijst[i].scale
    distributions[i] = MultilevelEstimators.Uniform(k,l)
    end
    return distributions
    end
    tol=5e-4
    elapsed = time() - start;
    display(elapsed)
    for i = 1:1
    Random.seed!(2000+i);
    distributions = distribution_QMC_levensduur(G,T);
    estimator = MultilevelEstimators.Estimator(ML(), QMC(), QMC_ML_levensduur, distributions, name="system20_uniform_MLQMC", max_index_set_param=length(Levels)-1,nb_of_warm_up_samples = 100,save_samples=true);#,robustify_bias_estimate = false);
    h2 = run(estimator, tol);
    println("finished with estimating");
    estimator = MultilevelEstimators.Estimator(ML(), MC(), level_levensduur, distributions, name="system20_uniform_MLMC", max_index_set_param=length(Levels)-1,nb_of_warm_up_samples = 100,save_samples=true);#,robustify_bias_estimate = false);
    h1 = run(estimator, tol);
    println("finished with estimating");
    end
    systeem = Dict(:systeem => G,:levels => Levels, :lijst => T,:invloed => varT,:setup_time => elapsed)
    @save "system20_uniform_system.jld2" systeem
