module MultilevelSystems
using Random,Systems,MultilevelEstimators,Distributions,Printf
export MC_levensduur, multilevel_levensduur, sample_levensduur, level_levensduur, distribution_levensduur, selectLevel
export QMC_ML_levensduur,levels_for_QMC,QMC_SL_levensduur
function MC_levensduur(G::System,N::Int)
    set,a = vindCutset(G,minimal=true,eigenschap="systeem");
    l = [];
    for i = 1:N
        sampleLevensduur(G);
        push!(l,levensduur(set));
    end
    return mean(l),var(l);
    end
function MC_levensduur(set,N::Int)
    l = [];
    for i = 1:N
        sampleLevensduur(G);
        push!(l,levensduur(set));
    end
    return mean(l),var(l);
    end
function MC_levensduur(G::System,v::Real;error::Real = 5e-2)
    set,a = vindCutset(G,minimal=true,eigenschap="systeem");
    l = [];
    N = ceil(4^2*v*(error/2)^(-2));
    for i = 1:N
        sampleLevensduur(G);
        push!(l,levensduur(set));
    end
    return mean(l),var(l);
    end
function MC_levensduur(G::System;error::Real = 5e-2)
    set,a = vindCutset(G,minimal=true,eigenschap="systeem");
    l = [];
    for i = 1:100
        sampleLevensduur(G);
        push!(l,levensduur(set));
    end
    v = var(l);
    return levensduur(G,set,v,error);
    end
function MC_levensduur(G::System,set,N::Int )
    l = [];
    for i = 1:N
        sampleLevensduur(G);
        push!(l,levensduur(set));
    end
    return mean(l),var(l);
    end
function MC_levensduur(G::System,set;error::Real = 5e-2)
    l = [];
    for i = 1:100
        sampleLevensduur(G);
        push!(l,levensduur(set));
    end
    v = var(l);
    return levensduur(G,set,v,error);
    end
function MC_levensduur(G::System,set,v::Real;error::Real= 5e-2)
    l = [];
    N = ceil(4^2*v*(error/2)^(-2));
    for i = 1:N
        sampleLevensduur(G);
        push!(l,levensduur(set));
    end
    return mean(l),var(l);
    end
function evaluateMC(G::System )
    gem = []
    variantie = []
    cost = []
    e = []
    l = []
    n = []
    a = []
    m = []
    for j = 1:15
    N = 2^j;
    MCmean = [];
    MCvar = [];
    MCcost = [];
    (mn,vr),t = @timed MC_levensduur(G,N);
    for i = 1:10
        (mn,vr),t = @timed MC_levensduur(G,N);
        push!(MCmean,mn);push!(MCvar,vr);push!(MCcost,t);
    end
    gemiddelde = mean(MCmean);
    error = (MCmean.-gemiddelde).^2;
    push!(e,mean(error));
    push!(gem,mean(MCmean));
    push!(variantie,var(MCmean));
    push!(l,mean(MCvar));
    push!(cost,mean(MCcost));
    push!(n,N);
    push!(a,1/(N));
    end
    f = plot(log2.(n),log2.(a),label="reference");plot!(f,log2.(n),log2.(variantie),label="MC variantie",xlabel ="log2(N)",ylabel="log2(variantie)");
        display(f)
    g = plot(log2.(n),log2.(n),label="reference");plot!(g,log2.(n),log2.(cost),label="MC cost",xlabel ="log2(N)",ylabel="log2(cost)");
    display(g)
    return hcat(gem,variantie,l,n,e,cost)
    end
selectLevel(G::System) = selecteerLevel(G)
function selecteerLevel(G::System)
    set,a = vindCutset(G,minimal=true,eigenschap="systeem");
    L = Int(max(floor(log2(a)-2),3));
    nr = Int(ceil(a/(2^L)));
    N = 100; # pilot MC
    pilot = [levensduur(set,extend=true) for i = 1:N];
    pilot = reshape(pilot,length(set),N);
    l = mean(pilot,dims=2); # level 0
    lst = sort(l,dims=1)[nr];
    ix = findall(l.<=lst);
    level_prev = set[ix];
    ix = findall(l.>lst);
    ctrial = set[ix];
    T_prev,Tvar = MC_levensduur(G,level_prev,100);
    levels = [level_prev];
    for l = 1:L # level 1 tot L
      l==L ? nr = size(set,1)-size(level_prev,1) : nr = Int(ceil(a/2^(L-l)))-Int(ceil(a/2^(L-(l-1))))
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
    return levels;
    end
function multilevel_levensduur(G::System,levels;err=5e-2)
  L_init = 3;
  N = 100*ones(1,L_init);
  N_prev = zeros(1,L_init);
  V = reshape([],1,0);
  Y = reshape([],1,0);
  y = [[]]
  for i = 2:L_init
    push!(y,[]);
  end
  for maxL = L_init:length(levels)
      println("samplen tot level ", maxL-1)
    while any(N .> N_prev.*1.01)
        Y = reshape([],1,0);
        V = reshape([],1,0);
        N_new = reshape([],1,0);
        println("Level:    N:      E:           V: ")
        for j = 1:maxL
          for i = N_prev[j]+1:N[j]
              w = [rand(Distributions.Weibull(G.lijst[i].shape,G.lijst[i].scale)) for i = 1:length(G.lijst)]
              dQ, T = level_levensduur(Level(j),w,G,levels)
              push!(y[j],dQ);
          end
            V = hcat(V,var(y[j]));
            Y = hcat(Y,mean(y[j]));
            @printf "%4.0f %10.0f    %2.5e     %2.5e \n" j-1 N[j] Y[j] V[j]
        end
        som = 0;
        for i = 1:maxL
            som += sqrt(V[i]*2^(i-1));
        end
        for i = 1:maxL
            N_new = hcat(N_new,ceil(2*err^(-2)*sqrt(V[i]*2.0^(-(i-1)))*som));
        end
        N_prev = N;
        N = maximum(vcat(N_prev,N_new),dims=1);
        println("aantal bijkomende samples:")
        for i = 1:length(N)
            println("Level $(i-1): ", N[i]-N_prev[i])
        end
    end
    abs(Y[end]) < 1/sqrt(2)*err ? println("stop hier maar met samplen") : println("sample op een hoger level")
      V = hcat(V,V[end]/2);
      som=0;
      N_prev = [N 0];
      N = reshape([],1,0);
    #  println("Variantie per level : ", V)
    for i = 1:size(V,2)
        som += sqrt(V[i]*2^(i-1));
    end
    for i = 1:size(V,2)
          N = hcat(N,max(N_prev[i],ceil(2*err^(-2)*sqrt(V[i]*2.0^(-(i-1)))*som)));
    end
    #println("aantal N  per level : ", N)
    println("aantal bijkomende samples:")
    for i = 1:length(N)
        println("Level $(i-1): ", N[i]-N_prev[i])
    end
    push!(y,[])
  end
  return sum(Y);
  end
multilevel_levensduur(G::System;err=5e-2)=multilevel_levensduur(G,selectLevel(G),err=err)
function lvl_levensduur(level::Level,ω::Vector{<:Real},G::System,Levels,lijst::Vector{<:System})
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
level_levensduur(level::Level,ω::Vector{<:Real},G::System,Levels) = lvl_levensduur(level,ω,G,Levels,component_list(G))
QMC_ML_levensduur(level::Level,ω::Vector{<:Real},G::System,Levels,lijst) = lvl_levensduur(level,ω,G,Levels,lijst)
function distribution_levensduur(G::System)
    lijst = component_list(G);
    distributions = [MultilevelEstimators.Weibull() for i in 1:length(lijst)]
    for i = 1:length(lijst)
    k = lijst[i].shape
    l = lijst[i].scale
    distributions[i] = MultilevelEstimators.Weibull(k,l)
    end
    return distributions
    end
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
     ω = [rand(Distributions.Weibull(G.lijst[i].shape,G.lijst[i].scale)) for i = 1:length(G.lijst)]
     A = hcat(A,ω);
     B = hcat(B,[rand(Distributions.Weibull(G.lijst[i].shape,G.lijst[i].scale)) for i = 1:length(G.lijst)]);
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
     ω = [rand(Distributions.Weibull(G.lijst[i].shape,G.lijst[i].scale)) for i = 1:length(G.lijst)]
     A = hcat(A,ω);
     B = hcat(B,[rand(Distributions.Weibull(G.lijst[i].shape,G.lijst[i].scale)) for i = 1:length(G.lijst)]);
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
end
