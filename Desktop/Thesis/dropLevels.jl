using Systems, MultilevelSystems,Distributions,MultilevelEstimators
function drop_level(G::System,levels)
    new_levels = []
    println("aantal levels: ", length(levels));
    pilot_samples = 100;
    j = 2;
    while j < length(levels)
        levels_nongeom = vcat(levels[1:j-1],levels[(j+1):end]);
        Ql = []
        Qk = []
        Tl = []
        Tk = []
        Qg = []
        for i = 1:pilot_samples
            w = [rand(Distributions.Weibull(G.lijst[i].shape,G.lijst[i].scale)) for i = 1:length(G.lijst)]
            dQ,T = @timed MultilevelSystems.level_levensduur(Level(j),w,G,levels)
            push!(Ql,dQ[1])
            push!(Tl,T)
            dQ,T = @timed MultilevelSystems.level_levensduur(Level(j+1),w,G,levels)
            push!(Qk,dQ[1])
            push!(Tk,T)
            dQ,T = MultilevelSystems.level_levensduur(Level(j),w,G,levels_nongeom)
            push!(Qg,dQ)
        end
        vl = var(Ql)
        vk = var(Qk)
        cl = mean(Tl)
        ck = mean(Tk)
        rhs = vk*(1+sqrt(vl*cl/(vk*ck)))^2
        if var(Qg) < rhs
            levels = levels_nongeom
        else
            j=j+1
        end
    end
    println("Aantal levels gehouden: ",length(levels))
    return levels
end
