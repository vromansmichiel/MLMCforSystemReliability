using Plots

function levensduur_repair(G::System;extend = false)
    cutset,aantal = vindCutset(G;minimal = true, eigenschap="systeem")
    l = 1.05;
    set = [];
    for i = 1:size(cutset,1)
        defect = Bool.(zeros(1,size(cutset[i],2)));
        failtime = [];
        for j = 1:size(cutset[i],2);
            push!(failtime,cutset[i][j].failuretime);
        end
        while true
        ix = argmin(failtime);
        if defect[ix] == false
            defect[ix] = true;
            if all(defect)
                push!(set,failtime[ix])
                break
            else
                failtime[ix] += rand(Exponential(l))
            end
        elseif defect[ix] == true
            if all(defect)
                push!(set,failtime[ix])
                break
            else
                failtime[ix] += levensduur(cutset[i][ix])
            end
            defect[ix] = false
        end
        end
    end
    if extend == true
        return set
    else
        return minimum(set)
    end
end
function levensduur_repair(cutset::Array{Array{System}};extend = false)
    # functie voor het vinden van de levensduur van een cutset van een systeem
    #
    # input
    # cutset       : cutset van componenten, niet de namen of levensduur
    # output
    # levensduur :  levensduur van de cutset (extend == false)
    #               array van de levensduur van de cuts in de cutsets(extend == true)
    # functie voor het vinden van de levensduur van een cutset van een systeem
    #
    # input
    # cutset       : cutset van componenten, niet de namen of levensduur
    # output
    # levensduur :  levensduur van de cutset (extend == false)
    #               array van de levensduur van de cuts in de cutsets(extend == true)
    l = 1.05;
    set = [];
    for i = 1:size(cutset,1)
        defect = Bool.(zeros(1,size(cutset[i],2)));
        failtime = [];
        for j = 1:size(cutset[i],2);
            push!(failtime,cutset[i][j].failuretime);
        end
        while true
        ix = argmin(failtime);
        if defect[ix] == false
            defect[ix] = true;
            if all(defect)
                push!(set,failtime[ix])
                break
            else
                failtime[ix] += rand(Exponential(l))
            end
        elseif defect[ix] == true
            if all(defect)
                push!(set,failtime[ix])
                break
            else
                failtime[ix] += levensduur(cutset[i][ix])
            end
            defect[ix] = false
        end
        end
    end
    if extend == true
        return set
    else
        return minimum(set)
    end
end

function levensduur_pilot(cutset::Array{Array{System}};extend = false)
    # functie voor het vinden van de levensduur van een cutset van een systeem
    #
    # input
    # cutset       : cutset van componenten, niet de namen of levensduur
    # output
    # levensduur :  levensduur van de cutset (extend == false)
    #               array van de levensduur van de cuts in de cutsets(extend == true)
    failuretime = 0;
    l = 1.05;
    set = reshape([],0,2);
    aantal = size(cutset,1);
    defect = [];
    failtime = [];
    for i = 1:aantal
    push!(defect,Bool.(zeros(1,size(cutset[i],2))));
    push!(failtime,[]);
    set = vcat(set,[false levensduur(cutset[i:i])]);
        for j = 1:size(cutset[i],2);
            push!(failtime[i],cutset[i][j].failuretime);
        end
    end
    werkend = true
    k=1;
    while werkend
        for i = 1:aantal
        ix = argmin(failtime[i]);
        if defect[i][ix] == false
            defect[i][ix] = true;
            if all(defect[i])
                set[i,:] = [true failtime[i][ix]]
                failuretime = set[i,2];
            else
                failtime[i][ix] += rand(Exponential(l))
            end
        elseif defect[i][ix] == true
            if all(defect[i])
                werkend = false
            else
                failtime[i][ix] += levensduur(cutset[i][ix])
                set[i] = [true failtime[i][ix]]
                failuretime = set[i,2];
            end
            defect[i][ix] = false
        end

    end
    index_min = argmin(set[:,2])
    prod(set[index_min,1])==1 ? werkend = false : werkend = true
    end
    if extend == true
        return set
    else
        return minimum(set[:,2])
    end
end

function levensduur_repair(cutset::Array{Array{System}},cutset_prev::Array{Array{System}};extend = false)
    # functie voor het vinden van de levensduur van een cutset van een systeem
    #
    # input
    # cutset       : cutset van componenten, niet de namen of levensduur
    # output
    # levensduur :  levensduur van de cutset (extend == false)
    #               array van de levensduur van de cuts in de cutsets(extend == true)
    # functie voor het vinden van de levensduur van een cutset van een systeem
    #
    # input
    # cutset       : cutset van componenten, niet de namen of levensduur
    # output
    # levensduur :  levensduur van de cutset (extend == false)
    #               array van de levensduur van de cuts in de cutsets(extend == true)
    l = 1.05;
    set = [];
    for i = 1:size(cutset,1)
        defect = Bool.(zeros(1,size(cutset[i],2)));
        failtime = [];
        for j = 1:size(cutset[i],2);
            push!(failtime,cutset[i][j].failuretime);
        end
        while true
        ix = argmin(failtime);
        if defect[ix] == false
            defect[ix] = true;
            if all(defect)
                push!(set,failtime[ix])
                break
            else
                failtime[ix] += rand(Exponential(l))
            end
        elseif defect[ix] == true
            if all(defect)
                push!(set,failtime[ix])
                break
            else
                failtime[ix] += levensduur(cutset[i][ix])
            end
            defect[ix] = false
        end
        end
    end
    if extend == true
        return set
    else
        return minimum(set),minimum(set)-minimum(set[1:length(cutset_prev)])
    end
end
