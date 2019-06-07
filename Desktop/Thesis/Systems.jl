module Systems
# This module contains the code developed for the construction of
# systems which can be used for reliability analysis of the failure time
# based on cutsets using Monte Carlo and Multilevel Monte Carlo
using Random,Distributions
export Netwerk,Parallel,Serie,System, Component,changeComponent, randNode, Node, vindCutset, geefCircuit, geefFoutenboom, numberOfComp, generateSystem, component_list, sampleLevensduur, geefComponent, levensduur,selectLevel,fullMLMC
global range
ran = [2,10]
abstract type System end
abstract type Structuur end
struct Parallel <: Structuur end
struct Serie <: Structuur end
mutable struct Component <: System
    naam::String;
    shape::Real;
    scale::Real;
    failuretime::Real;
    dist::UnionAll
    function Component(naam::String,shape::Real;dist::UnionAll=Weibull)
        l = rand(Uniform(ran[1],ran[2]));
        lvnsd = levensduur(shape,l,distr=dist);
        new(naam,shape,l,lvnsd,dist);
    end
    function Component(naam::String,shape::Real,scale::Real;dist::UnionAll=Weibull)
        lvnsd = levensduur(shape,scale,distr=dist);
        new(naam,shape,scale,lvnsd,dist);
    end
end
mutable struct Node <: System
    naam::String
    structuur::Structuur
    dochter::Array{System,1}
    function Node(naam::String,structuur::Structuur,dochter)
        for i=1:size(dochter,1)
            if typeof(dochter[i]) == Node
                if dochter[i].structuur == structuur
                    for j = 1:size(dochter[i].dochter,1)
                        insert!(dochter,i+j, dochter[i].dochter[j]);
                    end
                    deleteat!(dochter,i);
                end
            end
        end
            new(naam,structuur,dochter)
    end
end
mutable struct Netwerk
    naam::String
    nodes::Array{Component,1}
    vertices::Array{Component,2}
end
function vindPaden(N::Netwerk,s::Component,d::Component;visited= Vector{Bool}(undef, length(N.nodes)) .= false)
    paden = []
    if s == d
        error("Start is gelijk aan eind")
    else
        ix = findall(x->x==s, N.vertices[1,:])
        display(ix)
        #if length(ix) != 0
            for i = 1:length(ix)
                if N.vertices[2,ix[i]] == d
                    push!(paden,hcat(s, d))
                else
                    display(visited[findall(x->x == N.vertices[2,ix[i]],N.nodes)])
                    display(N.vertices[2,ix[i]])
                    if visited[findall(x->x == N.vertices[2,ix[i]],N.nodes)] == [false]
                        visitedpad = copy(visited)
                        visitedpad[findall(x->x == N.vertices[2,ix[i]],N.nodes)] = [true]
                        A = vindPaden(N,N.vertices[2,ix[i]],d,visited=visitedpad)
                        for i = 1:length(A)
                            push!(paden,hcat(s,A[i]))
                        end
                    else
                        continue
                    end
                end
            end
        #else
        #    return reshape(Component[],1,0)
        #end
    end
    return paden
end
function Base.getproperty(x::System, sym::Symbol)
           if sym == :levensduur
               levensduur(x)
           elseif sym == :circuit
               geefCircuit(x)
           elseif sym == :foutenboom
               geefFoutenboom(x)
           elseif sym == :cutset
                vindCutset(x)
           elseif sym == :aantal
               numberOfComp(x)
           elseif sym == :lijst
               component_list(x)
           else
               # fallback to getfield
               getfield(x, sym)
           end
       end
function System(naam::String,shape::Real,scale::Real;dist::UnionAll=Weibull)
        return Component(naam,shape,scale,dist=dist);
end
function System(naam::String,shape::Real;dist::UnionAll=Weibull)
        l = rand(Uniform(ran[1],ran[2]));
        return Component(naam,shape,l,dist=dist);
end
function System(naam::String,structuur::Structuur,dochter)
        return Node(naam,structuur,dochter)
end
function System(netwerk::Netwerk,s::Component,d::Component)
    vist= Vector{Bool}(undef, length(netwerk.nodes)) .= false
    vist[1] = true
    A = vindPaden(netwerk,s,d,visited=vist);
    subs = [System("subs$(i)",Serie(),vec(A[i])) for i = 1:length(A)];
    sys = System("S",Parallel(),subs)
    return sys
end
function levensduur(shapeparam::Real;aantal::Int=1,distr::UnionAll = Weibull)
    l = rand(Uniform(ran[1],ran[2]));
    k = shapeparam;
    if aantal == 1
        return rand(distr(k,l))
    else
    c = []
    for i = 1:aantal
        push!(c,rand(distr(k,l)))
    end
    return c
    end
    # Bepalen van de levensduur op basis van de shape parameter
end
function levensduur(shapeparam::Real,scaleparam::Real;aantal::Int=1,distr::UnionAll = Weibull)
    # Bepalen van de levensduur op basis van de shape parameter en scale parameter
    aantal == 1 ? c = rand(distr(shapeparam,scaleparam)) : c = [rand(distr(shapeparam,scaleparam)) for i = 1:aantal]
    return c
end
levensduur(component::Component;aantal::Int=1) = levensduur(component.shape,component.scale,aantal=aantal,distr=component.dist)

function levensduur(G::System;extend = false,repair = false)
    # functie voor het vinden van de levensduur van een systeem
    # input
    # G       : systeem waarvan de cutset gevonden moet worden
    # output
    # levensduur :  levensduur van het systeem (extend == false)
    #               array van de levensduur van de cuts (extend == true)
    cutset,aantal = vindCutset(G;minimal = true, eigenschap="levensduur")
    set = [];
    for i = 1:size(cutset,1)
       push!(set,maximum(cutset[i]))
    end
    if extend == true
        return set
    else
        return minimum(set)
    end
end
function levensduur(cutset::Array{Array{System}};extend = false )
    # functie voor het vinden van de levensduur van een cutset van een systeem
    #
    # input
    # cutset       : cutset van componenten, niet de namen of levensduur
    # output
    # levensduur :  levensduur van de cutset (extend == false)
    #               array van de levensduur van de cuts in de cutsets(extend == true)
    failuretime = Vector{Array{Real}}();
    for i = 1:size(cutset,1)
        cuttime = zeros(Real,1,size(cutset[i],2));
        for j = 1:size(cutset[i],2)
        cuttime[1,j] = cutset[i][j].failuretime;
        end
        push!(failuretime,cuttime);
    end
    set = [];
    for i = 1:size(failuretime,1)
       push!(set,maximum(failuretime[i]))
    end
    if extend == true
        return set
    else
        return minimum(set)
    end
end
function levensduur(cutset::Array{Array{System}},cutset_prev::Array{Array{System}};extend = false)
    if !issubset(cutset_prev,cutset[1:size(cutset_prev,1)])
    error("cutset2 is niet het voorgaande level van cutset1")
    else
    # functie voor het vinden van de levensduur van een cutset van een systeem
    #
    # input
    # cutset       : cutset van componenten, niet de namen of levensduur
    # output
    # levensduur :  levensduur van de cutset (extend == false)
    #               array van de levensduur van de cuts in de cutsets(extend == true)
    failuretime = Vector{Array{Real}}();
    for i = 1:size(cutset,1)
        cuttime = zeros(Real,1,size(cutset[i],2));
        for j = 1:size(cutset[i],2)
        cuttime[1,j] = cutset[i][j].failuretime;
        end
        push!(failuretime,cuttime);
    end
    set = [];
    for i = 1:size(failuretime,1)
       push!(set,maximum(failuretime[i]))
    end
        return minimum(set)-minimum(set[1:size(cutset_prev,1)])
    end
end
function vindCutset(G::System;minimal = false,eigenschap = "naam")
    # functie voor het vinden van de cutset
    # input
    # G       : systeem waarvan de cutset gevonden moet worden
    # minimal : variable die aangeeft of de cutset minimaal moet zijn of
    #               dat de volledige set als return is gevraagd.
    # eigenschap: welke eigenschap in de cutset moet staan
    # output
    # cutset : (minimale) cutset
    # aantal : aantal cutsets in het circuit
    if eigenschap == "naam"
    cutset = Vector{Array{String}}()
    elseif eigenschap == "levensduur"
    cutset = Vector{Array{Real}}()
    elseif eigenschap == "systeem"
        cutset = Vector{Array{System}}()
    else
        error("verkeerde eigenschap");
    end
    if typeof(G) == Component
        if eigenschap == "naam"
            set = [G.naam]
        elseif eigenschap == "levensduur"
            set = [G.failuretime]
        elseif eigenschap == "systeem"
            set = [G]
        end
        cutset = push!(cutset,set)
    elseif G.structuur == Serie()
        for i = 1 : size(G.dochter,1)
            set,Z = vindCutset(G.dochter[i];eigenschap = eigenschap);
            rows = size(set,1)
            cols = size(set,2)
            for i = 1:size(set,1)
            push!(cutset,set[i])
            end
        end
    elseif G.structuur == Parallel()
        for i = 1:size(G.dochter,1)
            cut,A = vindCutset(G.dochter[i];eigenschap = eigenschap);
            row = size(cut,1)
            col = size(cut,2)
            if i == 1
                cutset = cut
            else
                cRows = size(cutset,1)
                cCols = size(cutset,2)
                cutset = repeat(cutset,outer =(row,1))
                cut = repeat(cut,inner =(cRows,1))
                for j = 1:size(cutset,1)
                    cutset[j] = hcat(cutset[j] ,cut[j])
                end
            end
        end
    end
    cutset = vec(cutset);
    typeof(cutset) == String ? aantal = 1 : aantal = size(cutset,1);
    if minimal == true
        cutset = minimalSet(cutset);
        aantal = size(cutset,1);
    end
    return cutset,aantal
end
function vindPathset(G::System;minimal = false,eigenschap = "naam")
    # functie voor het vinden van de cutset
    # input
    # G       : systeem waarvan de pathset gevonden moet worden
    # minimal : variable die aangeeft of de pathset minimaal moet zijn of
    #               dat de volledige set als return is gevraagd.
    # eigenschap: welke eigenschap in de cutset moet staan
    # output
    # pathset : (minimale) pathset
    # aantal : aantal pathsets in het circuit
    if eigenschap == "naam"
    pathset = Vector{Array{String}}()
    elseif eigenschap == "levensduur"
    pathset = Vector{Array{Real}}()
    elseif eigenschap == "systeem"
        pathset = Vector{Array{System}}()
    else
        error("verkeerde eigenschap");
    end
    if typeof(G) == Component
        if eigenschap == "naam"
            set = [G.naam]
        elseif eigenschap == "levensduur"
            set = [G.failuretime]
        elseif eigenschap == "systeem"
            set = [G]
        end
        push!(pathset,set)
    elseif G.structuur == Parallel()
        for i = 1 : size(G.dochter,1)
            set,Z = vindPathset(G.dochter[i];eigenschap = eigenschap);
            rows = size(set,1)
            cols = size(set,2)
            for i = 1:size(set,1)
            push!(pathset,set[i])
            end
        end
    elseif G.structuur == Serie()
        for i = 1:size(G.dochter,1)
            cut,A = vindCutset(G.dochter[i];eigenschap = eigenschap);
            row = size(cut,1)
            col = size(cut,2)
            if i == 1
                pathset = cut
            else
                cRows = size(cutset,1)
                cCols = size(cutset,2)
                cutset = repeat(cutset,outer =(row,1))
                cut = repeat(cut,inner =(cRows,1))
                for j = 1:size(cutset,1)
                    pathset[j] = hcat(pathset[j] ,cut[j])
                end
            end
        end
    end
    pathset = vec(pathset);
    typeof(pathset) == String ? aantal = 1 : aantal = size(pathset,1);
    if minimal == true
        pathset = minimalSet(pathset);
        aantal = size(pathset,1);
    end
    return pathset,aantal
end
function minimalSet(s;arg = length)
    # functie die een bepaalde set vereenvoudigd, de redundante sets worden verwijderd.
    # redundante sets zijn sets waarvan een subset een cutset is van het circuit.
    # input
    #   s   : array van cutsets
    # output
    #   s   : array zonder redundante sets
    for i = 1:length(s)
    s[i] = reshape(fastuniq(s[i][1:end]),1,:)
    end
    arg == length ? sort!(s,by=arg) : sort!(s)
    j = 1
    while j <= size(s,1)
            max = size(s,1)
            i = j+1
        while i <= size(s,1)
            issubset(s[j],s[i]) ? deleteat!(s,i) : i = i+1
        end
        j=j+1
    end
        return s
end
function geefCircuit(systeem::System)
    # deze functie vertaalt een circuit naar een string waarbij
    # Componenten in serie worden verbonden met een & teken en
    # Componenten in parallel met een | teken.
    # input
    # systeem   : te bekijken circuit
    # output
    # s         : formule van het systeem dat het circuit weergeeft
    if typeof(systeem) == Component
        return systeem.naam
    elseif systeem.structuur == Serie()
        s = "("*geefCircuit(systeem.dochter[1])
        for i = 2:size(systeem.dochter,1)
        s = s*"&"*geefCircuit(systeem.dochter[i])
        end
    elseif systeem.structuur == Parallel()
        s = "("*geefCircuit(systeem.dochter[1])
        for i = 2:size(systeem.dochter,1)
        s = s*"|"*geefCircuit(systeem.dochter[i])
        end
    end
    s = s*")"
    return s
end
function geefFoutenboom(systeem::System)
    # deze functie vertaalt een circuit naar een foutenboom
    # een & teken is een AND poort en een | teken is een OF poort.
    # input
    # systeem   : te bekijken circuit
    # output
    # s         : formule van de foutenboom van het circuit
    if typeof(systeem) == Component
        return systeem.naam
    elseif systeem.structuur == Serie()
        s = "("*geefFoutenboom(systeem.dochter[1])
        for i = 2:size(systeem.dochter,1)
        s = s*"|"*geefFoutenboom(systeem.dochter[i])
        end
    elseif systeem.structuur == Parallel()
        s = "("*geefFoutenboom(systeem.dochter[1])
        for i = 2:size(systeem.dochter,1)
        s = s*"&"*geefFoutenboom(systeem.dochter[i])
        end
    end
    s = s*")"
    return s
end
function randNode(naam::String,kansS::Real,kansP::Real,shape::Real;random = MersenneTwister())
    # deze functie creeert een random Node
    # Input
    # kansP : kans op parallelschakeling
    # kansS : kans op serieschakeling
    # shape : de shape parameter van de componenten
    # Output
    # systeem: struct Node
    r = rand(random)
    if r < kansS
        s1 = System("$(naam)0",shape)
        s2 = System("$(naam)1",shape)
        systeem = System(naam,Serie(),System[s1,s2])
    elseif r < kansS+kansP && r > kansS
        p1 = System("$(naam)0",shape)
        p2 = System("$(naam)1",shape)
        systeem = System(naam,Parallel(),System[p1,p2])
    else
        systeem = System("$(naam)b",shape)
    end
    return systeem
end
function numberOfComp(systeem::System)
    # functie die bepaald hieveel Componenten er in het systeem zitten
    aantal = 0
    if typeof(systeem) == Component
        aantal = 1
    else
        for i = 1:size(systeem.dochter,1)
            aantal += numberOfComp(systeem.dochter[i]);
        end
    end
    return aantal
end
function generateSystem(naam::String,aantalComp::Int,kS::Real,kP::Real,shapeparam::Real;random = MersenneTwister())
    #functie die een random systeem genereert met aantalComp Componenten.
    #input
    # naam : naam van het systeem
    # aantalComp : hoeveel Componenten er in het systeem aanwezig moeten zijn
    # kansS, kansP: kans op een serie, resp. parallelschakeling
    # shapeparam : shapeparameter van de componenten
    #opmerking: als kansS+kansP < 1 dan is er een kans dat er een bridge gemaakt wordt tussen 2 componenten
    kansS = kS/(kS+kP);
    kansP = kP/(kS+kP)
  if aantalComp == 0
        error("leeg systeem")
  elseif aantalComp == 1
        systeem = System(naam,shapeparam);
  else
        systeem = randNode(naam,kansS,kansP,shapeparam;random=random);
        for i = 3:aantalComp
          r = rand(1:numberOfComp(systeem));
          nm = geefComponent(systeem,r,eigenschap="naam");
          if rand(Uniform(0,1)) > kS+kP
              systeem.lijst[rand(1:systeem.aantal)] = systeem.lijst[rand(1:systeem.aantal)]
          end
          A = randNode(nm ,kansS,kansP,shapeparam;random=random);
          comp = changeComponent(systeem,A,r);
        end
  end
  return systeem
end
function changeComponent(G::System,S::System,naam::String)
  for i = 1:size(G.dochter,1)
    if(typeof(G.dochter[i])==Component)
        if(G.dochter[i].naam == naam)
            if G.structuur == S.structuur
                deleteat!(G.dochter,i)
                for i = 1:size(S.dochter,1)
                push!(G.dochter,S.dochter[i]);
                end
            else
                G.dochter[i] = S;
            end
       end
    elseif (typeof(G.dochter[i])==Node)
        changeComponent(G.dochter[i],S,naam);
    end
  end
end
function changeComponent(G::System,S::System,nummer::Int)
  if nummer > numberOfComp(G)
    error("nummer is groter dan het systeem")
  end
  minus = 0;
  for i = 1:size(G.dochter,1)
    if nummer <= minus+numberOfComp(G.dochter[i]);
      if numberOfComp(G.dochter[i]) == 1
         if G.structuur == S.structuur
              deleteat!(G.dochter,i)
              for i = 1:size(S.dochter,1)
                  push!(G.dochter,S.dochter[i]);
              end
          else
              G.dochter[i] = S;
          end
      else
          changeComponent(G.dochter[i],S,nummer-minus);
      end
      break
    end
  minus+=numberOfComp(G.dochter[i]);
  end
end
function sampleLevensduur(cutset)
    for i = 1:size(cutset,1)
        for j = 1:size(cutset[i],2)
        cutset[i][j].failuretime = levensduur(cutset[i][j]);
        end
    end
end
function sampleLevensduur(G::System)
    if typeof(G) == Component
        G.failuretime = levensduur(G.shape,G.scale,distr=G.dist);
    else
        for i = 1:size(G.dochter,1)
            sampleLevensduur(G.dochter[i]);
        end
    end
end
function sampleLevensduur(component_lijst::Vector{Component},levensduur_lijst::Vector{Float64})
    size(component_lijst) != size(levensduur_lijst) ? error("Lengte van de vectoren komen niet overeen") : nothing;
    for i = 1:length(component_lijst)
        component_lijst[i].failuretime = levensduur_lijst[i];
    end
end
function geefComponent(G::System,nummer::Int;eigenschap = "systeem")
  C = 0;
  nummer > numberOfComp(G) ? error("nummer is groter dan het systeem") : nothing
  minus = 0;
  for i = 1:size(G.dochter,1)
    if nummer <= minus+numberOfComp(G.dochter[i]) ;
      if numberOfComp(G.dochter[i]) == 1
          if eigenschap == "levensduur"
              return G.dochter[i].failuretime;
          elseif eigenschap == "naam"
              return G.dochter[i].naam;
          else
              return G.dochter[i];
          end
      else
            return geefComponent(G.dochter[i],nummer-minus,eigenschap = eigenschap);
      end
    end
  minus+=numberOfComp(G.dochter[i]);
  end
end
function component_list(systeem::System)
    s = Vector{Component}()
    if typeof(systeem) == Component
        s = systeem
    else
        for i = 1:size(systeem.dochter,1)
            s = vcat(s,component_list(systeem.dochter[i]))
        end
        s = fastuniq(s);
    end
    return s;
end
import Base.isless
function isless(component1::Component,component2::Component)
    return isless(component1.naam,component2.naam)
end
#Dan Getz : https://stackoverflow.com/questions/36517553/remove-consecutive-duplicates-in-julia
function fastuniq(v)
  sort!(v)
  v1 = Vector{eltype(v)}()
  if length(v)>0
    laste = v[1]
    push!(v1,laste)
    for e in v
      if e != laste
        laste = e
        push!(v1,laste)
      end
    end
  end
  return v1
end
end
