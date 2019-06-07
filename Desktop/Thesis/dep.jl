mutable struct Hazard
    distr
    failuretime::Real;
    function Hazard(distr)
        new(distr,rand(distr));
    end
end


function dependence(G,dependency)
       for i = 1:size(dependency,1)
           dependency[i,1].failuretime = minimum(x->x.failuretime,dependency[i,:])
           #end
       end
end
