using Random
function creeerToyExample(shape::Real; ncomp = 7)
    Random.seed!(2019)
    if ncomp == 20
        scale = rand(MersenneTwister(20),2.:10.,ncomp)
        C = [System("C[$i]",shape,scale[i]) for i = 1:ncomp];
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
    elseif ncomp == 50
        scale = rand(MersenneTwister(50),2.0:3.0,ncomp)
        C = [System("C[$i]",shape,scale[i]) for i = 1:ncomp];
        G1 = System("G1",Serie(),System[C[1],C[2],C[3]])
        G2 = System("G2",Parallel(),System[G1,C[4],C[5],C[6],C[7],C[8]])
        G3 = System("G3",Parallel(),System[C[9],C[10],C[11],C[12]])
        G4 = System("G4",Parallel(),System[C[13],C[14]])
        G5 = System("G5",Serie(),System[G4,C[15]])
        G6 = System("G6",Parallel(),System[G5,C[16]])
        G7 = System("G7",Parallel(),System[C[17],C[18]])
        G8 = System("G8",Serie(),System[C[19],G7])
        G9 = System("G9",Parallel(),System[C[20],G8])
        G10 = System("G10",Serie(),System[C[21],C[22],C[23]])
        G11 = System("G11",Serie(),System[C[24],C[25],C[26]])
        G12 = System("G12",Serie(),System[C[27],C[28]])
        G13 = System("G13",Parallel(),System[C[29],G10,G11,C[30],G12,C[31]])
        G14 = System("G14",Serie(),System[C[32],C[33]])
        G15 = System("G15",Serie(),System[C[34],C[35]])
        G16 = System("G16",Serie(),System[C[36],C[37]])
        G17 = System("G17",Parallel(),System[C[38],G14,G15,G16])
        G18 = System("G18",Parallel(),System[C[39],C[40],C[41]])
        G19 = System("G19",Serie(),System[C[42],G18])
        G20 = System("G20",Serie(),System[C[43],C[44]])
        G21 = System("G21",Serie(),System[G13,G17])
        G22 = System("G22",Serie(),System[C[45],C[46]])
        G23 = System("G23",Serie(),System[G22,C[47]])
        G24 = System("G24",Serie(),System[G2,G3,G6,C[48],G9])
        G = System("G",Parallel(),System[G24,G21,C[49],G19,G20,G23,C[50]])
    elseif ncomp == 30
        scale = rand(MersenneTwister(50),2.:3.,ncomp)
        C = [System("C[$i]",shape,scale[i]) for i = 1:ncomp];
        G16 = System("G16",Parallel(),System[C[29],C[30]]);
        G15 = System("G15",Parallel(),System[C[27],C[28]]);
        G14 = System("G14",Serie(),System[G15,G16,C[26]]);
        G13 = System("G13",Serie(),System[C[24], C[25]]);
        G12 = System("G12",Parallel(),System[G13,G14]);
        G11 = System("G11",Parallel(),System[C[21],C[22],C[23]]);
        G10 = System("G10",Serie(),System[C[20],G11,G12]);
        G9 = System("G9",Serie(),System[C[18],C[19]]);
        G8 = System("G8",Parallel(),System[C[14],C[15]]);
        G7 = System("G7",Serie(),System[C[13],G8]);
        G6 = System("G6",Serie(),System[C[9],C[10],C[11],C[12]]);
        G5 = System("G5",Serie(),System[C[6],C[7],C[8]]);
        G4 = System("G4",Parallel(),System[C[3],C[4]]);
        G3 = System("G3",Parallel(),System[C[16],C[17],G9]);
        G2 = System("G2",Parallel(),System[C[5],G5,G6,G7]);
        G1 = System("G1",Serie(),System[G2,G3,G4]);
        G = System("G",Parallel(),[G1,G10,C[1],C[2]]);
    else
        C1 = System("C1",shape);
        C2 = System("C2",shape);
        C3 = System("C3",shape);
        C4 = System("C4",shape);
        C5 = System("C5",shape);
        C6 = System("C6",shape);
        C7 = System("C7",shape);
        C8 = System("C8",shape);
        G6 = System("G6",Serie(),System[C6,C8]);
        G5 = System("G5",Serie(),System[C6,C7]);
        G4 = System("G4",Serie(),System[C4,C5]);
        G3 = System("G3",Serie(),System[C3,C6]);
        G2 = System("G2",Parallel(),System[G4,G5]);
        G1 = System("G1",Serie(),System[G2,G3]);
        subG = System("subG",Parallel(),System[G3,G4]);
        subG2 = System("subG2",Serie(),System[G2,C1]);
        subG3 = System("subG3",Serie(),System[C3,G2,C1]);
        G = System("G",Serie(),[C1,G1,C2]);
    end
    return G
end
