using Systems,Distributions
function brandalarm(k)
      DC = Component("DC",k,10.)
      SR = Component("SR",k,10.)
      PS = Component("PS",k,10.)
      OP = Component("OP",k,15.)
      MS = Component("MS",k,8.)
      M = Node("MA",Serie(),[OP,MS,PS])
      VU = Component("VU",k,8.)
      SD1 = Component("SD1",k,7.)
      SD2 = Component("SD2",k,7.)
      SD3 = Component("SD3",k,7.)
      COMB1 = Node("COMB1",Parallel(),[SD1,SD2])
      COMB2 = Node("COMB2",Parallel(),[SD1,SD3])
      COMB3 = Node("COMB3",Parallel(),[SD2,SD3])
      RM = Node("RM",Serie(),[COMB1,COMB2,COMB3,VU])
      FP1 = Component("FP1",k,7.)
      FP2 = Component("FP2",k,7.)
      FP3 = Component("FP3",k,7.)
      FP4 = Component("FP4",k,7.)
      FP = Node("FP",Parallel(),[FP1,FP2,FP3,FP4])
      WD = Node("WD",Serie(),[FP,PS])
      DET = Node("DET",Parallel(),[M,RM,WD])
      Systeem = Node("Systeem",Serie(),[DC,DET,SR])
      return Systeem
end
function hydrocentrale()
      CG = Component("CG",1.8,3.)
      TB1 = Component("TB1",2.3,1.8)
      TB2 = Component("TB2",2.3,1.8)
      T1 = Component("T1",3.,4.)
      T2 = Component("T2",3.,4.)
      G1 = Component("G1",2.6,2.1)
      G2 = Component("G2",2.6,2.1)
      CB1 = Component("CB1",1,4.)
      CB2 = Component("CB2",1,4.)
      C1 = Node("C1",Serie(),[TB1,T1,G1,CB1])
      C2 = Node("C1",Serie(),[TB2,T2,G2,CB2])
      C = Node("C",Parallel(),[C1,C2])
      TX1 = Component("TX1",1.,3.,dist=Gamma)
      TX2 = Component("TX2",1.,3.,dist=Gamma)
      TX = Node("C",Parallel(),[TX1,TX2])
      CB3 = Component("CB3",1,3.85)
      Systeem = Node("Systeem",Serie(),[CG,C,CB3,TX])
      hazard1 = Hazard(Distributions.Weibull(2.,1))
      hazard2 = Hazard(Distributions.Weibull(3.,1.44))
      dependency = [[CG  hazard1];[T2 hazard1];[hazard2 CB3];[TX2 hazard2];[CB3 hazard2];[CB1 CB3];[TX1 CB1]];
      return Systeem,dependency
end
function boeing747()
      A = reshape([Component("A[$i]",1.,1.) for i = 1:14],1,:);
      B = reshape([Component("B[$i]",1.,1.) for i = 1:16],1,:);
      C = reshape([Component("C[$i]",1.,1.) for i = 1:6],1,:);
      D = reshape([Component("D[$i]",1.,1.) for i = 1:20],1,:);
      E = reshape([Component("E[$i]",1.,1.) for i = 1:14],1,:);
      F = reshape([Component("F[$i]",1.,1.) for i = 1:10],1,:);
      S = Component("S",1.,1.)
      T = Component("T",1.,1.)
      vertices = [S S S S S S S S S S ; A[1] A[2] A[4] A[6] A[7] A[8] A[9] A[11] A[13] A[14] ]
      vertices = hcat(vertices,
      [A[1] A[2] A[5] A[1] A[2] A[4] A[5] A[6] A[7] A[8] A[9]  A[10] A[11] A[13] A[14] A[5] A[5] A[10] A[10] A[9] A[10] A[11] A[13] A[14];
       B[1] B[3] B[6] A[3] A[3] A[3] A[3] B[7] B[8] B[9] B[10] B[11] B[13] B[14] B[15] A[6] A[8] A[8]  A[9]  A[7] A[12] A[12] A[12] A[12]])
      vertices = hcat(vertices,
      [B[4] B[6] B[7] B[8] B[9] B[10] B[11] B[13];
       C[1] C[1] C[2] C[3] C[4] C[5]  C[6]  C[6]])
      vertices = hcat(vertices,
      [C[1] C[2] C[3] C[4]  C[4]  C[5]  C[6] B[3] B[14];
       D[2] D[8] D[9] D[10] D[11] D[13] D[4] D[3] D[17]])
      vertices = hcat(vertices, [F[5] E[5] E[8] F[10];T T T T])
      vertices = hcat(vertices,
      [F[1] F[1] F[4] F[2] F[3] F[6] F[8] F[8] F[7]  F[9];
       F[3] F[2] F[2] F[5] F[5] F[7] F[7] F[9] F[10] F[10]])
      vertices = hcat(vertices,
      [D[1] D[3] D[7] D[8] D[9] D[12] D[13] D[14] D[16] D[17] D[19];
       E[1] E[2] E[4] E[5] E[6] E[7]  E[8]  E[9]  E[10] E[11] E[12]])
      vertices = hcat(vertices,
      [E[1] E[2] E[3] E[4] E[4] E[9] E[9] E[10] E[11] E[12] E[1] E[2] E[3]  E[4]  E[9]  E[10] E[11] E[12];
       F[1] F[1] F[1] F[1] F[4] F[6] F[8] F[8]  F[8] F[8] E[14] E[14] E[14] E[14] E[13] E[13] E[13] E[13]])
      vertices = hcat(vertices,
      [E[4] E[4] E[5] E[6] E[7] E[8] E[5] E[7] E[6] E[8] E[9] E[9];
       E[5] E[7] E[6] E[8] E[9] E[9] E[4] E[4] E[5] E[6] E[7] E[8]])
      vertices = hcat(vertices,
      [ B[1]  B[1] B[1] B[1] B[3] B[3] B[3] B[4] B[4] B[4] B[6] B[6] B[6] B[11] B[11] B[11] B[11] B[13] B[13] B[13] B[14] B[14] B[14] B[15] B[15];
        B[11] B[3] B[2] B[5] B[2] B[5] B[4] B[2] B[5] B[6] B[2] B[5] B[7] B[10] B[12] B[16] B[13] B[12] B[16] B[14] B[12] B[16] B[15] B[12] B[16]
      ])
      vertices = hcat(vertices,
      [ B[2] B[2] B[5] B[5] B[5] B[7] B[8]  B[9]  B[9]  B[10] B[10];
        B[7] B[9] B[7] B[8] B[9] B[8] B[10] B[12] B[16] B[16] B[12]])
        vertices = hcat(vertices,
       [D[1] D[1] D[3] D[3] D[5] D[6] D[7] D[7] D[7] D[14] D[14] D[14] D[16] D[16] D[17] D[17] D[19] D[19];
        D[2] D[6] D[2] D[6] D[2] D[6] D[2] D[6] D[4] D[15] D[18] D[20] D[15] D[20] D[15] D[20] D[15] D[20]])
        vertices = hcat(vertices,
       [D[4] D[4]  D[6] D[6]  D[8] D[9]  D[10] D[10] D[10] D[10] D[12] D[12] D[13] D[13];
        D[8] D[12] D[8] D[12] D[9] D[13] D[11] D[12] D[18] D[15] D[13] D[11] D[15] D[18]])
        vist = Vector{Bool}(undef, 82) .= false
      N = Netwerk("N",vcat(S,vec(A),vec(B),vec(C),vec(D),vec(E),vec(F),T),vertices)
      return System(N,S,T)
end
