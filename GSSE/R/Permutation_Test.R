
#' @export

Permutation_Test = function(Grid, F_carr, F_non, OY, ODelta, Op0G, nperm = 1000){

# -------------------------- #
# -    Permutation Test    - #
# -------------------------- #

# nperm = 1000;
n = length(OY);

Per.value = matrix(0, nrow=nperm, ncol=1);

ind = order(OY);
Y = OY[ind]; Delta = ODelta[ind]; p0G = Op0G[ind];

#SieveNPMLE_result = Sieve_NPMLE_Switch(Y=Y, p0G=p0G, Delta=Delta, px=Grid, Grid=Grid, Knot=4, degree =3);
#Sieve_F1.Grid = 1 - exp( - SieveNPMLE_result$Lamb1_Grid );
#Sieve_F2.Grid = 1 - exp( - SieveNPMLE_result$Lamb2_Grid );

Test0.value = test_stat ( F1 = F_carr, F2 = F_non );

for(b in 1:nperm){
  PermuIndex = sample(seq(1,n), replace=F);
  OPer.Y = OY[PermuIndex]; P.ind = order(OPer.Y);
  OPer.Delta = ODelta[PermuIndex];
  Per.Y = OPer.Y[P.ind];
  Per.Delta = OPer.Delta[P.ind];
  Per.p0G = Op0G[P.ind];

  Per.SieveNPMLE_result = Sieve_NPMLE_Switch (Y=Per.Y, p0G=Per.p0G, Delta=Per.Delta, px=Grid, Grid=Grid, Knot=4, degree =3);
  Per.F1.Grid = 1 - exp( - Per.SieveNPMLE_result$Lamb1 );
  Per.F2.Grid = 1 - exp( - Per.SieveNPMLE_result$Lamb2 );

  Per.value[b,1] = test_stat ( F1 = Per.F1.Grid, F2 = Per.F2.Grid );
}

pvalues = sum(Per.value[,1] >= Test0.value)/nperm;

return(list(Test_Stat = Test0.value, Pvalues = pvalues, Permu.value=Per.value))

}

