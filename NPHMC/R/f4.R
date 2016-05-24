f4 <-
function(t,accrualtime,followuptime,accrualdist,beta0,gamma0,pi0,survdist,k,lambda0){
 m(t,beta0,gamma0,pi0,survdist,k,lambda0)*f2(t,accrualtime,followuptime,accrualdist,survdist,k,lambda0)
}

