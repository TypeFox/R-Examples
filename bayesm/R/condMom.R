condMom=
function(x,mu,sigi,i)
{
#
# revision history:
#    rossi modified allenby code 4/05
#
# purpose:compute moments of conditional distribution of ith element of normal given
# all others
#
# arguments:
#   x: vector of values to condition on
#   mu: mean vector of length(x)-dim MVN
#   sigi: inverse of covariance matrix
#   i: element to condition on
# 
# output:
#   list with conditional mean and variance
#
# Model: x ~MVN(mu,Sigma)
#   computes moments of x_i given x_{-1}
#
sig=1./sigi[i,i]
m=mu[i] - as.vector(x[-i]-mu[-i])%*%as.vector(sigi[-i,i])*sig
return(list(cmean=as.vector(m),cvar=sig))
}

