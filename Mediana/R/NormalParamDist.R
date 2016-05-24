######################################################################################################################

# Function: NormalParamDist.
# Argument: c, Common critical value
#       w, Vector of hypothesis weights (1 x m)
#       corr, Correlation matrix (m x m)
# Description: Multivariate normal distribution function used in the parametric multiple testing procedure based on a multivariate normal distribution

NormalParamDist = function(c, w, corr) {

  m = dim(corr)[1]
  prob = mvtnorm::pmvnorm(lower = rep(-Inf, m), upper = c/(w * m), mean=rep(0, m), corr = corr)
  return(1 - prob[1])
}
# End of NormalParamDist