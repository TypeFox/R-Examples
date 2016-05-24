######################################################################
# Wrapper function for calling the Fortran based calculations
# for one lambda value
######################################################################

liuSamSize <- function(alpha, d, lambda) {
  #Call the Code bei Wei Liu
  res <- .Fortran("liusamsize",alpha=as.double(alpha), d=as.double(d), lambda=as.double(lambda), nstar=as.double(0), n0=as.double(0), cp=as.double(0), package="binomSamSiz",PACKAGE="binomSamSize")
  
  #Manipulate output where necessary, i.e. round n0
  res$n0 <- ceiling(res$n0)

  #Done and return list with the three important entries
  return(unlist(res[c("nstar","cp")]))
}

######################################################################
# Calculate sample size for a binomial parameter based on a
# confidence interval width specification.
#
# The objective is to find the minimum sample size nn so that
# the minimum coverage probability of the confidence interval
# for the binomial parameters based on [c(x) +- d] based on equation (3.1)
# in Liu & Bailey (2002) and expanded to the full length
# of 2d, is no less than 1-alpha.
#
# Parameters:
#   alpha - an (1-alpha/2)*100% confidence interval is computed
#   d     - half width of the confidence interval
#   lambda.grid - range of lambda values to try
#
# Details:
# if a single lambda value is givn then the sample size is calculated
# for this value. Otherwise, and following the suggestion in the paper
# all lambda values in lambda.grid going from 0 to 30 are tried
# and the lambda values resulting in the lowest sample size is used.
# Alternatively, one can use ((2*d)^(-1)+2*d)/2 as value.
#
# Returns:
# a vector containing the following four elements
#  n0     - sample size based on normal approximation
#  nstar  - sample size at most favorable lambda value in lambda.grid
#  cp     - coverage probability
#  lambda - value in lambda.grid giving the lowest nstar value
######################################################################

ciss.liubailey <- function(alpha, d, lambda.grid=0:30) {
  res <- sapply(lambda.grid, function(lambda) {
    liuSamSize(alpha=alpha, d=d, lambda=lambda)
  })

  minNstarIdx <- which.min(res["nstar",])
  c(res[,minNstarIdx],lambda=lambda.grid[minNstarIdx])
}

