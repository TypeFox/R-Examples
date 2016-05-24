#############
#### Parameters simulator
#############
# Given a mean and a covariance generates nsim multivariate normal rand variables
# and check constraints.
#
# 
# ARGS:
# theMean = (numeric) mean vector
# covar   = (matrix) covariance
# nsim    = (scalar integer) number of simulations
# constr  = (named list) of 3 elements
#           [["indexes"]] = (numeric integers) indexes of the elements to check
#           [["upper"]]   = (numeric) upper bounds for the elements in "indexes"
#           [["lower"]]   = (numeric) lower bounds for the elements in "indexes"
#
# OUT:
# A (nsim by length(theMean)) matrix were each row has been simulate from a MVN
# and is inside the constraints.

.paramsSimulator <- function(theMean, covar, nsim, constr = list())
{
  stopifnot( is.vector(theMean), is.matrix(covar), length(theMean) == ncol(covar) )
  
  # If a parameter has variance 0 in the proposal we save it's index in "fixPar"
  # we modidy the covariance and we save the initial covariance
  fixPar <- which( diag(covar) == 0 )
  anyFix <- ( length(fixPar) > 0 )
  if(anyFix) diag(covar)[fixPar] <- 0.1
  
  # Simulating the parameters, equivalent of "output <- mvrnorm(n = nsim, mu = theMean, Sigma = covar)"
  # but this code is no good because difference in covariance of order 1e-16 make it simulate different
  # random variables
  
  nPar <- length(theMean)
  cholFact <- chol(covar)
  output <- .rmvn(nsim, mu = theMean, sigma = cholFact, isChol = TRUE) 
  
  # (Optionally) check contraints and re-simulate parameter vectors that fall outside them
  if( length(constr) )
  {
    upper <- constr$upper
    lower <- constr$lower
    indexes <- constr$indexes
    
    stopifnot( length(upper) == length(lower), length(lower) <= nPar, length(indexes) == length(lower), 
               all(upper > lower) )
    
    # Calls C++ function that loops through output, and checks the constrainst
    output <- .Call("checkBoundsCpp", 
                    theMean_ = theMean, cholFact_ = t(cholFact), indexes_ = indexes,            
                    upper_ = upper, lower_ = lower, output_ = output, 
                    PACKAGE = "synlik")
  }
  
  # Resetting the parameters that are fixed.
  if(anyFix) output[ , fixPar] <- matrix(theMean[fixPar], nrow(output), length(fixPar), byrow = TRUE)
  
  return(output)
}


####
# Test
####
# 
# currPar <- c(0.5, 1.5)
# nsim <- 1000
# constr <- list("indexes" = c(1, 2), "upper" = c(1, 2), "lower" = c(0, 1))
# covar <- matrix(c(2, 1, 1, 2), 2, 2)
# 
# out1 <-  synlik:::.rmvn(nsim, currPar, covar)
# output <- synlik:::.paramsSimulator(theMean =  currPar, covar, nsim, constr)
# 
# plot(out1)
# points(output, col = 2)

