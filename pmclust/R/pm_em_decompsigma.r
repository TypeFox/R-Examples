### This file contains functions for decomposition of SIGMA of MVN.
### These will be majorly used in m.step() to update PARAM$U and PARAM$U.check.

decompsigma <- function(SIGMA){
  U <- chol(SIGMA)
  u <- diag(U)
  U.check <- TRUE
  if(any(u < .pmclustEnv$CONTROL$U.min | u > .pmclustEnv$CONTROL$U.max)){
    U.check <- FALSE
  }

  list(value = U, check = U.check)
} # End of decompsigma().

