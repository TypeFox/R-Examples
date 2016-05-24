momSE <- function(order = 4, n, mom)
{
  ## Purpose: Calculate standard errors of sample moments
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## order    order of the moment whose SE is required
  ## n        sample size
  ## mom      moments of the distribution from which sample is drawn
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date:  4 Feb 2010, 15:13

  if (!is.vector(mom))
    stop("A vector of moments must be supplied")
  if (!all.equal(length(mom), 2*order))
    stop("2*order moments must be supplied")

  if (order == 2){
    stErr <- sqrt((mom[4] - mom[2]^2)/n)
  } else if (order == 3){
    stErr <- sqrt((mom[6] - mom[3]^2 - 6*mom[4]*mom[2] + 9*mom[2]^3)/n)
  } else if (order == 4){
    stErr <-
      sqrt((mom[8] - mom[4]^2 - 8*mom[5]*mom[3] + 16*mom[2]*mom[3]^2)/n)
  }
    
  return(stErr) 
}
