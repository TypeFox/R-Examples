#################################################################################
####This function calculates the empricial cdf given a vector of data and an evaluation point
################################################################################
empirical.cdf <-
function(xx,c0)
  {

    #xx <- sort(xx)##no need for sorting to obtain Pr(xx<=c0)
    nobs <- length(xx)
    prob0 <- sum(xx<=c0)/nobs
    return(prob0)
  }

