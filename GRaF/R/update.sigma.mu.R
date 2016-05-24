update.sigma.mu <-
  function(Sigma, ds2, i, tnu) {
    # this takes 70% of total time in gpml Matlab code
    # tried re-coding in C++, but all the cost is the matrx algebra,
    # so little improvement. The approach is inherently slow.
    
    # get column i
    si <- Sigma[, i]
    # recalculate Sigma
    Sigma <- Sigma - ds2 / (1 + ds2 * si[i]) * si %*% t(si)
    # and mu
    mu <- Sigma %*% tnu
  
    return (list(Sigma = Sigma,
                 mu = mu))
}