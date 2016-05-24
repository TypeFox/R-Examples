# Bayesian Bootstrap
# Version:       0.1-6
# Date:     2011-02-24
# Author:         F.M.

BayesBoot <- function(ind.obs,...)
  {
    n.obs <- length(ind.obs)
    draw <- runif(n.obs-1,0,1)
    ## n.obs-1 random draws from a [0,1]uniform distribution
    diff.a <- diff.b <- c()
    diff.a[1:(n.obs-1)] <- sort(draw)
    diff.a[n.obs] <- 1
    diff.b[1] <- 0
    diff.b[2:n.obs] <- diff.a[1:(n.obs-1)]
    ## this creates two lists: list A has 1 as n.obs-th observation and  
    ## list B has 0 as first observation. The differences give a list of 
    ## n.obs probabilities which sum up to 1.
    p.draw <- diff.a - diff.b
    d.w.repl.obs <- rmultinom(1, size = n.obs, prob = p.draw)
    BB.ind.obs <- rep(ind.obs, d.w.repl.obs)
    return(BB.ind.obs)
  }
