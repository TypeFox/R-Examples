dglm.constant <- function(y,family,weights=1) {
  #  Constant term appearing in glm log-likelihood
  #  Used by dglm
  #  GKS  6 Jan 98, 4 Jul 98, 23 Sep 99.
  # "Binomial" changed to "binomial": PKD 05 Sep 2006
  #
  #  Exact cases (in binomial case, exact for phi near 1)
  #
  
  if ( family$family == "Tweedie") { 
    tweedie.p <- get("tweedie.p", envir = parent.frame() )
  }
  const <- switch(family$family[1],
                  "Gaussian" = length(y)*log(2*pi),
                  "Poisson" = 2*sum(y-y*ifelse(y>0,log(y),0)+lgamma(y+1)),
                  "Gamma" = 2*sum(log(y)),
                  "Inverse Gaussian" = sum(log(2*pi*y^3)),
                  "Tweedie" = switch(match( tweedie.p,c(0,1,2,3),nomatch=0),
                                     length(y)*log(2*pi),
                                     2*sum(y-y*ifelse(y>0,log(y),0)+lgamma(y+1)),
                                     2*sum(log(y)),
                                     sum(log(2*pi*y^3)) ),
                  "binomial" = -2*sum(lgamma(weights+1)-lgamma(weights*y+1)-
                                        lgamma(weights*(1-y)+1)+weights*(y*ifelse(y>0,log(y),0)+
                                                                           (1-y)*ifelse(1-y>0,log(1-y),0)))+sum(log(weights))
  )
  #
  #  Saddle-point approximation
  #
  if (is.null(const)) {
    V <- family$variance(y)
    if (any(V == 0)) V[V == 0] <- family$variance(y[V == 0]+1/6)
    const <- sum(log(2*pi*V))
    if (length(V) == 1 && length(y) > 1) const <- length(y)*const
  } ### END: if (is.null(const))
  const
} ### END: dglm.constant()

