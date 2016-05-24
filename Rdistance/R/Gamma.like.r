Gamma.like <- function(a, dist, w.lo=0, w.hi=max(dist), series="cosine", expansions=0, scale = TRUE){
#
#   Compute gamma likelihood.
#
#   Input:
#       a = vector of model coefficient and r gamma shape parameter
#       dist = vector of distance observations
#       w = strip half width
#       series = type of expansion to apply, does not apply to gamma yet
#       expansions = number of expansion terms, does not apply to gamma yet
#       scale = whether to scale the likelihood.
#
#   See comments in uniform.like.r
#

    r <- a[1]
    lam <- a[2]
    

    if( r <= 1 ) warning("Shape parameter of gamma likelihood invalid (<= 1).")
    if( lam < 0 ) warning("Scale parameter of gamma likelihood invalid (< 0).")

    dist[ (dist < w.lo) | (dist > w.hi) ] <- NA
    eta <- 0  # when have no covariates


    # This was all part of Becker and Quan's code
    #bb <- (1/gamma(r)) * (((r - 1)/exp(1))^(r - 1))
    #J <- dim(X.)[2]
    #eta <- X. %*% matrix(b2[1:(length(b2)-1)],ncol=1)
    #eta <- 0   # use if no covariates ?
    #lam <- exp(eta)
    #w1b <- wb/(lam * bb)
    #v1 <- dist/(lam * bb)
    #w1 <- w/(lam * bb)
    #loglik <- (r - 1) * log(v1) - v1 - lgamma(r) - log(bb) - eta - log(pgamma(w1, r))
    #if (USE.PGAMMA.WB){
    #    loglik <- (r - 1) * log(v1) - v1 - lgamma(r) - log(bb) - eta - log(pgamma(w1, r)-pgamma(w1b,r))
    #} else {
    #    #computationaly more efficient form of likelihood than eq. 8 in Becker and Quang 2008
    #    #note that when wb is zero pgamma(w1b,r) will also be zero so
    #    #     no blind strip in intergration
    #    loglik <- (r - 1) * log(v1) - v1 - lgamma(r) - log(bb) - eta - log(pgamma(w1, r))
    #}


    #   In the following, I am trying to avoid taking gamma(r), which could be huge   ## DON'T KNOW WHETER THIS WORKS
    #log.b <- -lgamma(r) + (r - 1)*(log(r - 1) - 1)
    #log.v1 <- log(dist) - log(lam) - log.b
    #v1 <- exp( log.v1 )
    #loglik <- (r - 1) * log.v1 - v1 - lgamma(r) - log.b - eta
    #like <- exp(loglik)

    #   In the following, I assume I can evaluate gamma(r)  ## THIS WORKS
    #b <- (1/gamma(r)) * (((r - 1)/exp(1))^(r - 1))
    #v1 <- dist/(lam * b)
    #w1 <- w1/(lam * b)
    #loglik <- (r - 1) * log(v1) - v1 - lgamma(r) - log(b) - eta # - log(pgamma(w1, r))  # this last term is meant to be integration constant, but it's off, and we use numerical integration below
    #like <- exp(loglik)

    #   In the following, I use the built in R function  ## THIS ALSO WORKS
    #   Note: I think Quan is missing an extra 1/lam in front of his density equation.
    b <- (1/gamma(r)) * (((r - 1)/exp(1))^(r - 1))
    like <- dgamma( dist, shape=r, scale=lam*b )

    if( scale ){
         like = like / integration.constant(Gamma.like, w.lo=w.lo, w.hi=w.hi, a=a,series=series,expansions=expansions)
    }

    like

}
