## rwish delivers a pseudo-random Wishart deviate
##
## USAGE:
##
##   A <- rwish(v, S)
##
## INPUT:
##
##   v    degrees of freedom
##
##   S    Scale matrix
##
## OUTPUT:
##
##  A     a pseudo-random Wishart deviate
##
## Based on code originally posted by Bill Venables to S-news
## on 6/11/1998
##
## KQ on 2/5/2001
##
## COPIED: verbatim from MCMCpack fro the monomvn package


"rwish" <-
  function(v, S) {
    if (!is.matrix(S))
      S <- matrix(S)
    if (nrow(S) != ncol(S)) {
      stop(message="S not square in rwish().\n")
    }
    if (v < nrow(S)) {
      stop(message="v is less than the dimension of S in rwish().\n")
    }
    p <- nrow(S)
    CC <- chol(S)
    Z <- matrix(0, p, p)
    diag(Z) <- sqrt(rchisq(p, v:(v-p+1)))
    if(p > 1) {
      pseq <- 1:(p-1)
      Z[rep(p*pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p*(p-1)/2)
    }
    return(crossprod(Z %*% CC))
  }
                                                                                                 
