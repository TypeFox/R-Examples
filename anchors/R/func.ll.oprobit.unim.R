ll.oprobit.unim <- function(Y,Xb,se,tau,n.cat,do.test=0,get.prob=FALSE)
{
#######################################################################
## Function: ll.oprobit.unim()
## Author  : Jonathan Wand (jwand@latte.harvard.edu)
##
## Calculate -LL of ordered probit using grouped data
##
## Created:   2002-05-29
## Modified:  $Date: 2002/10/17 14:44:48 $
## Revision:  $Revision: 1.2 $
## RCS-ID:    $Id: func.ll.oprobit.unim.R,v 1.2 2002/10/17 14:44:48 jwand Exp $
##
##
## INPUT: 
##   Y:    matrix (n x n.cat) with counts of choices
##   Xb:   vector of means
##   se:   standard deviation of normal
##   tau:  matrix (n x (n.cat-1)) with cutpoints
##   n.cat:number of categories per question
##
## OUTPUT:
##   returns scalar, sum of log-likelihood
#######################################################################

  n  <- length(Xb)
  if (n > 1) {
    taus <- row.cumsum(tau)
  } else {
    taus <- matrix(cumsum(tau),nrow=1)
  }
  
  pmat  <- NULL
  p.lag <- rep(0,n)

  if (do.test > 0)
    cat("ll.oprobit.unim: calc prob\n")
  for (jj in 1:(n.cat-1)) {
    p.tmp <- pnorm(taus[,jj],Xb,se)
    pmat  <- cbind( pmat, p.tmp-p.lag )
    p.lag <- p.tmp
  }
  ## CAUTION! going to reuse 'p.lag' from 2 lines above...
  if (do.test > 0)
    cat("ll.oprobit.unim: cbind\n")
  pmat <- cbind(pmat, 1-p.lag)
  if (do.test > 0)
    cat("ll.oprobit.unim: penalty\n")
  ## Y should have same dimensions as pmat!
  pmat[pmat==0] <- exp(-100)
  tmp <- log( pmat  ) * Y
  llik  <- sum(tmp)

  if (do.test > 1) {
    cat("ll.oprobit.unim(): TEST\n")
    cat("LLIK",llik,"SE",se,"N.CAT",n.cat,"\n")
    cat("XB\n")
    print(Xb)
    print("TAUs\n")
    print(taus)
    cat("PMAT2\n")
    print(pmat)
    print("Y")
    print(Y)
    print("LLi")
    print(tmp)
  }
  
  if (!get.prob) {
    return( llik )
  } else {
    return( pmat)
  }
}
