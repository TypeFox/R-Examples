ll.self.ai <- function(y,xb,sigma.y,sigma.re,tau,n.cat,vlist=NULL,
                       do.print=0) {

#######################################################################
## Function: ll.self.ai()
## Author:   Jonathan Wand
##
## Calculate -LL of sequence of ordered probits
##   with common means which share a random effect
##   but different cutpoints
##   => RE integrated out using adaptive integration
##
## Created:  2002-05-29
## Modified:  $Date: 2003/12/12 20:42:33 $
## Revision:  $Revision: 1.3 $
## RCS-ID:    $Id: func.ll.self.ai.R,v 1.3 2003/12/12 20:42:33 jwand Exp $
##
## INPUT:
##   Y:      n-vector of category chosen by each case
##   xb:     n-vector of means
##   sigma.y:  variance for normal distributions
##   sigma.re: variance for random effect
##   tau:    matrix (n x n.self*(n.cat-1)) with cut points 
##   n:      number of cases
##   n.self: number of ordered probits
##   n.cat:  number of categories per oprobit
##
## OUTPUT:
##   returns scalar, sum of log-likelihood
##
#######################################################################  
  #cat("Entering ll.self.ai()...\n")
  dd <- dim(y)
  n     <- dd[1]
  n.self<- dd[2]
        
  se     <- sigma.y[1]
  su     <- sigma.re

  if (is.null(vlist)) {
    int.lb <- -Inf
    int.ub <- Inf
    int.tol<- 1e-12
    int.sub<- 100
  } else {
    int.lb <- vlist$int.range[1]
    int.ub <- vlist$int.range[2]
    int.tol<- vlist$int.tol
    int.sub<- vlist$int.sub
  }
  
  pp <- rep(NA,n)
  ## loop over cases
  for (i in 1:n) {
    ## make a function to that takes one argument (a vector)
    ## and returns a vector of same length
    ## !! integrate passes multiple points at which to evaluate each call !!
    f <- function(eta) {
      p <- ll.self.ai.prob(y[i,],xb[i,]+eta,se,su,tau[i,],n.self,n.cat)
      return( dnorm(eta,0,su) * apply(p,1,prod) )
    }
    ## go get expected probability
    pp[i] <- integrate(f,int.lb,int.ub,subdivisions=int.sub,
                       rel.tol=int.tol,abs.tol=int.tol     )$value
  }
  
  #cat("Exiting ll.self.ai()...\n")
  return(sum(log(pp)))
}


ll.self.ai.prob <- function(y,xb,se,su,tau,n.self,n.cat) {
#######################################################################
##
## Function: ll.self.ai.prob()
## Author:   Jonathan Wand
##
## Calculate probability of choosing a category.
##
## Created:   2002-05-29
## Modified:  $Date: 2003/12/12 20:42:33 $
## Revision:  $Revision: 1.3 $
## RCS-ID:    $Id: func.ll.self.ai.R,v 1.3 2003/12/12 20:42:33 jwand Exp $
##
## Input:
##   Y:      n-vector of category chosen by each case
##   xb:     n-vector of means
##   se:     oprobit sd
##   su:     RE sd
##   tau:    matrix (n x (n.cat-1)) with cutpoints
##   n:      number of cases
##   n.self: number of ordered probits
##   n.cat:  number of categories per oprobit
##
## Output:
##   returns matrix of probabilities
##
#######################################################################
  
  tidx <- function(i.self,i.cat) {
    (n.cat-1)*(i.self-1)+i.cat
  }

  p <- NULL
  for (i in 1:n.self){

    taus <- cumsum( tau[ tidx(i,1):tidx(i,(n.cat-1)) ] )
    
    j <- y[i]
    if (j==1|j>(n.cat-1)) {
      if (j==1)
        p <- cbind(p,pnorm(taus[1],xb,se))
      else 
        p <- cbind(p,pnorm(taus[n.cat-1],xb,se,lower.tail=FALSE))
    }
    else
      p <- cbind(p, pnorm(taus[j],xb,se)-pnorm(taus[j-1],xb,se))
  }
  return(p)
}
