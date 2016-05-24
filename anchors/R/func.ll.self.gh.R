#######################################################################  
## Function: ll.self.gh()
## Author:   Jonathan Wand
## Date:     2002-05-29, updated 06-04
##
## Purpose: - Calculate -LL of sequence of ordered probits
##          with common means which share a random effect
##          but different cutpoints
##          - RE integrated out using fixed point gaussian hermite integration
##
## INPUT:
##   y0:      matrix of category counts associated with each unique obs
##   xb:      vector of means 
##   ln.sigma.y:  variance for normal distributions
##   ln.sigma.re: variance for random effect
##   tau:     matrix (n x n.self*(n.cat-1)) with cut points 
##   gh:      dataframe with GH values and weights
##   n.cat:   number of categories per oprobit
##
## OUTPUT:
##   returns scalar, sum of log-likelihood
##
#######################################################################

ll.self.gh <- function(y0,xb,sigma.y,sigma.re,taus,gh,n.cat,do.print=0) {

  tidx <- function(i.self,i.cat) {
    (n.cat-1)*(i.self-1)+i.cat
  }
  
#  Rprof("self.gh",append=TRUE,interval=0.005)

  se     <- sigma.y[1]
  su     <- sigma.re
  rho    <- su^2 / (se^2+su^2)
  s2su   <- sqrt(2*rho/(1-rho));

  dd <- dim(y0)
  n     <- dd[1]
  n.self<- dd[2]
  llik  <- rep(0, n)

  
  if (sum(!is.finite(se)) > 0) {
    cat("BAD se:",se,"\n")
    return( -200 * n)
  }
  if (sum(!is.finite(su)) > 0) {
    cat("BAD su:",su,"\n")
    return( -200 * n)
  }
  if (sum(!is.finite(s2su)) > 0) {
    cat("BAD s2su:",s2su,"\n")
    return( -200 * n)
  }

  opfunc <- function(n,n.self,n.cat,s2su,se,penalty,y0,xb,taus,llik,gh) {
    #cat("LENGTHS xb",length(as.double(xb)),
    #    "taus",length(as.double(taus)),
    #    "y0",length(as.double(y0)),"\n")
    .C("opllgh",as.integer(n), as.integer(dim(gh)[1]) , as.integer(n.self), as.integer(n.cat),
       as.double(s2su), as.double(se), as.double(penalty),
       as.double(gh[,1]), as.double(gh[,2]),
       as.integer(y0),
       as.double(xb), as.double(taus),
       llik = as.double(llik),
       PACKAGE = "anchors"
       )$llik
  }
  
  llik <- opfunc(n,n.self,n.cat,s2su,se,-200,y0,xb,taus,llik,gh)
  rv <- sum(llik)
#  cat("rv",rv,"\n")
#  print(llik)
  
#  Rprof(NULL)
  return( rv )
}
