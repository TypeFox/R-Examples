#######################################################################
## Function: ll.oprobit.unim()
## Author  : Jonathan Wand (jwand@latte.harvard.edu)
##
## Calculate -LL of ordered probit
##
## Created:   2002-05-29
## Modified:  $Date: 2004/07/07 04:08:44 $
## Revision:  $Revision: 1.4 $
## RCS-ID:    $Id: func.ll.oprobit.R,v 1.4 2004/07/07 04:08:44 jwand Exp $
##
##
## INPUT: 
##   Y:    matrix (n x n.cat) with T/F for selecting cases
##   Xb:   vector of means
##   se:   standard deviation of normal
##   tau:  matrix (n x (n.cat-1)) with cutpoints
##   n.cat:number of categories per question
##
## OUTPUT:
##   prob=FALSE: returns scalar, sum of log-likelihood
##   prob=TRUE : vector of log-likelihood values
##
#######################################################################

ll.oprobit <- function(y,xb,se,taus,n.cat,debug=0)
{
  n  <- length(xb)
  llik  <- rep(1, n)
  if (sum(!is.finite(se)) > 0) {
    cat("ll.oprobit BAD se:",se,"\n")
    return( -200 * n)
  }

  
  opfunc <- function(n,n.cat,se,penalty,y,xb,taus,llik) {

#    cat("ll.oprobit\n")
#    print( cbind(y,xb))
#    print(se)
#    cat("taus\n")
#    print(taus)
#    cat("N",n,n.cat,"\n")
  
    .C("opll",as.integer(n), as.integer(n.cat),
       as.double(se), as.double(penalty),
       as.integer(y),
       as.double(xb), as.double(taus),
       llik = as.double(llik),
       PACKAGE = "anchors"
       )$llik
  }
  
  opfunc(n,n.cat,se,-200,y,xb,taus,llik)

}
