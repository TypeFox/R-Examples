#######################################################################
## Function: gr.oprobit()
## Author  : Jonathan Wand (jwand@latte.harvard.edu)
##
## Calculate gradient of ordered probit
##
## Created:   2003-05-05
## Modified:  $Date: 2004/07/07 04:08:44 $
## Revision:  $Revision: 1.2 $
## RCS-ID:    $Id: func.gr.oprobit.self.C.R,v 1.2 2004/07/07 04:08:44 jwand Exp $
##
##
## INPUT: 
##   Y:    matrix (n x n.cat) with counts of choices
##   Xb:   vector of means
##   se:   standard deviation of normal
##   taus: matrix (n x (n.cat-1)) with cumulative cutpoints
##   n.cat:number of categories per question
##
## OUTPUT:
##   vector of gradients
#######################################################################

gr.oprobit.self.C <- function(Y,Xb,se,taus,V,V1,X,n.cat,n.self,do.test=0,verbose=FALSE)
{

  n    <- length(Xb)
  nvarX <- NCOL(X)
  nvarV <- NCOL(V)
  nvarV1<- NCOL(V1)

  if (verbose) {
    cat("GR.oprobit.self\n")
    cat("n.cat",n.cat,
        "se",se,
        "length taus",length(taus),
        "n",n,
        "nvarV",nvarV,
        "nvarV1",nvarV1,
        "nvarX",nvarX,
        "\n")
  }

  dLdBeta  <- rep(0.0, nvarX)
  dLdGamma <- rep(0.0, nvarV*(n.cat-2)*n.self)
  dLdGamma1 <- rep(0.0, nvarV1)
  dLdSigma <- 0.0


  z <- .C("opllgrself",
            as.integer(n), 
            as.integer(n.cat),  
            as.integer(nvarX),  
            as.integer(nvarV),
            as.integer(nvarV1),
            as.integer(n.self),
            as.double(se), 
            as.integer(Y), 
            as.double(Xb), 
            as.double(taus), 
            as.double(V),
            as.double(V1),
            as.double(X),
          dLdSigma = as.double(dLdSigma),
          dLdBeta  = as.double(dLdBeta),
          dLdGamma = as.double(dLdGamma),
          dLdGamma1 = as.double(dLdGamma1),
          PACKAGE = "anchors"
            );


  return( list( beta=z$dLdBeta ,
               sigma=z$dLdSigma,
               gamma=z$dLdGamma,
               gamma1=z$dLdGamma1))

}
