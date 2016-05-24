gr.oprobit.vign <- function(Y,Xb,se,taus,V,n.cat,do.test=0,verbose=FALSE)
{
#######################################################################
## Function: gr.oprobit()
## Author  : Jonathan Wand (jwand@latte.harvard.edu)
##
## Calculate gradient of ordered probit
##
## Created:   2003-05-05
## Modified:  $Date: 2004/07/07 04:08:44 $
## Revision:  $Revision: 1.2 $
## RCS-ID:    $Id: func.gr.oprobit.vign.R,v 1.2 2004/07/07 04:08:44 jwand Exp $
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

  options(warn=1)
  
  n    <- length(Xb)
  nvar <- NCOL(V)

  taus <- matrix(taus,nrow=n)

  if (verbose) {
    cat("GR.oprobit.vign\n")
    cat("n.cat",n.cat,
        "se",se,
        "dim taus",dim(taus),
        "n",n,
        "nvar",nvar,
        "dim V",dim(V),
        "\n")
  }
  
  vmat  <- rep(0,n)
  pmat  <- rep(0,n)
  dmat  <- rep(0,n)
  for (k in 1:(n.cat-1)) {
    pmat  <- cbind( pmat, pnorm((taus[,k]-Xb)/se))
    dmat  <- cbind( dmat, dnorm((taus[,k]-Xb)/se))
    vmat  <- cbind( vmat,       (taus[,k]-Xb)/se )
  }
  vmat  <- cbind( vmat, rep(0,n))
  pmat  <- cbind( pmat, rep(1,n))
  dmat  <- cbind( dmat, rep(0,n))

  dLdGamma <- matrix(0, n, (n.cat-1)*nvar)
  dLdSigma <- dLdTheta <- 0
  for (i in 1:n) {
    j   <- Y[i] + 1
    if (j == 1)
      next

    denom <- pmat[i,j] - pmat[i,j-1]

    dLdTheta <- dLdTheta + (dmat[i,j-1] - dmat[i,j]) / denom
    dLdSigma <- dLdSigma + (dmat[i,j-1] * vmat[i,j-1] - dmat[i,j] * vmat[i,j])/ denom

#    if (0) {
#      cat("DD",dmat[i,j-1],dmat[i,j],"\n")
#      cat("PP",pmat[i,j-1],pmat[i,j],"\n")
#      cat("RR",rrj,rrj1,"\n")
#      cat("V\n")
#      print(V[i,])
#    }

    m1 <- 1
    if (Y[i] == n.cat) 
      m1 <- 0
    for (k in 1:(n.cat-1)) {
      idx <- 1:nvar + (k-1)*nvar
      if (k <= (j-1)) {
        dLdGamma[i,idx] <- 
          (m1 * V[i,] * dmat[i,j] -
           as.numeric(k <= (j-2)) * V[i,] * dmat[i,j-1]) / denom
        
      }
    }
  }
  dLdTheta <- sum(dLdTheta) * 1/se
  dLdSigma <- sum(dLdSigma) * 1/se
  dLdGamma <- apply(dLdGamma,2,sum) * 1/se

  
  return(list(theta=dLdTheta,
         sigma=dLdSigma,
         gamma=dLdGamma))

}
