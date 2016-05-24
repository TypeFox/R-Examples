gr.oprobit.self <- function(Y,Xb,se,taus,V,X,n.cat,do.test=0,verbose=FALSE)
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
## RCS-ID:    $Id: func.gr.oprobit.self.R,v 1.2 2004/07/07 04:08:44 jwand Exp $
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
  nvarX <- NCOL(X)
  nvarV <- NCOL(V)

  taus <- matrix(taus,nrow=n)

  if (verbose) {
    cat("GR.oprobit.self\n")
    cat("n.cat",n.cat,
        "se",se,
        "dim taus",dim(taus),
        "n",n,
        "nvarV",nvarV,
        "nvarX",nvarX,
        "\n")
  }

  vmat  <- rep(0,n)
  pmat  <- rep(0,n)
  dmat  <- rep(0,n)
  for (k in 1:(n.cat-1)) {
    pmat  <- cbind( pmat, pnorm((taus[,k]-Xb)/se))
    dmat  <- cbind( dmat, dnorm((taus[,k]-Xb)/se))
    vmat  <- cbind( vmat,       (taus[,k]-Xb)/se )
#    pmat  <- cbind( pmat, dnorm(taus[,k],Xb,se))
#    dmat  <- cbind( dmat, pnorm(taus[,k],Xb,se))
#    vmat  <- cbind( vmat,      (taus[,k]-Xb)/se)
  }
  vmat  <- cbind( vmat, rep(0,n))
  pmat  <- cbind( pmat, rep(1,n))
  dmat  <- cbind( dmat, rep(0,n))

  print( cbind(as.double(pmat),
               as.double(dmat),
               as.double(vmat),
               as.double(Xb)))
  
  dLdBeta  <- matrix(0, n, nvarX)
  dLdGamma <- matrix(0, n, (n.cat-1)*nvarV)
  dLdSigma <- 0
  for (i in 1:n) {
    j   <- Y[i] + 1
    if (j == 1)
      next

    denom <- pmat[i,j] - pmat[i,j-1]
    cat("DENOM",i,denom,"\n")

    dLdSigma <- dLdSigma + (dmat[i,j-1] * vmat[i,j-1] - dmat[i,j] * vmat[i,j])/ denom

    dLdBeta[i,] <- X[i,] * ( dmat[i,j-1] - dmat[i,j]) / denom


    m1 <- 1
    if (Y[i] == n.cat) 
      m1 <- 0
    for (k in 1:(n.cat-1)) {
      idx <- 1:nvarV + (k-1)*nvarV
      if (k <= (j-1)) {
        dLdGamma[i,idx] <- 
          (m1 * V[i,] * dmat[i,j] -
           as.numeric(k <= (j-2)) * V[i,] * dmat[i,j-1]) / denom
        
      }
    }
    
  }

  dLdSigma <- sum(dLdSigma)         * 1/se
  dLdGamma <- apply(dLdGamma,2,sum) * 1/se
  dLdBeta  <- apply(dLdBeta ,2,sum) * 1/se

  
  return( list(beta=dLdBeta ,
         sigma=dLdSigma,
         gamma=dLdGamma))

}
