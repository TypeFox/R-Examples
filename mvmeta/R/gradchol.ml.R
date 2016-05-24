###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2014
#
gradchol.ml <-
function(par,U,ind1,ind2,invSigmalist,reslist,nalist,k,m) {
#
################################################################################
# FUNCTION TO COMPUTE THE MATRIX DERIVATIVES IN TERMS OF PARAMETERS OF
# THE CHOLESKY DECOMPOSITION
#
  grad <- sapply(seq(length(par)),function(i) {
    # COMPUTE THE DERIVATIVE OF Psi IN TERMS OF ITS CHOLESKY DECOMPOSITION U
    A <- B <- C <- diag(0,k)
    A[ind2[i],] <- B[,ind2[i]] <- U[ind1[i],]
    C[ind2[i],] <- C[,ind2[i]] <- 1
    D <- C*A+C*B
    # COMPUTE THE GRADIENT
    gr <- sum(mapply(function(invSigma,res,na) {
      E <- crossprod(res,invSigma)%*%D[!na,!na]%*%invSigma%*%res
      F <- sum(diag(invSigma%*%D[!na,!na]))
      return(as.numeric(0.5*(E-F)))},invSigmalist,reslist,nalist))
  })
#
  grad
}

#
