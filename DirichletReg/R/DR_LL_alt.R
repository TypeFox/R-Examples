DReg.repar <- function(x, logY, X, ncolX, Z, ncolZ, n, d, k, w, base, npar, bi, bx, gi, NR, h_dims, h_vars){
################################################################################
### PREPARATION ################################################################
################################################################################

  B <- matrix(0.0, nrow = k, ncol = d)
  B[bi] <- x[bx]

  eps  <- apply(B, 2L, function(b){ exp(X %*% b) })

  mu <- eps / .rowSums(eps, n, d, FALSE)

  phi <- as.numeric(exp( Z %*% matrix(x[gi], ncol=1L) ))

  A <- mu * phi

  digamma_A <- digamma(A)
  trigamma_A <- trigamma(A)
  digamma_phi <- digamma(phi)



################################################################################
### LOG-LIKELIHOOD & GRADIENT ##################################################
################################################################################

  LL <- .Call("wght_LL_grad_alternative", logY, A, mu, phi, digamma_A, digamma_phi, X, ncolX, Z, ncolZ, n, d, base, npar, w)



################################################################################
### HESSIAN ####################################################################
################################################################################

  if(NR){

  trigamma_phi <- trigamma(phi)

  hessian <- matrix(NA_real_, nrow=npar, ncol=npar)

  for(hess.j in seq_len(npar)){
    for(hess.i in seq_len(npar)){
      if(hess.i < hess.j){ next }

      v1 <- h_vars[hess.i]
      v2 <- h_vars[hess.j]

      derv <- h_dims[c(hess.i, hess.j)]
      d1 <- derv[1L]
      d2 <- derv[2L]

      ##########################################################################
      ############################################### BETAs - SAME RESPONSES ###
      if((derv[1L] == derv[2L]) & all(derv != -1L)) {
        derv <- derv[1L]

        hessian[hess.i, hess.j] <- hessian[hess.j, hess.i] <-
          sum(w*(
            X[,v1] * X[,v2] * A[,derv] * (
              (2.0*mu[,derv] - 1.0) * (
                rowSums(mu[,-derv,drop=FALSE]*(logY[,-derv,drop=FALSE]-digamma_A[,-derv,drop=FALSE]))
              - (1.0 - mu[,derv]) * (logY[,derv]-digamma_A[,derv])
              )
              -
              A[,derv]*(
                (1.0 - mu[,derv])^2 * trigamma_A[,derv]
              + rowSums(mu[,-derv,drop=FALSE]^2*trigamma_A[,-derv,drop=FALSE])
              )
            )
          ))
      ##########################################################################
      ########################################## BETAs - DIFFERENT RESPONSES ###
      } else if((derv[1L] != derv[2L]) & all(derv != -1L)) {
        hessian[hess.i, hess.j] <- hessian[hess.j, hess.i] <-
          sum(w*(
            X[,v1] * X[,v2] * mu[,d1] * mu[,d2] * phi * (
              rowSums(
                mu[,-derv,drop=FALSE] * (
                  2.0 * (logY[,-derv,drop=FALSE] - digamma_A[,-derv,drop=FALSE])
                - A[,-derv,drop=FALSE] * trigamma_A[,-derv,drop=FALSE]
                )
              )
            +
              rowSums(
                (2*mu[,derv,drop=FALSE] - 1.0) * (logY[,derv,drop=FALSE] - digamma_A[,derv,drop=FALSE])
              - A[,derv,drop=FALSE] * (mu[,derv,drop=FALSE] - 1.0) * trigamma_A[,derv,drop=FALSE]
              )
            )
          ))
      ##########################################################################
      ######################################################### BETA / GAMMA ###
      } else if(any(derv != -1L) & any(derv == -1L)) {
        derv <- derv[which(derv != -1L)]

        hessian[hess.i, hess.j] <- hessian[hess.j, hess.i] <-
          sum(w*(
            Z[,v1] * X[,v2] * A[,derv] * (
              rowSums(mu[,-derv,drop=FALSE] * (
                digamma_A[,-derv,drop=FALSE] + A[,-derv,drop=FALSE]*trigamma_A[,-derv,drop=FALSE] - logY[,-derv,drop=FALSE]
              ))
            +
              (mu[,derv] - 1.0) * (
                digamma_A[,derv] + A[,derv]*trigamma_A[,derv] - logY[,derv]
              )
            )
          ))
      ##########################################################################
      ############################################################### GAMMAs ###
      } else if(all(derv == -1)){
        hessian[hess.i, hess.j] <- hessian[hess.j, hess.i] <-
          sum(w*(
            Z[,v1] * Z[,v2] * phi * ( digamma_phi + phi * trigamma_phi + rowSums(mu * (logY - digamma_A - A * trigamma_A)) )
          ))
      }
    }
  }

  attr(LL, "hessian") <- hessian

  }

  return(LL)

}
