DReg <- function(x, logY, X, ncolX, n, d, k, w, npar, seq_along_d, bx, NR, h_dims, h_vars){

  B <- lapply(bx, function(b_ind){ x[b_ind] })

  A <- matrix(unlist(lapply(seq_along_d, function(i){ exp(X[[i]] %*% B[[i]]) })), nrow = n, ncol = d)
  Aplus <- .rowSums(A, n, d, FALSE)

  digamma_A  <- digamma(A)
  trigamma_A <- trigamma(A)
  digamma_Aplus  <- digamma(Aplus)
  trigamma_Aplus <- trigamma(Aplus)



################################################################################
### LOG-LIKELIHOOD & GRADIENT ##################################################
################################################################################

  LL <- .Call("wght_LL_grad_common", logY, A, Aplus, digamma_A, digamma_Aplus, X, ncolX, c(n, d), npar, w)



################################################################################
### HESSIAN ####################################################################
################################################################################

  if(NR){

    hessian <- matrix(NA_real_, nrow = npar, ncol = npar)

    for(hess.j in seq_len(npar)){
      for(hess.i in seq_len(npar)){
        if(hess.i < hess.j) next

        derv <- h_dims[c(hess.i, hess.j)]

        vars <- h_vars[c(hess.i, hess.j)]

        ########################################################################
        ##################################################### SAME RESPONSES ###
        if(derv[1L] == derv[2L]) {
          derv <- derv[1L]

          hessian[hess.i, hess.j] <- hessian[hess.j, hess.i] <-
          sum(w*(
            X[[derv]][,vars[1L]] * X[[derv]][,vars[2L]] * A[,derv] * (
            logY[,derv] + digamma_Aplus - digamma_A[,derv] + A[,derv] * (
              trigamma_Aplus - trigamma_A[,derv]
              )
            )
          ))
        ########################################################################
        ################################################ DIFFERENT RESPONSES ###
        } else {
          hessian[hess.i, hess.j] <- hessian[hess.j, hess.i] <-
          sum(w*(
            X[[derv[1L]]][,vars[1L]]*X[[derv[2L]]][,vars[2L]]*A[,derv[1L]]*A[,derv[2L]]*trigamma_Aplus
          ))
        }
      }
    }

    attr(LL, "hessian") <- hessian

  }

  return(LL)

}
