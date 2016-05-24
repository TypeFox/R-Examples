tsHessian <- function(param, fun, ...)
{
  ## Description:
  ##   Computes two sided (TS) approximated Hessian

  ## Source:
  ## This code was borrowed from the fBasics
  ## package, utils-hessian.R with slight modification
  ## A function borrowed from Kevin Sheppard's Matlab garch toolbox
  ## as implemented by Alexios Ghalanos in his rgarch package

  ## Notes:
  ##   optionally requires package Matrix (added as suggestion)


  ## FUNCTION:

  ## Settings:
  n <- length(param)
  fx <- fun(param, ...)
  eps <- .Machine$double.eps

  ## Compute the stepsize:
  h <- eps^(1/3) *
    apply( as.data.frame(param), 1, FUN = function(z) max(abs(z), 1.0e-2))
  paramh <- param + h
  h <- paramh - param
  ee <- diag(h) # Matrix(diag(h), sparse = TRUE)

  ## Compute forward and backward steps:
  gp <- vector(mode = "numeric", length = n)
  for(i in 1:n) gp[i] <- fun(param + ee[, i])
  gm <- vector(mode = "numeric", length = n)
  for(i in 1:n) gm[i] <- fun(param - ee[, i])
  H <- h %*% t(h)
  Hm <- H
  Hp <- H

  ## Compute double forward and backward steps:
  for(i in 1:n){
    for(j in  i:n){
      Hp[i, j] <- fun(param + ee[, i] + ee[, j])
      Hp[j, i] <- Hp[i, j]
      Hm[i, j] <- fun(param - ee[, i] - ee[, j])
      Hm[j, i] <- Hm[i, j]
    }
  }

  ## Compute the Hessian:
  for(i in 1:n){
    for(j in  i:n){
      H[i, j] <- ((Hp[i, j] - gp[i] - gp[j] + fx + fx - gm[i] - gm[j] +
                   Hm[i, j]) / H[i, j]) / 2
      H[j, i] <- H[i, j]
    }
  }

  ## Return Value:
  return(H)
}
