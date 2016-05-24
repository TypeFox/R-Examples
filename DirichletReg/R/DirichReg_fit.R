DirichReg_fit <- function(Y, X, Z, sv, d, k, w, ctls, repar, base, vrb){

  n <- nrow(Y)
  npar <- length(sv)

  ops <- options(warn = -1L)
  on.exit(options(ops))

  ##############################################################################
  ############################################## alternative parametrization ###
  if(repar){
    k <- k[1L]
    beta_ind <- as.integer(matrix(seq_len(d*k), nrow = k, ncol = d)[,-base])
    beta_x_ind <- seq_len((d-1L)*k)
    gamma_ind <- seq.int((d-1L)*k + 1L, npar)
    ncolX <- ncol(X[[1L]])
    ncolZ <- ncol(Z)
    hessian.ind <- rbind(as.matrix(expand.grid(seq_len(k), seq_len(d)[-base])[,c(2L, 1L)]), cbind(-1L, seq_len(ncolZ)))

#    ## first a couple of iterations only for precisions to accelerate
#      bfgs1 <- maxBFGS(fn=DReg.repar,
#        start=sv, fixed=1:((d-1)*ncol(X[[1]])),
#        finalHessian=FALSE, iterlim=5L, tol=1e-2, reltol=1e-2,
#        logY=log(Y), X=X[[1L]], ncolX=ncolX, Z=Z, ncolZ=ncolZ, n=n, d=d, k=k, w=w, base=base, npar=npar, bi=beta_ind, bx=beta_x_ind, gi=gamma_ind, NR=FALSE)

      bfgs <- maxBFGS(fn=DReg.repar,
        start=sv,#bfgs1$estimate,
        finalHessian=FALSE, iterlim=ctls$iterlim, tol=ctls$tol1, reltol=ctls$tol1, print.level=ifelse(vrb == 0, 0, vrb - 1),
        logY=log(Y), X=X[[1L]], ncolX=ncolX, Z=Z, ncolZ=ncolZ, n=n, d=d, k=k, w=w, base=base, npar=npar, bi=beta_ind, bx=beta_x_ind, gi=gamma_ind, NR=FALSE)

      res <- maxNR(fn=DReg.repar,
        start=bfgs$estimate,
        iterlim=ctls$iterlim, tol=ctls$tol2, reltol=ctls$tol2, print.level=ifelse(vrb == 0, 0, vrb - 1),
        logY=log(Y), X=X[[1L]], ncolX=ncolX, Z=Z, ncolZ=ncolZ, n=n, d=d, k=k, w=w, base=base, npar=npar, bi=beta_ind, bx=beta_x_ind, gi=gamma_ind, NR=TRUE, h_dims=hessian.ind[,1L], h_vars=hessian.ind[,2L])

  ##############################################################################
  ################################################### common parametrization ###
  } else {
    seq_along_d <- seq_len(d)
    ncolX <- unlist(lapply(X, ncol))
    beta_x_ind <- lapply(seq_along_d, function(i){ seq.int(cumsum(c(0L, k))[i] + 1L, cumsum(k)[i]) })
    hessian.ind <- cbind(rep(seq_along_d, k), unlist(lapply(k, function(i){ seq_len(i) })))
    
      bfgs <- maxBFGS(fn=DReg,
        start=sv,
        finalHessian=FALSE, iterlim=ctls$iterlim, tol=ctls$tol1, reltol=ctls$tol1, print.level=ifelse(vrb == 0, 0, vrb - 1),
        logY=log(Y), X=X, ncolX=ncolX, n=n, d=d, k=k, w=w, npar=npar, seq_along_d=seq_along_d, bx=beta_x_ind, NR=FALSE)

      res <- maxNR(fn=DReg,
        start=bfgs$estimate,
        iterlim=ctls$iterlim, tol=ctls$tol2, reltol=ctls$tol2, print.level=ifelse(vrb == 0, 0, vrb - 1),
        logY=log(Y), X=X, ncolX=ncolX, n=n, d=d, k=k, w=w, npar=npar, seq_along_d=seq_along_d, bx=beta_x_ind, NR=TRUE, h_dims=hessian.ind[,1L], h_vars=hessian.ind[,2L])
  }

  res$bfgs.it <- bfgs$iterations

  return(res)

}
