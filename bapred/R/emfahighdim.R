emfahighdim <-
function(eps, nbf, minerr = 1e-06, maxiter=100) {

  n <- nrow(eps)
  p <- ncol(eps)

  Sdiag <- apply(eps, 2, var)*(1 - 1/n)

  eig <- svd((1/sqrt((n - 1))) * t(eps))

  if(nbf <= n) {
    evectors <- eig$u[, 1:nbf]
    evalues <- eig$d^2
    if (nbf > 1) 
      B <- evectors[, 1:nbf] %*% diag(sqrt(evalues[1:nbf]))
    if (nbf == 1) 
      B <- matrix(evectors, nrow = p, ncol = 1) * sqrt(evalues[1])
  }
  else
     B <- matrix(nrow=p, ncol=nbf, data=rnorm(p*nbf))

  Psi <- Sdiag - rowSums(B^2)
  crit <- 1
  "++" <- function(x, ...) if (nargs() == 1) x else x + Recall(...)


  itercount <- 0
  
  while ((crit > minerr) & (itercount < maxiter)) {
  
    # E Step:

    tBPsiinv <- t(apply(t(B), 1, function(x) x*(1/Psi)))
    G <- solve(diag(nbf) + tBPsiinv%*%B)

    Z <- G%*%tBPsiinv%*%t(eps)
    S <- lapply(as.data.frame(Z), function(x) G + x%*%t(x))
    names(S) <- NULL

    # M Step:

    Ssum <- do.call("++", S) 

    epsZtlist <- mapply(function(x, y) x%*%t(y), as.data.frame(t(eps)), as.data.frame(Z), SIMPLIFY=FALSE)
    names(epsZtlist) <- NULL
    epsZtsum <- do.call("++", epsZtlist)
    rm(epsZtlist)

    B <- epsZtsum%*%solve(Ssum)
    Psinew <- Sdiag - colSums(t(B)*t(epsZtsum))/n

    if(all(Psinew!=0)) {
      crit <- mean((Psi - Psinew)^2)
      Psi <- Psinew
	}
    else {
      warning(paste("Iteration aborted at crit = ", crit, " > ", minerr, " because there occured zero values in the estimated residual variances.", sep=""))
      crit <- 0
    }	
	
	itercount <- itercount + 1
	
	if(itercount==maxiter)
	  warning(paste("No convergence reached, estimation aborted after maxiter = ", maxiter, " iterations.", sep=""))
	  
  }
  
  # cat(paste("numberiterations: ", itercount, sep=""), "\n")
  
  predfac = function(epsnew)
    t(G%*%tBPsiinv%*%t(epsnew))

  res <- list(B = B, Psi = Psi, Factors=t(Z), predfac=predfac)
  return(res)

}
