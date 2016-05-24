# Minimize_X { (1/2)||X - A||_F^2 + lam||P*X||_1} s.t. X >= del * I
# ...using ADMM

ProxADMM <- function(A, del, lam, P, rho=.1, tol=1e-6, maxiters=100, verb=FALSE) {
  # Minimize_X { (1/2)||X - A||_F^2 + lam||P*X||_1} s.t. X >= del * I
  #
  # ADMM approach

  # first, check if simple soft-thesholding works... if so, skip the ADMM!
  soft <- SoftThreshold(A, lam * P)
  minev <- min(eigen(soft, symmetric=T, only.values=T)$val)
  if (minev >= del) {
    return(list(X=soft, Z=soft, obj=ComputeProxObjective(soft, A, lam, P)))
  }
  
  p <- nrow(A)
  obj <- NULL
  
  # initialize Z, Y
  Z <- soft
  Y <- matrix(0, p, p)

  # main loop
  for (i in seq(maxiters)) {    
    # update X:
    B <- (A + rho * Z - Y) / (1 + rho)
    if (min(eigen(B, symmetric=T, only.values=T)$val) < del) {
      # note: even though eigen is called twice, only.values=T is
      #       much faster, making this worthwhile.
      eig <- eigen(B, symmetric=T)
      X <- eig$vec %*% diag(pmax(eig$val, del)) %*% t(eig$vec)
    }
    else {
      X <- B
    }
    # check for convergence:
    obj <- c(obj, ComputeProxObjective(X, A, lam, P))
    if (verb)
      cat(" ", obj[i], fill=T)
    if (i > 1)
      if (obj[i] > obj[i - 1] - tol) {
        if (verb)
          cat(" ADMM converged after ", i, " steps.", fill=T)
        break
      }
      
    # update Z:
    Z <- SoftThreshold(X + Y / rho, lam * P / rho)

    # update Y:
    Y <- Y + rho * (X - Z)
  }

  list(X=X, Z=Z, obj=obj)
}

SoftThreshold <- function(x, lam) {
  # note: this works also if lam is a matrix of the same size as x.
  sign(x) * (abs(x) - lam) * (abs(x) > lam)
}

ComputeProxObjective <- function(X, A, lam, P) {
  sum((X-A)^2) / 2 + lam * sum(abs(P*X))
}
