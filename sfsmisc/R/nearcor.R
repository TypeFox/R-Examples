#### Copyright (2007) Jens Oehlschlägel
#### GPL licence, no warranty, use at your own risk

### NOTA BENE: nearPD() in package Matrix is a new version, slightly more elegant
###            ^^^^^^^^ also using Matrix-builtin functionality

nearcor <- function(  # Computes the nearest correlation matrix to an approximate correlation matrix, i.e. not positive semidefinite.
  R                   # n-by-n approx correlation matrix
, eig.tol   = 1.0e-6  # defines relative positiveness of eigenvalues compared to largest
, conv.tol  = 1.0e-7  # convergence tolerance for algorithm
, posd.tol  = 1.0e-8  # tolerance for enforcing positive definiteness
, maxits    = 100     # maximum number of iterations allowed
, verbose   = FALSE   # set to TRUE to verbose convergence

                      # RETURNS list of class nearcor with components cor, iterations, converged
){
  if (!(is.numeric(R) && is.matrix(R) && identical(R,t(R))))
    stop('Error: Input matrix R must be square and symmetric')

  # Inf norm
  inorm <- function(x)max(rowSums(abs(x)))
  # Froebenius norm
  fnorm <- function(x)sqrt(sum(diag(t(x) %*% x)))

  n <- ncol(R)
  U <- matrix(0, n, n)
  Y <- R
  iter <- 0

  while (TRUE){
      T <- Y - U

      # PROJECT ONTO PSD MATRICES
      e <- eigen(Y, symmetric=TRUE)
      Q <- e$vectors
      d <- e$values
      D <- diag(d)

      # create mask from relative positive eigenvalues
      p <- (d>eig.tol*d[1]);

      # use p mask to only compute 'positive' part
      X <- Q[,p,drop=FALSE] %*% D[p,p,drop=FALSE] %*% t(Q[,p,drop=FALSE])

      # UPDATE DYKSTRA'S CORRECTION
      U <- X - T

      # PROJECT ONTO UNIT DIAG MATRICES
      X <- (X + t(X))/2
      diag(X) <- 1

      conv <- inorm(Y-X) / inorm(Y)
      iter <- iter + 1
      if (verbose)
        cat("iter=", iter, "  conv=", conv, "\n", sep="")

      if (conv <= conv.tol){
        converged <- TRUE
        break
      }else if (iter==maxits){
        warning(paste("nearcor did not converge in", iter, "iterations"))
        converged <- FALSE
        break
      }
      Y <- X
  }
  X <- (X + t(X))/2
  # begin from posdefify(sfsmisc)
  e <- eigen(X, symmetric = TRUE)
  d <- e$values
  Eps <- posd.tol * abs(d[1])
  if (d[n] < Eps) {
      d[d < Eps] <- Eps
      Q <- e$vectors
      o.diag <- diag(X)
      X <- Q %*% (d * t(Q))
      D <- sqrt(pmax(Eps, o.diag)/diag(X))
      X[] <- D * X * rep(D, each = n)
      ## force symmetry
      X <- (X + t(X))/2
  }
  # end from posdefify(sfsmisc)
  diag(X) <- 1
  ret <- list(cor=X, fnorm=fnorm(R-X), iterations=iter, converged=converged)
  class(ret) <- "nearcor"
  ret
}
