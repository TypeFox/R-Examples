nearPD <-
function (x,
                    corr = FALSE,
                    keepDiag = FALSE,
                    do2eigen = TRUE,
                    doSym = FALSE,
                    doDykstra = TRUE,
                    only.values = FALSE,
                    only.matrix = TRUE,
                    eig.tol = 1e-06,
                    conv.tol = 1e-07,
                    posd.tol = 1e-08,
                    maxit = 100,
                    trace = FALSE)
{
  n <- ncol(x)
  if (keepDiag)
    diagX0 <- diag(x)
  if (doDykstra)
  {
    D_S <- x
    D_S[] <- 0
  }
  X <- x
  iter <- 0
  converged <- FALSE
  conv <- Inf
  while (iter < maxit && !converged)
  {
    Y <- X
    if (doDykstra)
      R <- Y - D_S
    e <- eigen(
      if (doDykstra)
        R
      else
        Y, symmetric = TRUE)
    Q <- e$vectors
    d <- e$values
    p <- d > eig.tol * d[1]
    if (!any(p))
      stop("Matrix seems negative semi-definite")
    Q <- Q[, p, drop = FALSE]
    X <- tcrossprod(Q * rep(d[p], each = nrow(Q)), Q)
    if (doDykstra)
      D_S <- X - R
    if (doSym)
      X <- (X + t(X))/2
    if (corr)
      diag(X) <- 1
    else
      if (keepDiag)
        diag(X) <- diagX0
    conv <- norm(Y - X, "I")/norm(Y, "I")
    iter <- iter + 1
    if (trace)
      cat(sprintf("iter %3d : #{p}=%d, ||Y-X|| / ||Y||= %11g\n",
                  iter, sum(p), conv))
    converged <- (conv <= conv.tol)
  }
  if (!converged)
  {
    warning("nearPD() did not converge in ", iter, " iterations")
  }
  if (do2eigen || only.values)
  {
    e <- eigen(X, symmetric = TRUE)
    d <- e$values
    Eps <- posd.tol * abs(d[1])
    if (d[n] < Eps)
    {
      d[d < Eps] <- Eps
      if (!only.values)
      {
        Q <- e$vectors
        o.diag <- diag(X)
        X <- Q %*% (d * t(Q))
        D <- sqrt(pmax(Eps, o.diag)/diag(X))
        X[] <- D * X * rep(D, each = n)
      }
    }
    if (only.values)
      return(d)
    if (only.matrix)
      return(matrix(X, ncol = n, nrow = n))
    if (corr)
      diag(X) <- 1
    else
      if (keepDiag)
        diag(X) <- diagX0
  }
  structure(list(mat = matrix(data = X, nrow = n, ncol = n,
                              dimnames = dimnames(x)), eigenvalues = d, corr = corr,
                 normF = norm(x - X, "F"), iterations = iter, rel.tol = conv,
                 converged = converged), class = "nearPD")
}
