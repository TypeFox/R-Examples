graf.fit.laplace <-
  function (y, x, mn, l, wt, e = NULL, tol  = 10 ^ -6, itmax = 50,
            verbose = FALSE) {
    
    if (is.vector(x)) x <- as.matrix(x)
    mn <- qnorm(mn)
    n <- length(y)
    
    # create the covariance matrix
    K <- cov.SE(x1 = x, e1 = e, e2 = NULL, l = l)
    
    # an identity matrix for the calculations
    eye <- diag(n) 
    # initialise
    a <- rep(0, n)
    f <- mn
    obj.old <- Inf
    obj <- -sum(wt * d0(f, y))
    it <- 0
    
    # start newton iterations
    while (obj.old - obj > tol & it < itmax) {
      it <- it + 1
      obj.old <- obj
      W <- -(wt * d2(f, y))
      rW <- sqrt(W)
      cf <- f - mn
      mat1 <- rW %*% t(rW) * K + eye
      L <- tryCatch(chol(mat1),
                    error = function(x) return(NULL))
      b <- W * cf + wt * d1(f, y)
      mat2 <- rW * (K %*% b)
      adiff <- b - rW * backsolve(L, forwardsolve(t(L), mat2)) - a 
      dim(adiff) <- NULL
      
      # find optimum step size using Brent's method
      res <- optimise(psiline, c(0, 2), adiff, a, as.matrix(K), y, d0, mn, wt)
      a <- a + res$minimum * adiff
      f <- K %*% a + mn
      obj <- psi(a, f, mn, y, d0, wt)
      
    }
    
    # return marginal negative log-likelihood
    lp <- sum(wt * d0(f, y))
    mnll <- (a %*% cf)[1, 1] / 2 - lp + sum(log(diag(L)))
    
    if(verbose ) cat(paste("  ", it, "Laplace iterations\n"))
    if(it == itmax) print("timed out, don't trust the inference!")
    return(list(y = y, x = x, MAP = f, ls = l, a = a, W = W, L = L, K = K,
                e = e, obsx = x, obsy = y, mnll = mnll, wt = wt))
  }