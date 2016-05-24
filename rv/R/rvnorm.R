
# ========================================================================
# rvnorm  -  Generate variates from a normal sampling model
# ========================================================================

rvnorm <- function (n=1, mean=0, sd=1, var=NULL, precision) {
  if (!missing(precision)) {
    if (!is.null(dim(precision))) {
      var <- solve(precision) # matrix inverse
    } else {
      sd <- 1/sqrt(precision)
    }
  }
  if (!is.null(var)) {
    if (!is.null(dim(var))) {
       return(.rvmvnorm(n=n, mean=mean, Sigma=var))
    }
    sd <- sqrt(var)
  }
  rvvapply(stats::rnorm, n.=n, mean=mean, sd=sd)
}

.mvrnorm <- function (n = 1, mu, Sigma, tol = 1e-06, empirical = FALSE) { 
  #
  # from the MASS package
  #
    p <- length(mu)
    if (!all(dim(Sigma) == c(p, p)))
        stop("incompatible arguments")
    eS <- eigen(Sigma, symmetric = TRUE)
    if (any(is.na(eS$vectors))) {
      stop("Some eigenvectors of the variance-covariance matrix are missing...")
    }
    ev <- eS$values
    if (!all(ev >= -tol * abs(ev[1])))
        stop("Variance-covariance matrix is not positive definite")
    X <- matrix(rnorm(p * n), n)
    if (empirical) {
        X <- scale(X, TRUE, FALSE)
        X <- X %*% svd(X, nu = 0)$v
        X <- scale(X, FALSE, TRUE)
    }
    X <- mu + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
    nm <- names(mu)
    if (is.null(nm) && !is.null(dn <- dimnames(Sigma)))
        nm <- dn[[1]]
    dimnames(X) <- list(nm, NULL)
    if (n == 1)
        drop(X)
    else t(X)
}

# ========================================================================
# .rvmvnorm  -  Generate multivariate normal random variables.
# ========================================================================
#
# TODO: 
#  - preserve dimensionality when arrays of n,mean,sd are given.
#  - dim(n) , dim(mean), dim(sd), in this order? at least dim(n) 
#    should be the most important.
#  - parameter cumulative = FALSE
#  - what to do with n rv? conditional d'n given n?
#  - what to do with the parameter 'n'? Should repeat, yes, how to implement?
#

.rvmvnorm <- function (n=1, mean, Sigma) {
  if (is.null(dim(Sigma))) {
    Sigma <- diag(Sigma)
  } else if (!is.numeric(Sigma)) 
    stop("Invalid (non-numeric) covariance matrix Sigma")
  else if (nrow(Sigma) != ncol(Sigma)) 
     stop("Invalid (nonsquare) covariance matrix Sigma")
  #
  if (length(mean)==1) {
    mean <- rep(mean, length.out=nrow(Sigma))
  }
  if (length(mean) != nrow(Sigma)) {
    stop("Length of mean vector != nrow of Sigma.")
  }
  #
  nm <- names(mean)
  if (is.rv(n)) stop("n cannot be random (to be implemented)")
  n <- n[1]
  if (is.rv(Sigma)) {
    r <- rvmapply(.mvrnorm, n=n, mu=mean, Sigma=Sigma)
    if (n==1) r <- drop(r)
    if (!is.null(dim(r))) {
      dimnames(r)[[2]] <- nm
    }
    return(r)
  }
  n.sims <- getnsims()
  dim.mean <- dim(mean)
  if (is.list(Sigma)) {
    Sigma <- .expand.as.matrix(Sigma)
  }
  ##if (length(mean) > nrow(Sigma)) {
  ##  # DEBUG: ??
  ##  n.Sigmas <- (length(mean)%/%nrow(Sigma))+(length(mean)%%nrow(Sigma)!=0)
  ##  Sigma <- .expand.as.matrix(rep(rv(Sigma), n.Sigmas))
  ##}
  mean <- rep(mean, length.out = nrow(Sigma))
  mu <- rep(0, length(mean))
  n.all <- n*n.sims
  X <- .mvrnorm(n = n.all, mu = mu, Sigma = Sigma)
  if (n==1) {
    r <- rvsims(X)
  } else {
    dim(X) <- c(n.sims, length(X) %/% n.sims)
    r <- rvsims(X)
    dim(r) <- c(length(mean), n)
  }
  r <- t(r + mean) # Must be in this order to preserve dimensions
  if (n==1) r <- drop(r)
  if (is.null(dim(r))) {
    if (length(r)==length(nm)) names(r) <- nm
  } else {
    dimnames(r) <- list(NULL, nm)
  }
  return(r)
}
