## binary.R -- Part of the bayesGDS package 
## Copyright (C) 2013-2015 Michael Braun
## See LICENSE file for details.

## Functions to compute objective function, gradient, and Hessian for
## binary choice example, and some unit tests.


#' @name binary
#' @title Binary choice example
#' @description Functions for binary choice example in the vignette.
#' @param P Numeric vector of length (N+1)*k.  First N*k elements are heterogeneous coefficients. The remaining k elements are population parameters.
#' @param data List of data matrices Y and X, and choice count integer T
#' @param priors List of named matrices inv.Omega and inv.Sigma
#' @return Log posterior density, gradient and Hessian.
#' @details Hessian is sparse, and returned as a dgcMatrix object
#' @rdname binary
#' @export
binary.f <- function(P, data, priors) {

    N <- length(data$Y)
    k <- NROW(data$X)
    beta <- matrix(P[1:(N*k)], k, N)
    mu <- P[(N*k+1):(N*k+k)]
    
    bx <- colSums(data$X * beta)
    
    log.p <- bx - log1p(exp(bx))
    log.p1 <- -log1p(exp(bx))
    
    data.LL <- sum(data$Y*log.p + (data$T-data$Y)*log.p1)
    
    Bmu <- apply(beta, 2, "-", mu)
    
    prior <- -0.5 * sum(diag(tcrossprod(Bmu) %*% priors$inv.Sigma))
    hyp <- -0.5 * t(mu) %*% priors$inv.Omega %*% mu
    res <- data.LL + prior + hyp
    return(as.numeric(res))
    
}

#' @rdname binary
#' @export
binary.grad <- function(P, data, priors) {

    Y <- data$Y
    X <- data$X
    T <- data$T
    inv.Sigma <- priors$inv.Sigma
    inv.Omega <- priors$inv.Omega
    
  q1 <- .dlog.f.db(P, Y, X, T, inv.Omega, inv.Sigma)
  q2 <- .dlog.f.dmu(P, Y, X, T, inv.Omega, inv.Sigma)
  res <- c(q1, q2)
  return(res)
}

#' @rdname binary
#' @export
binary.hess <- function(P, data, priors) {

    Y <- data$Y
    X <- data$X
    T <- data$T
    inv.Sigma <- priors$inv.Sigma
    inv.Omega <- priors$inv.Omega
    N <- length(Y)
    k <- NROW(X)
        
    SX <- Matrix(inv.Sigma)
    XO <- Matrix(inv.Omega)
    B1 <- .d2.db(P, Y, X, T, SX)
    
    cross <- .d2.cross(N, SX)    
    Bmu <- .d2.dmu(N,SX, XO)
    res <- rBind(cBind(B1, Matrix::t(cross)),cBind(cross, Bmu))
    
    return(res)
}



.dlog.f.db <- function(P, Y, X, T, inv.Omega, inv.Sigma) {

    N <- length(Y)
    k <- NROW(X)
    
    
    beta <- matrix(P[1:(N*k)], k, N)
    mu <- P[(N*k+1):length(P)]
    bx <- colSums(X * beta)
    
    p <- exp(bx)/(1+exp(bx))
    
    tmp <- Y - T*p
    
    dLL.db <- apply(X,1,"*",tmp)
    
    Bmu <- apply(beta, 2, "-", mu)
    dprior <- -inv.Sigma %*% Bmu
    
    res <- t(dLL.db) + dprior
    
    return(as.vector(res))
    
}

.dlog.f.dmu <- function(P, Y, X, T, inv.Omega, inv.Sigma) {

    N <- length(Y)
    k <- NROW(X)
    
    beta <- matrix(P[1:(N*k)], k, N)
    mu <- P[(N*k+1):length(P)]
    Bmu <- apply(beta, 2, "-", mu)
    
    res <- inv.Sigma %*% (rowSums(Bmu)) -  inv.Omega %*% mu
    return(res)
}



.d2.db <- function(P, Y, X, T, inv.Sigma) {

    N <- length(Y)
    k <- NROW(X)
    
    beta <- matrix(P[1:(N*k)], k, N)
    mu <- P[(N*k+1):length(P)]
    ebx <- exp(colSums(X * beta))
    
    p <- ebx/(1+ebx)
    
    q <- vector("list",length=N)
    for (i in 1:N) {
        q[[i]] <- -T*p[i]*(1-p[i])*tcrossprod(X[,i]) - inv.Sigma
    }
    B <- bdiag(q)
    return(B)    
}


.d2.dmu <- function(N, inv.Sigma, inv.Omega) {
  return(-N*inv.Sigma-inv.Omega)
}

.d2.cross <- function(N, inv.Sigma) {
  res <- kronecker(Matrix(rep(1,N),nrow=1),inv.Sigma)
  return(res)
}




