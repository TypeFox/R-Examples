# Method gJADE
gJADE <- function(X,...) UseMethod("gJADE")

# main function for gJADE
gJADE.default <- function(X, k = 0:12, eps = 1e-06, maxiter = 100, method = "frjd", ...)
{
    nk <- length(k)
    method <- match.arg(method, c("rjd", "djd", "frjd"))
    MEAN <- colMeans(X)
    COV <- cov(X)
    p <- ncol(X)
    EVD <- eigen(COV, symmetric = TRUE)
    COV.sqrt.i <- EVD$vectors %*% tcrossprod(diag(EVD$values^(-0.5)), EVD$vectors)
    X.C <- sweep(X, 2, MEAN, "-")
    Y <- tcrossprod(X.C, COV.sqrt.i)
    
    ccks <- CCK(Y, k)
    U <- switch(method, frjd = {
      frjd(ccks, eps = eps, maxiter = maxiter, ...)$V
    }, rjd = {
      rjd(ccks, eps = eps, maxiter = maxiter, ...)$V
    }, djd = {
      djd(ccks, eps = eps, maxiter = maxiter, ...)
    })
    W <- crossprod(U, COV.sqrt.i)
    S <- tcrossprod(X.C, W)
    S <- ts(S, names = paste("Series", 1:p))
    RES <- list(W = W, k = k, S = S)
    class(RES) <- "bss"
    RES
}


gJADE.ts <- function(X, ...)
{
  x <- as.matrix(X)
  RES <- gJADE.default(x,...)
  S <- RES$S
  attr(S, "tsp") <- attr(X, "tsp")
  RES$S <- S
  RES
}
