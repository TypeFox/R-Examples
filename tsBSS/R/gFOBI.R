# Method gFOBI
gFOBI <- function(X,...) UseMethod("gFOBI")

# main function for gFOBI
gFOBI.default <- function (X, k = 0:12, eps = 1e-06, maxiter = 100, method = "frjd", ...) 
{
  nk <- length(k)
  method <- match.arg(method, c("rjd", "djd", "frjd"))
  MEAN <- colMeans(X)
  COV <- cov(X)
  EVD <- eigen(COV, symmetric = TRUE)
  COV.sqrt.i <- EVD$vectors %*% tcrossprod(diag(EVD$values^(-0.5)), EVD$vectors)
  X.C <- sweep(X, 2, MEAN, "-") #Centered
  Y <- tcrossprod(X.C, COV.sqrt.i) #Scaled
  p <- ncol(X)
  R <- array(0, dim = c(p, p, nk))
  n <- nrow(X)
  for (i in 1:nk) {
    Yt <- Y[1:(n - k[i]), ]
    Yti <- Y[(1 + k[i]):n, ]
    r <- sqrt(rowSums(Yt^2))
    Yu <- r*Yti
    Ri <- crossprod(Yu)/nrow(Yt)
    R[, , i] <- Ri
  }
  JD <- switch(method, frjd = {
    frjd(R, eps = eps, maxiter = maxiter, ...)$V
  }, rjd = {
    rjd(R, eps = eps, maxiter = maxiter, ...)$V
  }, djd = {
    djd(R, eps = eps, maxiter = maxiter, ...)
  })
  W <- crossprod(JD, COV.sqrt.i)
  W <- diag(sign(rowMeans(W))) %*% W
  S <- tcrossprod(X.C, W)
  S <- ts(S, names = paste("Series", 1:p))
  RES <- list(W = W, k = k, S = S)
  class(RES) <- "bss"
  RES
}

gFOBI.ts <- function(X, ...)
{
  x <- as.matrix(X)
  RES <- gFOBI.default(x,...)
  S <- RES$S
  attr(S, "tsp") <- attr(X, "tsp")
  RES$S <- S
  RES
}

