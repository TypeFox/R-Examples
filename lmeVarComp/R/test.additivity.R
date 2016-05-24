test.additivity <-
function(x, y, type = c("RLR", "pseudo", "score"),
  nbasis = 10L, kernel = c("gaussian", "polynomial", "spline"),
  nsim = 5000L, seed = 130623L)
{
  type <- match.arg(type)
  kernel <- match.arg(kernel)
  stopifnot((min(x) >= 0) && (max(x) <= 1))
  
  n <- nrow(x)
  p <- ncol(x)
  Y <- y
  
  # linear effects
  X <- cbind(1, x)
  
  # additive effects
  cubic <- additive.cubic.spline(x, nbasis)
  Z1 <- cubic$Z
  V1 <- cbind(X, Z1)
  Sigma1 <- sinv(cubic$S)

  # joint effects with orthogonality constraint imposed
  R <- pairwise.product(x)
  Z2 <- R - V1 %*% mnls(V1, R)
  V1 <- cbind(V1, R)
  Sigma2 <- diag(ncol(Z2))
  
  R <- switch(kernel,
    "gaussian" = gaussian.kernel(x),
    "polynomial" = polynomial.kernel(x),
    "spline" = spline.kernel(x))
  Z3 <- qr.Q(qr(V1, LAPACK = TRUE), 
    complete = TRUE)[, (ncol(V1) + 1L) : n, drop = FALSE]
  Sigma3 <- crossprod(Z3, R %*% Z3)
  Sigma3 <- (Sigma3 + t(Sigma3)) / 2

  # tests zero variance component
  Z <- list(Z1, Z2, Z3)
  Sigma <- list(Sigma1, Sigma2, Sigma3)
  if (type == "RLR") {
    result <- rlr.test(Y, X, Z, Sigma, 1L, nsim, seed)
    result <- c(result$RLRT, result$GFT)
    names(result) <- c("RLRT.stat.obs", "RLRT.p.value", 
      "GFT.stat.obs", "GFT.p.value")
  } else if (type == "pseudo") {
    result <- pseudo.rlr.test(Y, X, Z, Sigma, 1L, nsim, seed)
  } else if (type == "score") {
    result <- score.test(Y, X, Z, Sigma, 1L)
  }

  # returns results
  result
}
