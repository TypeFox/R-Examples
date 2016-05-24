sample_mu <-
function (z, X, covlist) {
  nsp <- ncol(z)
  mu <- matrix(NA, nrow(z), nsp)
  blis <- NULL
  for (i in 1:nsp) {
    x <- matrix(X[, covlist[[i]]], ncol = length(covlist[[i]]))
    b <- sample_B(matrix(z[, i], ncol = 1), x)
    mu[, i] <- x %*% b
    blis[[i]] <- b
  }
  list(mu, blis)
}