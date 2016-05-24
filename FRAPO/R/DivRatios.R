##
## Diversification Ratio
##
dr <- function(weights, Sigma){
  Sigma <- as.matrix(Sigma)
  if(!isSymmetric(Sigma)){
    stop("Object provided for 'Sigma' is not a symmetric matrix.\n")
  }
  w <- as.vector(weights)
  if(length(w) != ncol(Sigma)){
    stop("Length of 'weights' vector differs from row/column dimension of 'Sigma'.\n")
  }
  nom <- w %*% sqrt(diag(Sigma))
  denom <- sqrt(t(w) %*% Sigma %*% w)
  res <- as.numeric(nom / denom)
  return(res)
}
##
## Concentration Ratio
##
cr <- function (weights, Sigma){
  Sigma <- as.matrix(Sigma)
  if (!isSymmetric(Sigma)) {
    stop("Object provided for 'Sigma' is not a symmetric matrix.\n")
  }
  w <- as.vector(weights)
  if (length(w) != ncol(Sigma)) {
    stop("Length of 'weights' vector differs from row/column dimension of 'Sigma'.\n")
  }
  prod <- weights * sqrt(diag(Sigma))
  nom <- sum(prod^2)
  denom <- sum(prod)^2
  res <- as.numeric(nom / denom)
  return(res)
}
##
## Volatility weighted correlation
##
rhow <- function(weights, Sigma){
  Sigma <- as.matrix(Sigma)
  if (!isSymmetric(Sigma)) {
    stop("Object provided for 'Sigma' is not a symmetric matrix.\n")
  }
  w <- as.vector(weights)
  if (length(w) != ncol(Sigma)) {
    stop("Length of 'weights' vector differs from row/column dimension of 'Sigma'.\n")
  }
  idx <- which(upper.tri(Sigma), arr.ind = TRUE)
  S <- sqrt(diag(Sigma))
  C <- cov2cor(Sigma)
  prod <- w[idx[, 1]] * S[idx[, 1]] * w[idx[, 2]] * S[idx[, 2]]
  nom <- sum(prod * C[idx])
  denom <- sum(prod)
  res <- as.numeric(nom / denom)
  return(res)
}
