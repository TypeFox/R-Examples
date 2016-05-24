### Computes Tukey's biweight rho function with constant `c` for all
### values in the vector `x`.
rho.biweight <- function(x, c) {
  hulp <- x^2/2 - x^4/(2 * c^2) + x^6/(6 * c^4)
  ifelse(abs(x) < c, hulp, c^2/6)
}


## Computes the first derivative of Tukey's biweight psi function with
## constant `c` for all values in vector `x`.
psi.prime <- function(x, c) {
  psi.bisquare(x, c, deriv=TRUE)
}


IF.obs <- function(x, B, const, c) {
  p <- length(x)
  norm.x <- sqrt(sum(x^2))
  ifs0 <- ((x %*% t(x)/norm.x^2) -
           diag(p)/p) * psi.bisquare(norm.x, c) * norm.x^2 * p/const$gamma1 +
           (2/const$gamma3) * (rho.biweight(norm.x, c) - const$b0) * diag(p)
  ifs <- B %*% ifs0 %*% t(B)
  ifs
}


## Performs matrix multiplication of each slices the 'x' array
## along the Z-axis and B
multarray <- function(x, B) {
  result <- apply(x, 3, '%*%', B)
  dim(result) <- c(dim(x)[1], ncol(B), dim(x)[3])
  result
}


## Performs matrix multiplication of B and each slices the 'x' array
## along the Z-axis
tmultarray <- function(x, B) {
  result <- apply(x, 3, function(x.slice) B %*% x.slice)
  dim(result) <- c(nrow(B), dim(x)[2], dim(x)[3])
  result
}

