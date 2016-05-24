#' Simulate from an LSBLCUST model
#' 
#' This is a simple functiomn to simulate data from an LSBCLUST model.
#' 
#' @param N   Integer giving the number of observations required.
#' @param nclust Integer giving the number of clusters required.
#' @param J Integer giving the number of rows to sample.
#' @param K Integer giving the number of columns ot sample.
#' @param ndim Integer giving the desired dimensionality of the sampled cluster means.
#' @param sd The sd argument of \code{\link{rnorm}}.
#' @export
sim.lsbclust <- function(N = 100, nclust = 5, J = 10, K = 8, ndim = 2, sd = 1) {
  Cs <- replicate(nclust, matrix(rnorm(J*ndim), J, ndim), simplify = FALSE)
  Ds <- replicate(nclust, matrix(rnorm(K*ndim), K, ndim), simplify = FALSE)
  CDs <- Map(tcrossprod, Cs, Ds)
  data <- array(rnorm(J*K*N, sd = sd), dim = c(J, K, N))
  cluster <- sample.int(n = nclust, size = N, replace = TRUE)
  for (i in 1:N) data[, , i] <- data[, , i] + CDs[[cluster[i]]]
  return(list(data = data, cluster = cluster, Cs = Cs, Ds = Ds, CDs = CDs))
}