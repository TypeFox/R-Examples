#' Multivariable Random Uniform data
#'
#' Generate a simulated multivariable random uniformly distributed dataset
#' using the method of Cholesky Decomposition.
#'
#' @importFrom stats glm.fit binomial runif
#'
#' @param n   the number of rows of observations in the dataset
#' @param mu  a vector containing the column means of the dataset
#' @param Cov   a covariance matrix
#' @param Q   an optional orthogonal matrix
#' @return A simulated matrix of values based on the input parameters is returned.
#'
#' @examples
#'## Simulated data based on the iris dataset
#' mu <- c(rep(0, 4))
#' covmatr <- matrix(c(0.7, -0.04, 1.3, 0.5, -0.04, 0.2, -0.3, -0.1,
#' 1.3, -0.3, 3.1, 1.3, 0.5, -0.1, 1.3, 0.6), ncol = 4)
#' sim.dat <- randunif(n = 100, mu = mu, Cov = covmatr)
#' head(sim.dat)
#'
#' @references Rizzo M. L., \emph{"Statistical Computing with R", Chapman & Hall/CRC} (2007)

randunif<-function(n, mu, Cov, Q) {

  if (missing(Q)) Q <- diag(1, dim(Cov)[1])

  m <- dim(Cov)[1]
  mu <- matrix(rep(1, n), ncol = 1) %*% mu
  y <- runif(n * m, 0, 1)
  y <- 2 * sqrt(3) * (y - 0.5)
  y <- matrix(y, ncol = m)
  DQ <- eigen(Cov)
  QC <- DQ$vectors
  D <- diag(sqrt(DQ$values))
  E <- QC %*% D %*% t(QC)
  E <- Re(E)
  y <- y %*% Q
  unifdat <- y %*% E + mu
  return(unifdat)
}
