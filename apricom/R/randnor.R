#' Multivariable Random Normal data
#'
#' Generate a simulated multivariable random normally distributed dataset
#' using the method of Cholesky Decomposition.
#'
#' @importFrom stats glm.fit binomial rnorm
#'
#' @param n   the number of rows of observations in the dataset
#' @param mu  a vector of length m containing the column means of the dataset
#' @param Cov   an m x m covariance matrix
#' @return A simulated matrix of values based on the input parameters is returned.
#' @examples
#'## Simulated data based on the iris dataset
#' mu <- c(rep(0, 4))
#' covmatr <- matrix(c(0.7, -0.04, 1.3, 0.5, -0.04, 0.2, -0.3, -0.1,
#' 1.3, -0.3, 3.1, 1.3, 0.5, -0.1, 1.3, 0.6), ncol = 4)
#' sim.dat <- randnor(n = 100, mu = mu, Cov = covmatr)
#' head(sim.dat)
#'
#' @references Rizzo M. L., \emph{"Statistical Computing with R", Chapman & Hall/CRC} (2007)

randnor <- function(n, mu, Cov){

m <- dim(Cov)[1]
mu <- matrix(rep(1, n), ncol = 1) %*% mu
y <- rnorm(n * m, 0, 1)
y <- matrix(y, ncol = m)
DQ <- eigen(Cov)
Q <- DQ$vectors
D <- diag(sqrt(DQ$values))
E <- Q %*% D %*% t(Q)
E <- Re(E)
normdat <- y %*% E + mu
return(normdat)
}
