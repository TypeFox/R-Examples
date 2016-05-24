#' Generate PCA data and Calculates Correlation Matrices
#' 
#' Generate a response style data set from a specific correlation matrix, clean
#' the data with constrained dual scaling and report the original, cleaned and
#' contaminated correlation matrices in a list.
#' 
#' @param nr.indv Vector; number of individuals in each response style group.
#' It is passed to \code{\link{simpca}}.
#' @param m scalar; Number of items.
#' @param q scalar; Number of rating categories, such that the rating scale is
#' \code{1:q}.
#' @param r scalar; Rank of simulated correlation matrices.
#' @param err.coeff scalar; Standard deviation used in simulations that is
#' passed on to \code{\link{simpca}}.
#' @param alphamat matrix; Contains the spline parameters for the different
#' response styles that is passed to \code{\link{simpca}}.
#' @param randomize logical; See \code{\link{simpca}}.
#' @param ...  Arguments passed to \code{\link{cds}}.
#' @return A list with components: \item{Rsim}{Correlation matrix from which
#' the sample was generated} \item{Rclean}{Correlation matrix for the cleaned
#' data} \item{Rcont}{Correlation matrix for the contaminated data}
#' @author Pieter C. Schoonees
#' @export genPCA
genPCA <- function(nr.indv = rep(100,5), m = 10, q = 7, r = 3, err.coeff = 0.1, 
                   alphamat = rbind(c(0.5, 2, 4), c(10, 2, 10), c(1, 2, 1), 
                                    c(4, 2, 0.5), c(0.1, 2, 0.1))[1:length(nr.indv),], randomize = TRUE, ...){
  K <- length(nr.indv)
  L <- matrix(rnorm(m*r), nrow = m, ncol = r)
  L <- sweep(L, 1, sqrt(rowSums(L^2)), "/")
  R <- tcrossprod(L)
  dat <- simpca(nr.indv = nr.indv, m = m, q = q, R = list(R = R, L = L), err.coeff = err.coeff, 
                alphamat = alphamat, randomize = randomize)
  out <- cds(dat, K = K, parallel = TRUE, ...)
  clean <- clean.scales(data = dat, ds.obj = out)
  Rclean <- cor(clean$clean)
  Rcont <- cor(dat$data)
  list(Rsim = R, Rclean = Rclean, Rcont = Rcont)
}
