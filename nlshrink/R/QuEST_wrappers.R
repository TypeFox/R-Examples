#' Generate sample eigenvalues from population eigenvalues
#'
#' @param tau (Required) A non-negative numeric vector of population
#'   eigenvalues.
#' @param n (Required) A positive integer representing the number of datapoints
#'   of a hypothetical data matrix with dimension \code{c(n, p = length(tau))}.
#' @return A numeric vector of the same length as \code{tau}, containing the
#'   sample eigenvalue estimates, sorted in ascending order.
#' @description The Marcenko Pastur (MP) law relates the limiting distribution
#'   of the sample eigenvalues to that of the population eigenvalues. In the
#'   finite-dimensional case, the population spectral distribution (PSD) can be
#'   represented as a sum of point masses, and the empirical spectral
#'   distribution (ESD) can be obtained by solving the discretized MP equation.
#'   The QuEST function(see references), uses the quantile function of the ESD
#'   to compute the sample eigenvalues for any given ratio \eqn{c = p/n \in
#'   (0,\infty)}.
#' @references \itemize{ \item Ledoit, O. and Wolf, M. (2015). Spectrum
#'   estimation: a unified framework for covariance matrix estimation and PCA in
#'   large dimensions. Journal of Multivariate Analysis, 139(2) \item Ledoit, O.
#'   and Wolf, M. (2016). Numerical Implementation of the QuEST function.
#'   arXiv:1601.05870 [stat.CO] }
#' @examples
#' lambda_estimate(tau = rep(1,200), n = 300)
#' @export
lambda_estimate <- function(tau, n) {
  if (is.unsorted(tau)) {
    tausort <- sort(tau)
    tauorder <- order(tau)
  } else {
    tausort <- tau
    tauorder <- 1:length(tau)
  }
  Q <- QuEST(tausort, n)
  return (Q$lambda)
}

#' Compute the empirical spectral distribution (ESD) for a set of population
#' eigenvalues
#'
#' @param tau (Required) A non-negative numeric vector of population
#'   eigenvalues.
#' @param n (Required) A positive integer representing the number of datapoints
#'   of a hypothetical data matrix with dimension \code{c(n, p = length(tau))}.
#' @return A named numeric vector of containing points of the ESD. The names
#'   give the corresponding points on the x axis.
#' @description The Marcenko Pastur (MP) law relates the limiting distribution
#'   of the sample eigenvalues to that of the population eigenvalues. In the
#'   finite-dimensional case, the population spectral distribution (PSD) can be
#'   represented as a sum of point masses, and the empirical spectral
#'   distribution (ESD) can be obtained by solving the discretized MP equation.
#'   Theoretical and implementation details in the references.
#' @references \itemize{ \item Ledoit, O. and Wolf, M. (2015). Spectrum
#'   estimation: a unified framework for covariance matrix estimation and PCA in
#'   large dimensions. Journal of Multivariate Analysis, 139(2) \item Ledoit, O.
#'   and Wolf, M. (2016). Numerical Implementation of the QuEST function.
#'   arXiv:1601.05870 [stat.CO] }
#' @examples
#' tau_ESD <- ESD(tau = rep(1,200), n = 300)
#' plot(names(tau_ESD), tau_ESD, ylab="F(x)", xlab="x")
#' @export
ESD <- function(tau, n) {
  Q <- QuEST(tau,n)
  return( setNames(Reduce(c, Q$F), Reduce(c,Q$dis_x_list)) )
}
