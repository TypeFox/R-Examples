#' Package
#'
#' @description A package for estimating population eigenvalues and covariance
#'   matrices, based on publications by Ledoit and Wolf (2004, 2012, 2015,
#'   2016).
#'
#' @details A common assumption in statistics is that for a data matrix \eqn{X}
#'   of dimension \eqn{n \times p}, the number of predictor variables (p)
#'   vanishes relative to the number of datapoints (n) as \eqn{n \to \infty}.
#'   However, in modern datasets, it is often the case that p is comparable to
#'   or greater than n. In this scenario, a more appropriate asymptotic
#'   framework is to assume that the ratio \eqn{c := p/n} approaches a finite
#'   positive value as \eqn{n,p \to \infty}. In this case, the sample covariance
#'   matrix \eqn{S} is no longer a consistent estimator of the population
#'   covariance matrix \eqn{\Sigma}. Similarly, the sample eigenvalues deviate
#'   substantially from the population eigenvalues. This package contains
#'   implementations of Ledoit and Wolf's linear and non-linear shrinkage
#'   population eigenvalue and covariance estimation methods, based on their
#'   2016 publication and the accompanying MATLAB code. Theoretical and
#'   implementation details of these methods can be found in the following
#'   publications:
#'
#'   \itemize{ \item Ledoit, O. and Wolf, M. (2004). A well-conditioned
#'   estimator for large-dimensional covariance matrices. Journal of
#'   Multivariate Analysis, 88(2) \item Ledoit, O. and Wolf, M. (2012).
#'   Nonlinear shrinkage estimation of large-dimensional covariance matrices.
#'   Annals of Statistics, 40(2). \item Ledoit, O. and Wolf, M. (2015). Spectrum
#'   estimation: a unified framework for covariance matrix estimation and PCA in
#'   large dimensions. Journal of Multivariate Analysis, 139(2). \item Ledoit,
#'   O. and Wolf, M. (2016). Numerical Implementation of the QuEST function.
#'   arXiv:1601.05870 [stat.CO]. }
#'
#' @docType package
#' @name nlshrink
#'
#' @import stats
#' @importFrom MASS mvrnorm
#' @importFrom graphics legend
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom nloptr nloptr
NULL
