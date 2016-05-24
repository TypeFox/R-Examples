#' Demonstration of non-linear shrinkage estimator of population eigenvalues
#'
#' @param tau (Optional) Input population eigenvalues. Non-negative Numeric
#'   vector of length \code{p}.
#' @param n (Optional) Number of rows in simulated data matrix (Default 300).
#' @param p (Optional) Number of columns in simulated data matrix (Default 300).
#' @param method (Optional) The optimization routine called in
#'   \code{\link{tau_estimate}}. Choices are \code{nlminb} (default) and
#'   \code{nloptr}.
#' @param control (Optional) A list of control parameters. Must correspond to
#'   the selected optimization method. See \code{\link[stats]{nlminb}},
#'   \code{\link[nloptr]{nloptr}} for details.
#' @return NULL
#' @description This is a demonstration of the non-linear shrinkage method for
#'   estimating population eigenvalues. The inputted population eigenvalues are
#'   used to simulate a matrix of multivariate normal random variates with
#'   covariance matrix \code{diag(tau)}. This data matrix is then used to
#'   estimate the input population eigenvalues using the non-linear shrinkage
#'   method. The output plot shows the comparison of the various estimators
#'   (sample eigenvalues, linear shrinkage, non-linear shrinkage) to the true
#'   population eigenvalues.
#' @section NOTE: \code{nlminb} is usually robust and accurate, but does not
#'   allow equality constraints, so the sum of the estimated population
#'   eigenvalues is in general not equal to the sum of the sample eigenvalues.
#'   \code{nloptr} enforces an equality constraint to preserve the trace, but is
#'   substantially slower than \code{nlminb}. The default optimizer used for
#'   \code{nloptr} is the Augmented Lagrangian method with local optimization
#'   using LBFGS. These can be modified using the control parameter.
#'   \code{Rcgmin} does not enforce equality constraints, but may be more
#'   efficient for certain higher dimensional problems. The ideal optimization
#'   routine depends on the underlying structure of the population eigenvalues.
#' @references \itemize{\item Ledoit, O. and Wolf, M. (2015). Spectrum
#'   estimation: a unified framework for covariance matrix estimation and PCA in
#'   large dimensions. Journal of Multivariate Analysis, 139(2) \item Ledoit, O.
#'   and Wolf, M. (2016). Numerical Implementation of the QuEST function.
#'   arXiv:1601.05870 [stat.CO]}
#' @export
nlshrink_demo <- function(tau = NULL, n = 300, p = 300, method="nlminb", control = list()) {
  if (is.null(tau)) {
    z = seq(0,1,length.out=p)
    z1 <- z[z <= 1/2]
    expo = 3
    tau1=(1-(1-(2*z1)^expo)^(1/expo))/2;
    tau <- 1 + 9*c(tau1, (1-rev(tau1)))
    tauorder <- 1:p
  } else {
    if (length(tau) != p) stop("tau must be of length p")
    if (is.unsorted(tau)) {
      tauorder <- order(tau)
      tau <- sort(tau)
    } else tauorder <- 1:p
  }
  Sigma <- diag(tau)
  X <- mvrnorm(n, mu = rep(0,p), Sigma = Sigma)
  S <- 1/(n-1)*crossprod(X)
  lambda <- sort(eigen(S,only.values=TRUE)$values)
  tau_est <- tau_estimate(X, method=method, control = control)
  lambda_ls <- linshrink(X)
  plot(lambda[tauorder], col="red", ylab="Eigenvalues", cex=0.5)
  points(lambda_ls[tauorder], col = "green", cex=0.5)
  points(tau_est[tauorder], pch=16, cex=0.5)
  points(tau[tauorder], col="blue", cex=0.5)
  legend("topleft", legend=c("Sample", "Population", "Non-linear shrinkage estimate", "Linear shrinkage estimate"), pch=c(1,1,16), col=c("red","blue","black","green"))
}
