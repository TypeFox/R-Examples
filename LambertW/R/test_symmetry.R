#' @title Test symmetry based on Lambert W heavy tail(s)
#' 
#' @description
#' Performs a test for the null hypothesis of symmetry, \eqn{H_0: \delta_l =
#' \delta_r}, versus the alternative of asymmetry. This can be done using a Wald
#' test of the linear restriction \eqn{H_0: \delta_l - \delta_r = 0} or a
#' likelihood ratio test.
#' 
#' By default it uses \code{"Wald"} test since this only requires the Hessian of
#' the \code{"hh"} Lambert W fit.  The \code{"LR"} test requires the
#' log-likelihood values for both MLEs (type \code{"h"} and \code{"hh"}) and
#' thus takes longer to compute.
#'  
#' @param LambertW.fit an object of class \code{LambertW_fit} with \code{type =
#'     "hh"} or a numeric vector (observed data). If it is data, then an
#'     asymmetric Lambert W \eqn{\times} Gaussian distribution (\code{distname =
#'     "normal"}) with two tail parameters (\code{"hh"}) will be fit to the data
#'     internally and then used as the new \code{LambertW.fit}.
#' @param method test methodology: \code{"Wald"} (default) or a likelihood ratio
#'     \code{"LR"} test
#' @return A list of class \code{"htest"} containing: \item{statistic}{value of
#'     the test statistic,} \item{p.value }{p-value for the test,}
#'     \item{method}{character string describing the test,} \item{data.name}{a
#'     character string giving the name(s) of the data.}
#' @keywords htest
#' @export
#' @examples
#' 
#' 
#' # skewed
#' yy <- rLambertW(n = 500, theta = list(delta = c(0.1, 0.25), beta = c(2, 1)), 
#'                 distname = "normal")
#' fit.ml <- MLE_LambertW(yy, type = "hh", distname = "normal", 
#'                        hessian = TRUE)
#' summary(fit.ml)
#' test_symmetry(fit.ml, "LR")
#' test_symmetry(fit.ml, "Wald")
#' 
#' \dontrun{
#' # symmetric 
#' yy <- rLambertW(n = 500, theta = list(delta = c(0.2, 0.2), beta = c(2, 1)), 
#'                 distname = "normal")
#' fit.ml <- MLE_LambertW(yy, type = "hh", distname = "normal")
#' summary(fit.ml)
#' test_symmetry(fit.ml, "LR")
#' test_symmetry(fit.ml, "Wald")
#' }

test_symmetry <- function(LambertW.fit, method = c("Wald", "LR")) {
  
  method <- match.arg(method)
  if (is.numeric(LambertW.fit)) {
    obj <- MLE_LambertW(LambertW.fit, type = "hh", distname = "normal",
                        hessian = (method == "LR"))
  } else if (inherits(LambertW.fit, "LambertW_fit")) {
    if (LambertW.fit$type != "hh") {
      stop("Estimated LambertW.fit method must be of type 'hh'.")
    }
    obj <- LambertW.fit
  } 
  
  hessian <- obj$hessian
  theta.hat <- obj$params.hat
  KK <- length(theta.hat)
  
  if (method == "Wald") {
    var.theta.hat <- try(-solve(hessian), silent = TRUE)
    if (inherits(var.theta.hat, "try-error")) {
      warning("Hessian was singular or NA. Changed method to 'LR'.")
      method <- "LR"
    }
  }
  
  if (method == "Wald") {
    VV <- var.theta.hat * length(obj$data)
    RR <- matrix(0, ncol = KK, nrow = 1)
    RR[1, names(theta.hat) == c("delta_l", "delta_r")] <- c(1, -1)
    WW <- ((RR %*% VV %*% t(RR)))^(-1) * t(theta.hat) %*% theta.hat
    method <- "Wald test for symmetry (H_0: delta_l - delta_r = 0)"
    pval <- 1 - pchisq(WW, 1)
    statistic <- c(W = WW)
  } else if (method == "LR") {
    method <- "Likelihood ratio test"
    theta.init.sym <- obj$theta
    theta.init.sym$delta <- mean(theta.init.sym$delta)
    mod.h <- MLE_LambertW(obj$data, type = "h", distname = obj$distname,
                          hessian = FALSE,
                          theta.init = theta.init.sym)
    loglik.h <- mod.h$loglik
    loglik.hh <- obj$loglik
    
    lambda <- 2 * (loglik.hh - loglik.h)
    pval <- 1 - pchisq(lambda, 1)
    statistic <- c(D = lambda)
  }
  
  DNAME <- c(deparse(substitute(LambertW.fit)))
  RVAL <- list(statistic = statistic, p.value = pval, method = method, 
               data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
} 
