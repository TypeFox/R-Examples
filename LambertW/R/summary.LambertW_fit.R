#' @rdname LambertW_fit-methods
#' @description \code{summary.LambertW_fit} computes some auxiliary results from
#'     the estimate such as standard errors, theoretical support (only for
#'     \code{type="s"}), skewness tests (only for \code{type="hh"}), etc.  See
#'     \code{print.summary.LambertW_fit} for print out in the console.

#' @return
#' \code{summary} returns a list of class \code{summary.LambertW_fit} 
#' containing 
#' \item{call}{function call} 
#' \item{coefmat}{matrix with 4 columns: \eqn{\widehat{\theta}}, its standard errors, t-statistic, and
#' two-sided p-values} 
#' \item{distname}{see Arguments}
#' \item{n}{number of observations} 
#' \item{data}{original data (\code{y})} 
#' \item{input}{back-transformed input data}
#' \item{support}{support of output random variable Y}
#' \item{data.range}{empirical data range}
#' \item{method}{estimation method} 
#' \item{hessian}{Hessian at the optimum. Numerically obtained for \code{method = "MLE"}; 
#' for \code{method = "IGMM"} a diagonal-matrix approximation from covariance matrix
#' obtained by simulations for \eqn{n = 1000} samples in Goerg (2011).}
#' \item{p_m1, p_m1n}{Probability that one (or n) observation were caused by input 
#' from the non-principal branch (see \code{\link{p_m1}}); only for \code{type = "s"}.}
#' \item{symmetry.p.value}{p-value from Wald test of identical left and right tail parameters (see
#' \code{\link{test_symmetry}}); only for \code{type = "hh"}.}
#' @export
summary.LambertW_fit <- function(object, ...) {
  if (object$method == "IGMM") {
    object$params.hat <- object$tau
    object$distname.tmp <- "normal"
  } else {
    object$distname.tmp <- object$distname
  }
  
  try.inverse <- try(solve(-object$hessian), silent = TRUE)
  
  if (any(class(try.inverse) == "try-error")) {
    se <- rep(NA, length(object$params.hat))
  } else {
    se <- sqrt(diag(try.inverse))
  }
  tval <- object$params.hat/se
  TAB <- cbind(Esimate = object$params.hat, StdErr = se, t.value = tval, 
               p.value = 2 * (1 - pnorm(abs(tval))))
  dimnames(TAB) <- list(names(tval), 
                        c(" Estimate", " Std. Error", " t value", "Pr(>|t|)"))
    
  result <- list(call = object$call,
                 method = object$method,
                 output = data,
                 input = get_input(object$data, tau = object$tau))
  
  dist.family <- get_distname_family(object$distname.tmp)
  support.tmp <- get_support(object$tau, 
                             is.non.negative = dist.family$is.non.negative)    
  result <- c(result,
              list(coefmat = TAB,
                   hessian = object$hessian,
                   distname = object$distname,
                   type = object$type,
                   n = length(object$data),
                   support = support.tmp,
                   data.range = range(object$data),
                   theta = object$theta))
  
  switch(object$method,
         "IGMM" = {
           result$tau <- object$tau
         },
         "MLE" = {
           result$tau <- theta2tau(object$theta, distname = object$distname,
                                   use.mean.variance = object$use.mean.variance)
         })
  
  if (object$method == "MLE" && object$type == "s") {
    result$p_m1 <- p_m1(gamma = object$theta$gamma, distname = object$distname, 
                        beta = object$theta$beta, n = 1,
                        use.mean.variance = object$use.mean.variance)
    result$p_m1n <- p_m1(gamma = object$theta$gamma, distname = object$distname, 
                         beta = object$theta$beta, n = result$n,
                         use.mean.variance = object$use.mean.variance)
  } else {
    result$p_m1 <- result$p_m1n <- NA
  }
  if (object$type == "hh") {
    result$symmetry.p.value <- test_symmetry(object)$p.value
  }
  class(result) <- "summary.LambertW_fit"
  return(result)
} 
