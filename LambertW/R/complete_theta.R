#' @rdname theta-utils
#' @description
#' \code{complete_theta} completes missing values in a parameters list so users
#' don't have to specify everything in detail. If not supplied, then
#' \code{alpha = 1}, \code{gamma = 0}, and \code{delta = 0} will be set by default.
#' @param LambertW.input optional; if \code{beta} is missing in \code{theta},
#' \code{LambertW.input} (which has a \code{beta} element) must be specified.
#' 
#' @return
#' \code{complete_theta} returns a list containing:
#' \item{alpha}{ heavy tail exponent(s),}
#' \item{beta}{ named vector \eqn{\boldsymbol \beta} of the input distribution, }
#' \item{gamma}{ skewness parameter,} 
#' \item{delta}{ heavy-tail parameter(s).}
#' 
#' @export
complete_theta <- function(theta = list(), LambertW.input = NULL) {

  if (is.null(theta$beta)) {
    if (is.null(LambertW.input$beta)) {
      stop("If you don't specify a 'beta' in 'theta', you must pass 'LambertW.input'.")
    }
    theta$beta <- LambertW.input$beta
  }
  
  # type 's'
  if (is.null(theta$gamma)) {
    theta$gamma <- 0
  } else {
    names(theta$gamma) <- NULL
  }  
  
  if (is.null(theta$alpha)) {
    theta$alpha <- 1
  }
  if (is.null(theta$delta)) {
    theta$delta <- 0
  }
  num.alphas <- length(theta$alpha)
  num.deltas <- length(theta$delta)

  stopifnot(num.deltas <= 2,
            num.alphas <= 2)
  
  # match alpha and delta to same length
  if (num.deltas == 2 && num.alphas == 1) {
    theta$alpha <- rep(theta$alpha, 2)
  }
  if (num.alphas == 2 && num.deltas == 1) {
    theta$delta <- rep(theta$delta, 2)
  }
  
  if (length(theta$alpha) == 1) {
    names(theta$alpha) <- NULL
  } else {
    names(theta$alpha) <- c("alpha_l", "alpha_r")
  }
  
  if (length(theta$delta) == 1) {
    names(theta$delta) <- NULL
  } else {
    names(theta$delta) <- c("delta_l", "delta_r") 
  }
  return(theta)
} 