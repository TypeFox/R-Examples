#' @title Computes support for skewed Lambert W x F distributions
#' 
#' @description If the input \eqn{X \sim F} has support on the entire real line
#'     \eqn{(-\infty, \infty)}, then the skewed Lambert W \eqn{\times} F
#'     distribution has truncated support \eqn{[a,b]}, \eqn{a,b \in R \cup \pm
#'     \infty} depending on \eqn{\boldsymbol \beta} and (the sign of)
#'     \eqn{\gamma}.
#' 
#' For scale-families no truncation occurs.
#' 
#' @name get_support
#' @inheritParams common-arguments
#' @inheritParams loglik-LambertW-utils
#' @return A vector of length 2 with names \code{'lower'} and \code{'upper'}.
#' @details Half-open interval on the real line (if \eqn{\gamma \neq 0}) for
#'     input with support on the entire real line. For \eqn{\gamma = 0} the
#'     support of Y is the same as for X. Heavy-tail Lambert W RVs are not
#'     affected by truncated support (for \eqn{\delta \geq 0}); thus support is
#'     \code{c(lower = -Inf, upper = Inf)}.
#' @keywords math
#' @param input.bounds interval; the bounds of the input distribution. If
#'     \code{is.non.negative = FALSE}, then it will adjust it to \code{c(0,
#'     Inf)}; also useful for bounded input distributions, such as
#'     \code{"unif"}.
#' @export
#' @examples
#' 
#' get_support(c(mu_x = 0, sigma_x = 1, gamma = 0)) # as gamma = 0
#' # truncated on the left since gamma > 0
#' get_support(c(mu_x = 0, sigma_x = 1, gamma = 0.1)) 
#' 
#' # no truncation for heavy tail(s)
#' get_support(c(mu_x = 0, sigma_x = 1, delta = 0.1))

get_support <- function(tau, is.non.negative = FALSE, 
                        input.bounds = c(-Inf, Inf)) {
  stopifnot(is.logical(is.non.negative),
            is.numeric(input.bounds),
            length(input.bounds) == 2,
            input.bounds[1] <= input.bounds[2])
  
  tau <- complete_tau(tau)
  check_tau(tau)
  if (is.na(tau["gamma"])) {
    tau["gamma"] <- 0
  } 

  if (!(identical(input.bounds, c(-Inf, Inf)) || 
          identical(input.bounds, c(0, Inf)))) {
    rv.support <- get_output(input.bounds, tau)
  } else {
    if (is.non.negative) {
      input.bounds <- c(0, Inf)
    }
    if (is.non.negative) {
      # assumes that gamma, delta, alpha >= 0
      rv.support <- c(0, Inf)
    } else {
      if (tau["gamma"] == 0) {
        rv.support <- c(-Inf, Inf) 
      } else {
        bb <- normalize_by_tau(1/(-tau["gamma"] * exp(1)), tau, inverse = TRUE)
        if (tau["gamma"] > 0) {
          rv.support <- c(bb, Inf)
        } else if (tau["gamma"] < 0) {
          rv.support <- c(-Inf, bb)
        }
      }
    }
  }
  names(rv.support) <- c("lower", "upper")
  return(rv.support)
}
