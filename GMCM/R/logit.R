#' Logit and inverse logit transforms
#'
#' The logit transformation (i.e. the log of the odds) and its inverse (also
#' called expit).
#'
#' @aliases logit inv.logit
#' @param p A vector of probabilities.
#' @return \code{logit} returns a vector of the same length as \code{p} with
#'   the log odds of \code{p}.
#' @author Anders Ellern Bilgrau <anders.ellern.bilgrau@@gmail.com>
#' @seealso Used in \code{\link{tt}} and \code{\link{inv.tt}}.
#' @examples
#' p <- runif(100)
#' print(a <- GMCM:::logit(p))
#' p - GMCM:::inv.logit(a)
#' @keywords internal
logit <- function(p) {  # logit function
  log(p/(1-p))
}

