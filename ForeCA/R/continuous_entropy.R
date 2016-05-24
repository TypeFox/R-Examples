#' @title Shannon entropy for a continuous pdf
#' 
#' @description
#' Computes the Shannon entropy \eqn{\mathcal{H}(p)} for a continuous 
#' probability density function (pdf) \eqn{p(x)} using numerical integration.
#' 
#' @details
#' The Shannon entropy of a continuous random variable (RV) \eqn{X \sim p(x)} is defined as
#' \deqn{
#' \mathcal{H}(p) = -\int_{-\infty}^{\infty} p(x) \log p(x) d x.
#' }
#' 
#' Contrary to discrete RVs, continuous RVs can have negative entropy (see Examples).
#' 
#' @param pdf R function for the pdf \eqn{p(x)} of a RV \eqn{X \sim p(x)}. This function must
#' be non-negative and integrate to \eqn{1} over the interval [\code{lower}, \code{upper}].
#' @param lower,upper lower and upper integration limit. \code{pdf} must integrate to 
#' \code{1} on this interval.
#' @inheritParams common-arguments
#' @return 
#' scalar; entropy value (real). 
#' 
#' Since \code{continuous_entropy} uses numerical integration (\code{integrate()}) convergence
#' is not garantueed (even if integral in definition of \eqn{\mathcal{H}(p)} exists).
#' Issues a warning if \code{integrate()} does not converge.
#' @keywords math univar
#' @seealso \code{\link{discrete_entropy}}
#' @export
#' @examples
#' # entropy of U(a, b) = log(b - a). Thus not necessarily positive anymore, e.g.
#' continuous_entropy(function(x) dunif(x, 0, 0.5), 0, 0.5) # log2(0.5)
#' 
#' # Same, but for U(-1, 1)
#' my_density <- function(x){
#'   dunif(x, -1, 1)
#' }
#' continuous_entropy(my_density, -1, 1) # = log(upper - lower)
#' 
#' # a 'triangle' distribution
#' continuous_entropy(function(x) x, 0, sqrt(2))
#'
continuous_entropy <- function(pdf, lower, upper, base = 2) {
  
  stopifnot(length(lower) == 1,
            length(upper) == 1,
            is.numeric(lower),
            is.numeric(upper),
            lower < upper)
  
  if (upper == Inf) {
    # set it to some large value < Inf where pdf(x) != 0
    upper.v <- 10^(1:8)
    upper <- upper.v[which.min(pdf(upper.v))[1] - 1]
  }
  if (lower == -Inf) {
    # set it to some large value < Inf where pdf(x) != 0
    lower.v <- -10^(1:8)
    lower <- lower.v[which.min(pdf(lower.v))[1] + 1]
  }

  # check if it is non-negative (at least on a fine grid)
  if (any(pdf(seq(lower, upper, length = 1e3)) < 0)) {
    stop("'pdf' must be non-negative.")
  }
  
  # check if it integrates to 1
  pdf.int <- integrate(pdf, lower, upper)$value
  if (!isTRUE(all.equal(target = 1, current = pdf.int))) {
    stop("'pdf' must integrate to 1; it integrates to ", pdf.int, ".\n")
  }
  
  .pdf_times_log_pdf <- function(xx) {-pdf(xx) * log(pdf(xx), base = base)}

  entropy.int <- integrate(.pdf_times_log_pdf, lower, upper)
  if (entropy.int$message != "OK") {
    warning("Numerical integration didn't not converge properly. ",
            "Use results with caution.  \n ",
            "Message from 'integrate()': ", entropy.int$message)
  }
  entropy.eval <- entropy.int$value
  attr(entropy.eval, "base") <- as.character(base)
  return(entropy.eval)
}
