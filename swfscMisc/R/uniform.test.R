#' @title Uniform Distribution Test
#' @description Tests whether a histogram is significantly different from a uniform distribution.
#' 
#' @param hist.output output from a call to \code{hist}.
#' @param B number of replicates for chi-squared permutation.
#' 
#' @return result of chi-squared test.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' x.unif <- runif(100)
#' uniform.test(hist(x.unif), B = 1000)
#' x.lnorm <- rlnorm(100)
#' uniform.test(hist(x.lnorm), B = 1000)
#' 
#' @importFrom stats chisq.test
#' @export
#' 
uniform.test <- function(hist.output, B = NULL) {
  break.diff <- diff(hist.output$breaks)
  probs <- break.diff / sum(break.diff)
  if (is.null(B)) {
    chisq.test(x = hist.output$counts, p = probs)
  } else {
    chisq.test(x = hist.output$counts, p = probs, simulate.p.value = TRUE, B = B)
  }
}
