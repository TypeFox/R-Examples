#' @rdname logit
#' @param a A vector of real values.
#' @return \code{inv.logit} returns a vector of the same length as \code{a} of the
#'   inverse logit transformed values. This function is also known as the
#'   expit-function.
#' @keywords internal
inv.logit <- function(a) { # inverse logit function
  ans <- exp(a)/(1+exp(a))
  ans[is.nan(ans)] <- 1  # If a == Inf then ans should be 1
  return(ans)
}
