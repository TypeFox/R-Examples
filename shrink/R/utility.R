#' Test if \code{\link[rms]{rcs}} from the \code{rms} Package was Used in the Fit.
#'
#' @param fit a fit object of class \code{lm}, \code{glm}, \code{coxph}, or \code{mfp}.
#'
#' @return A list with two elements:
#' \tabular{ll}{
#'     [[1]] \tab \code{\link[rms]{rcs}} is included in formula (right-hand side), but not as a variable name.\cr
#'     [[2]] \tab variables with \code{\link[rms]{rcs}}\cr
#' }
#'
#' @keywords internal
#' @export
isrcs <-
  function(fit)
{
  list(any(all.names(fit$call$formula[3]) == "rcs") && !any(all.vars(fit$call$formula[3]) == "rcs"),
       regexpr(
         pattern = "rcs(", fixed = TRUE, text = names(fit$coefficients)
       ) == TRUE)
}
