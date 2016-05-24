namedList <- function(...) {
  L <- list(...)
  snm <- sapply(substitute(list(...)), deparse)[-1]
  if (is.null(nm <- names(L))) nm <- snm
  if (any(nonames <- nm == "")) nm[nonames] <- snm[nonames]
  setNames(L,nm)
}

#' @title Control of Joint Mean Covariance Model Fitting
#'
#' @description Construct control structures for joint mean covariance model
#' fitting
#'
#' @param trace whether or not the value of the objective function and the
#' parameters should be print on every trace'th iteration.
#' @param profile whether or not parameters should be estimated sequentially
#' using the idea of profile likelihood.
#' @param ignore.const.term whether or not the constant term should be considered 
#' when calculating log-likelihood and BIC.
#' function
#' @param errorMsg whether or not the error message should be print
#'
#' @export jmcmControl
jmcmControl <- function(trace = FALSE, profile = TRUE, 
                        ignore.const.term = TRUE, errorMsg = FALSE)
{
    structure(namedList(trace, profile, ignore.const.term, errorMsg),
              class = "jmcmControl")
}