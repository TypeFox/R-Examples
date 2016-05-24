#' Summary results from regsem.
#'
#' @param object An object from regsem.
#' @param ... Other arguments.
#' @export


summary.regsem <- function(object,...){
  fits = fit_indices(object,CV=FALSE)

  TAB <- cbind(convergence = object$convergence,
               df = object$df,
               fit=object$fit,
               rmsea = fits$fit["rmsea"],
               BIC = fits$fit["BIC"])

  ret <- list(call=object$call,
              estimates = object$coefficients,
              returnVals=TAB)

  class(ret) <- "summary.regsem"
  print(ret)
}
