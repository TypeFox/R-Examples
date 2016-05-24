#' @title Redo a funreg with different data (for internal use by permutation test)
#' @description For internal use by package functions.  Not
#'  intended to be directly called by the data analyst. 
#'  This function repeats the analysis for an object
#'  of class \code{funreg}, with the same settings 
#'  but with possibly different data.  This is 
#'  convenient in doing resampling techniques like
#'  bootstrapping or permutation testing.
#' @param object A \code{funreg} object.
#' @param id The new values for id (see the \code{funreg} 
#' function documentation)
#' @param response The new values for response 
#' @param time The new values for time 
#' @param other.covariates The new values for other.covariates (which may be NULL)
#' @param x The new values for the functional covariate.
#' @return The \code{funreg} object for the new fitted model.
#'@export
redo.funreg <- function(object,
                        id,
                        response,
                        time,
                        other.covariates,
                        x) {
    stopifnot(class(object)=="funreg");
    return(funreg(id=id,
                  response=response,
                  time=time,
                  x=x,
                  basis.method = object$settings$basis.method,
                  deg = object$settings$deg,
                  deg.penalty = object$settings$deg.penalty,
                  family = object$settings$family,
                  other.covariates=other.covariates,
                  num.bins = object$settings$num.bins,
                  preferred.num.eigenfunctions = object$settings$preferred.num.eigenfunctions,
                  preferred.num.knots.for.beta = object$settings$preferred.num.knots.for.beta,
                  se.method = object$settings$se.method,
                  smoothing.method = object$settings$smoothing.method,
                  times.for.fit.grid = object$times.for.fit.grid ));
}