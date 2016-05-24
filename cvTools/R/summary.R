# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Summarize cross-validation results
#' 
#' Produce a summary of results from (repeated) \eqn{K}-fold cross-validation.  
#' 
#' @method summary cv
#' 
#' @param object  an object inheriting from class \code{"cv"} or 
#' \code{"cvSelect"} that contains cross-validation results (note that the 
#' latter includes objects of class \code{"cvTuning"}).
#' @param \dots  currently ignored.
#' 
#' @return 
#' An object of class \code{"summary.cv"}, \code{"summary.cvSelect"} or 
#' \code{"summary.cvTuning"}, depending on the class of \code{object}.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{cvFit}}, \code{\link{cvSelect}}, 
#' \code{\link{cvTuning}}, \code{\link{summary}}
#' 
#' @example inst/doc/examples/example-summary.R
#' 
#' @keywords utilities
#' 
#' @export

summary.cv <- function(object, ...) {
    cv <- aggregate(object, summary)
    out <- list(n=object$n, K=object$K, R=object$R, cv=cv)
    class(out) <- "summary.cv"
    out
}


#' @rdname summary.cv
#' @method summary cvSelect
#' @export

summary.cvSelect <- function(object, ...) {
    cv <- aggregate(object, summary)
    out <- list(n=object$n, K=object$K, R=object$R, best=object$best, cv=cv)
    class(out) <- "summary.cvSelect"
    out
}


#' @rdname summary.cv
#' @method summary cvTuning
#' @export

summary.cvTuning <- function(object, ...) {
    out <- summary.cvSelect(object, ...)
    out <- list(n=out$n, K=out$K, R=out$R, tuning=object$tuning, 
        best=out$best, cv=out$cv)
    class(out) <- c("summary.cvTuning", class(out))
    out
}
