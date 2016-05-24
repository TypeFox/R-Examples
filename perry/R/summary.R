# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Summarize resampling-based prediction error results
#' 
#' Produce a summary of resampling-based prediction error results.  
#' 
#' @method summary perry
#' 
#' @param object  an object inheriting from class \code{"perry"} or 
#' \code{"perrySelect"} that contains prediction error results (note that the 
#' latter includes objects of class \code{"perryTuning"}).
#' @param \dots  currently ignored.
#' 
#' @return 
#' An object of class \code{"summary.perry"}, \code{"summary.perrySelect"} or 
#' \code{"summary.perryTuning"}, depending on the class of \code{object}.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{perryFit}}, \code{\link{perrySelect}}, 
#' \code{\link{perryTuning}}, \code{\link{summary}}
#' 
#' @example inst/doc/examples/example-summary.R
#' 
#' @keywords utilities
#' 
#' @export

summary.perry <- function(object, ...) {
    pe <- aggregate(object, summary)
    out <- list(pe=pe, splits=object$splits)
    class(out) <- "summary.perry"
    out
}


#' @rdname summary.perry
#' @method summary perrySelect
#' @export

summary.perrySelect <- function(object, ...) {
    pe <- aggregate(object, summary)
    out <- list(pe=pe, splits=object$splits, best=object$best)
    class(out) <- "summary.perrySelect"
    out
}


#' @rdname summary.perry
#' @method summary perryTuning
#' @export

summary.perryTuning <- function(object, ...) {
    out <- summary.perrySelect(object, ...)
    out <- list(pe=out$pe, splits=object$splits, 
        best=out$best, tuning=object$tuning)
    class(out) <- c("summary.perryTuning", class(out))
    out
}
