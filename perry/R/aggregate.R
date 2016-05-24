# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Aggregate resampling-based prediction error results
#' 
#' Compute summary statistics of resampling-based prediction error results.  
#' 
#' @method aggregate perry
#' 
#' @param x  an object inheriting from class \code{"perry"} or 
#' \code{"perrySelect"} that contains prediction error results (note that the 
#' latter includes objects of class \code{"perryTuning"}).
#' @param FUN  a function to compute the summary statistics.
#' @param select  a character, integer or logical vector indicating the columns 
#' of prediction error results for which to compute the summary statistics.
#' @param \dots  for the \code{"perryTuning"} method, additional arguments to 
#' be passed to the \code{"perrySelect"} method.  Otherwise additional 
#' arguments to be passed to \code{FUN}.
#' 
#' @return 
#' The \code{"perry"} method returns a vector or matrix of aggregated 
#' prediction error results, depending on whether \code{FUN} returns a single 
#' value or a vector.
#' 
#' For the other methods, a data frame containing the aggregated 
#' prediction error results for each model is returned.  In the case of the 
#' \code{"perryTuning"} method, the data frame contains the combinations of 
#' tuning parameters rather than a column describing the models.
#' 
#' @note Duplicate indices in \code{subset} or \code{select} are removed such 
#' that all models and prediction error results are unique.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{perryFit}}, \code{\link{perrySelect}}, 
#' \code{\link{perryTuning}}, \code{\link[stats]{aggregate}}
#' 
#' @example inst/doc/examples/example-aggregate.R
#' 
#' @keywords utilities
#' 
#' @export
#' @import stats

aggregate.perry <- function(x, FUN = mean, select = NULL, ...) {
    if(is.null(PE <- x$reps)) PE <- t(x$pe)  # matrix is required
    if(!is.null(select)) {
        select <- checkSelect(select, peNames(x))
        PE <- PE[, select, drop=FALSE]
    }
    apply(PE, 2, FUN=FUN, ...)
}


#' @rdname aggregate.perry
#' @method aggregate perrySelect
#' @export

aggregate.perrySelect <- function(x, FUN = mean, select = NULL, ...) {
    by <- "Fit"
    if(is.null(PE <- x$reps)) PE <- x$pe
    peNames <- peNames(x)
    select <- if(is.null(select)) peNames else checkSelect(select, peNames)
    aggregate(PE[, select, drop=FALSE], by=PE[, by, drop=FALSE], FUN=FUN, ...)
}


#' @rdname aggregate.perry
#' @method aggregate perryTuning
#' @export

aggregate.perryTuning <- function(x, ...) {
    # call method for class "perrySelect"
    out <- aggregate.perrySelect(x, ...)
    # replace column specifying the fit by grid of tuning parameters
    cbind(x$tuning, out[, -1, drop=FALSE])
}
