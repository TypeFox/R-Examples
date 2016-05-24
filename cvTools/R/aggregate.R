# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Aggregate cross-validation results
#' 
#' Compute summary statistics of results from repeated \eqn{K}-fold 
#' cross-validation.  
#' 
#' @method aggregate cv
#' 
#' @param x  an object inheriting from class \code{"cv"} or \code{"cvSelect"} 
#' that contains cross-validation results (note that the latter includes 
#' objects of class \code{"cvTuning"}).
#' @param FUN  a function to compute the summary statistics.
#' @param select  a character, integer or logical vector indicating the columns 
#' of cross-validation results for which to compute the summary statistics.
#' @param \dots  for the \code{"cvTuning"} method, additional arguments to be 
#' passed to the \code{"cvSelect"} method.  Otherwise additional arguments to 
#' be passed to \code{FUN}.
#' 
#' @return 
#' The \code{"cv"} method returns a vector or matrix of aggregated 
#' cross-validation results, depending on whether \code{FUN} returns a single 
#' value or a vector.
#' 
#' For the other methods, a data frame containing the aggregated 
#' cross-validation results for each model is returned.  In the case of the 
#' \code{"cvTuning"} method, the data frame contains the combinations of tuning 
#' parameters rather than a column describing the models.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{cvFit}}, \code{\link{cvSelect}}, 
#' \code{\link{cvTuning}}, \code{\link[stats]{aggregate}}
#' 
#' @example inst/doc/examples/example-aggregate.R
#' 
#' @keywords utilities
#' 
#' @export
#' @import stats

aggregate.cv <- function(x, FUN = mean, select = NULL, ...) {
    if(is.null(CV <- x$reps)) {
        CV <- t(x$cv)  # matrix is required
    }
    if(!is.null(select)) CV <- CV[, select, drop=FALSE]
    apply(CV, 2, FUN=FUN, ...)
}


#' @rdname aggregate.cv
#' @method aggregate cvSelect
#' @export

aggregate.cvSelect <- function(x, FUN = mean, select = NULL, ...) {
    by <- "Fit"
    if(is.null(CV <- x$reps)) CV <- x$cv
    cvNames <- cvNames(x)
    if(is.null(select)) {
        select <- cvNames
    } else if(!is.character(select)) select <- cvNames[select]
    aggregate(CV[, select, drop=FALSE], by=CV[, by, drop=FALSE], FUN=FUN, ...)
}


#' @rdname aggregate.cv
#' @method aggregate cvTuning
#' @export

aggregate.cvTuning <- function(x, ...) {
    # call method for class "cvSelect"
    out <- aggregate.cvSelect(x, ...)
    # replace column specifying the fit by grid of tuning parameters
    cbind(x$tuning, out[, -1, drop=FALSE])
}
