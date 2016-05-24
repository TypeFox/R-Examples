# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

#' Weighted median income
#' 
#' Compute the weighted median income.
#' 
#' The implementation strictly follows the Eurostat definition.
#' 
#' @param inc either a numeric vector giving the (equivalized disposable)
#' income, or (if \code{data} is not \code{NULL}) a character string, an integer
#' or a logical vector specifying the corresponding column of \code{data}.
#' @param weights optional; either a numeric vector giving the personal sample
#' weights, or (if \code{data} is not \code{NULL}) a character string, an
#' integer or a logical vector specifying the corresponding column of
#' \code{data}.
#' @param sort optional; either a numeric vector giving the personal IDs to be
#' used as tie-breakers for sorting, or (if \code{data} is not \code{NULL}) a
#' character string, an integer or a logical vector specifying the corresponding
#' column of \code{data}.
#' @param years optional; either a numeric vector giving the different years of
#' the survey, or (if \code{data} is not \code{NULL}) a character string, an
#' integer or a logical vector specifying the corresponding column of
#' \code{data}.  If supplied, values are computed for each year.
#' @param data an optional \code{data.frame}.
#' @param na.rm a logical indicating whether missing values should be removed.
#' 
#' @return A numeric vector containing the value(s) of the weighted median
#' income is returned.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{arpt}}, \code{\link{weightedMedian}}
#' 
#' @references Working group on Statistics on Income and Living Conditions
#' (2004) Common cross-sectional EU indicators based on EU-SILC; the gender pay
#' gap.  \emph{EU-SILC 131-rev/04}, Eurostat.
#' 
#' @keywords survey
#' 
#' @examples
#' data(eusilc)
#' incMedian("eqIncome", weights = "rb050", data = eusilc)
#' 
#' @export

incMedian <- function(inc, weights = NULL, sort = NULL, 
        years = NULL, data = NULL, na.rm = FALSE) {
	## initializations
    if(!is.null(data)) {
        inc <- data[, inc]
        if(!is.null(weights)) weights <- data[, weights]
        if(!is.null(sort)) sort <- data[, sort]
        if(!is.null(years)) years <- data[, years]
    }
    # check vectors
    if(!is.numeric(inc)) stop("'inc' must be a numeric vector")
    n <- length(inc)
    if(!is.null(weights) && !is.numeric(weights)) {
        stop("'weights' must be a numeric vector")
    }
    if(!is.null(sort) && !is.vector(sort) && !is.ordered(sort)) {
        stop("'sort' must be a vector or ordered factor")
    }
    if(!is.null(years) && !is.numeric(years)) {
        stop("'years' must be a numeric vector")
    }
    if(is.null(data)) {  # check vector lengths
        if(!is.null(weights) && length(weights) != n) {
            stop("'weights' must have the same length as 'x'")
        }
        if(!is.null(sort) && length(sort) != n) {
            stop("'sort' must have the same length as 'x'")
        }
        if(!is.null(years) && length(years) != n) {
            stop("'years' must have the same length as 'x'")
        }
    }
    ## sort values and weights
    order <- if(is.null(sort)) order(inc) else order(inc, sort)
    inc <- inc[order]
    weights <- weights[order]  # also works if 'weights' is NULL
    ## computations
    if(is.null(years)) {  # no breakdown
        med <- weightedMedian(inc, weights, sorted=TRUE, na.rm=na.rm)
    } else {  # breakdown by years
        years <- years[order]
        # define wrapper functions
        calcMedian <- function(y, inc, weights, years, na.rm) {
            i <- years == y
            weightedMedian(inc[i], weights[i], sorted=TRUE, na.rm=na.rm)
        }
        # apply wrapper function
        ys <- sort(unique(years))
        med <- sapply(ys, calcMedian, inc=inc, 
            weights=weights, years=years, na.rm=na.rm)
        names(med) <- ys  # use years as names
    }
    ## return results
    return(med)
}
