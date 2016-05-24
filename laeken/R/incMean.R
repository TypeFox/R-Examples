# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

#' Weighted mean income
#' 
#' Compute the weighted mean income.
#' 
#' @param inc either a numeric vector giving the (equivalized disposable)
#' income, or (if \code{data} is not \code{NULL}) a character string, an integer
#' or a logical vector specifying the corresponding column of \code{data}.
#' @param weights optional; either a numeric vector giving the personal sample
#' weights, or (if \code{data} is not \code{NULL}) a character string, an
#' integer or a logical vector specifying the corresponding column of
#' \code{data}.
#' @param years optional; either a numeric vector giving the different years of
#' the survey, or (if \code{data} is not \code{NULL}) a character string, an
#' integer or a logical vector specifying the corresponding column of
#' \code{data}.  If supplied, values are computed for each year.
#' @param data an optional \code{data.frame}.
#' @param na.rm a logical indicating whether missing values should be removed.
#' 
#' @return A numeric vector containing the value(s) of the weighted mean income
#' is returned.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{weightedMean}}
#' 
#' @keywords survey
#' 
#' @examples
#' data(eusilc)
#' incMean("eqIncome", weights = "rb050", data = eusilc)
#' 
#' @export

incMean <- function(inc, weights = NULL, years = NULL, 
        data = NULL, na.rm = FALSE) {
	## initializations
    if(!is.null(data)) {
        inc <- data[, inc]
        if(!is.null(weights)) weights <- data[, weights]
        if(!is.null(years)) years <- data[, years]
    }
    # check vectors
    if(!is.numeric(inc)) stop("'inc' must be a numeric vector")
    n <- length(inc)
    if(!is.null(weights) && !is.numeric(weights)) {
        stop("'weights' must be a numeric vector")
    }
    if(!is.null(years) && !is.numeric(years)) {
        stop("'years' must be a numeric vector")
    }
    if(is.null(data)) {  # check vector lengths
        if(!is.null(weights) && length(weights) != n) {
            stop("'weights' must have the same length as 'x'")
        }
        if(!is.null(years) && length(years) != n) {
            stop("'years' must have the same length as 'x'")
        }
    }
    ## computations
    if(is.null(years)) {  # no breakdown
        xn <- weightedMean(inc, weights, na.rm=na.rm)
    } else {  # breakdown by years
        # define wrapper functions
        calcMean <- function(y, inc, weights, years, na.rm) {
            i <- years == y
            weightedMean(inc[i], weights[i], na.rm=na.rm)
        }
        # apply wrapper function
        ys <- sort(unique(years))
        xn <- sapply(ys, calcMean, inc=inc, 
            weights=weights, years=years, na.rm=na.rm)
        names(xn) <- ys  # use years as names
    }
    ## return results
    return(xn)
}
