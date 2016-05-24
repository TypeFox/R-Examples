# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

#' Weighted income quintile
#' 
#' Compute weighted income quintiles.
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
#' @param k a vector of integers between 0 and 5 specifying the quintiles to be
#' computed (0 gives the minimum, 5 the maximum).
#' @param data an optional \code{data.frame}.
#' @param na.rm a logical indicating whether missing values should be removed.
#' 
#' @return A numeric vector (if \code{years} is \code{NULL}) or matrix (if
#' \code{years} is not \code{NULL}) containing the values of the weighted income
#' quintiles specified by \code{k} are returned.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{qsr}}, \code{\link{weightedQuantile}}
#' 
#' @references Working group on Statistics on Income and Living Conditions
#' (2004) Common cross-sectional EU indicators based on EU-SILC; the gender pay
#' gap.  \emph{EU-SILC 131-rev/04}, Eurostat.
#' 
#' @keywords survey
#' 
#' @examples
#' data(eusilc)
#' incQuintile("eqIncome", weights = "rb050", data = eusilc)
#' 
#' @export

incQuintile <- function(inc, weights = NULL, sort = NULL, 
        years = NULL, k = c(1, 4), data = NULL, na.rm = FALSE) {
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
    if(!is.numeric(k) || length(k) == 0 || any(k < -0.5 | k >= 5.5)) {
        stop("'k' must be a vector of integers between 0 and 5")
    } else k <- round(k)
    ## sort values and weights
    order <- if(is.null(sort)) order(inc) else order(inc, sort)
    inc <- inc[order]
    weights <- weights[order]  # also works if 'weights' is NULL
    ## computations
    if(is.null(years)) {  # no breakdown
        q <- weightedQuantile(inc, weights, probs=k/5, sorted=TRUE, na.rm=na.rm)
        names(q) <- k  # use quintile numbers as names
    } else {  # breakdown by years
        years <- years[order]
        # define wrapper functions
        calcQuantile <- function(y, inc, weights, years, k, na.rm) {
            i <- years == y
            weightedQuantile(inc[i], weights[i], 
                probs=k/5, sorted=TRUE, na.rm=na.rm)
        }
        # apply wrapper function
        ys <- sort(unique(years))
        q <- t(sapply(ys, calcQuantile, inc=inc, 
                weights=weights, years=years, k=k, na.rm=na.rm))
        rownames(q) <- ys  # use years as row names
        colnames(q) <- k   # use quintile numbers as column names
    }
    ## return results
    return(q)
}
