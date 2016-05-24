# ------------------------------------------
# Authors: Andreas Alfons and Matthias Templ
#          Vienna University of Technology
# ------------------------------------------

#' Weighted mean
#' 
#' Compute the weighted mean.
#' 
#' This is a simple wrapper function calling \code{\link[stats]{weighted.mean}}
#' if sample weights are supplied and \code{\link{mean}} otherwise.
#' 
#' @param x a numeric vector.
#' @param weights an optional numeric vector giving the sample weights.
#' @param na.rm a logical indicating whether missing values in \code{x} should
#' be omitted.
#' 
#' @return The weighted mean of values in \code{x} is returned.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{incMean}}
#' 
#' @keywords survey
#' 
#' @examples
#' data(eusilc)
#' weightedMean(eusilc$eqIncome, eusilc$rb050)
#' 
#' @export

weightedMean <- function(x, weights = NULL, na.rm = FALSE) {
    # initializations
    if (!is.numeric(x)) stop("'x' must be a numeric vector")
    if (is.null(weights)) mean(x, na.rm=na.rm)
    else {
        n <- length(x)
        if (!is.numeric(weights)) stop("'weights' must be a numeric vector")
        else if (length(weights) != n) {
            stop("'weights' must have the same length as 'x'")
        } else if (!all(is.finite(weights))) stop("missing or infinite weights")
        if (any(weights < 0)) warning("negative weights")
        weighted.mean(x, weights, na.rm=na.rm)
    }
}
