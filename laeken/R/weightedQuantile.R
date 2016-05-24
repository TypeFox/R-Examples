# ------------------------------------------
# Authors: Andreas Alfons and Matthias Templ
#          Vienna University of Technology
# ------------------------------------------

#' Weighted quantiles
#' 
#' Compute weighted quantiles (Eurostat definition).
#' 
#' The implementation strictly follows the Eurostat definition.
#' 
#' @param x a numeric vector.
#' @param weights an optional numeric vector giving the sample weights.
#' @param probs numeric vector of probabilities with values in \eqn{[0,1]}.
#' @param sorted a logical indicating whether the observations in \code{x} are
#' already sorted.
#' @param na.rm a logical indicating whether missing values in \code{x} should
#' be omitted.
#' 
#' @return A numeric vector containing the weighted quantiles of values in
#' \code{x} at probabilities \code{probs} is returned.  Unlike
#' \code{\link[stats]{quantile}}, this returns an unnamed vector.
#' 
#' @author Andreas Alfons and Matthias Templ
#' 
#' @seealso \code{\link{incQuintile}}, \code{\link{weightedMedian}}
#' 
#' @references Working group on Statistics on Income and Living Conditions
#' (2004) Common cross-sectional EU indicators based on EU-SILC; the gender pay
#' gap.  \emph{EU-SILC 131-rev/04}, Eurostat.
#' 
#' @keywords survey
#' 
#' @examples
#' data(eusilc)
#' weightedQuantile(eusilc$eqIncome, eusilc$rb050)
#' 
#' @export

weightedQuantile <- function(x, weights = NULL, probs = seq(0, 1, 0.25), 
        sorted = FALSE, na.rm = FALSE) {
    # initializations
    if (!is.numeric(x)) stop("'x' must be a numeric vector")
    n <- length(x)
    if (n == 0 || (!isTRUE(na.rm) && any(is.na(x)))) {
        # zero length or missing values
        return(rep.int(NA, length(probs)))
    }
    if (!is.null(weights)) {
        if (!is.numeric(weights)) stop("'weights' must be a numeric vector")
        else if (length(weights) != n) {
            stop("'weights' must have the same length as 'x'")
        } else if (!all(is.finite(weights))) stop("missing or infinite weights")
        if (any(weights < 0)) warning("negative weights")
        if (!is.numeric(probs) || all(is.na(probs)) || 
            isTRUE(any(probs < 0 | probs > 1))) {
            stop("'probs' must be a numeric vector with values in [0,1]")
        }
        if (all(weights == 0)) { # all zero weights
            warning("all weights equal to zero")
            return(rep.int(0, length(probs)))
        }
    }
    # remove NAs (if requested)
    if(isTRUE(na.rm)){
        indices <- !is.na(x)
        x <- x[indices]
        if(!is.null(weights)) weights <- weights[indices]
    }
    # sort values and weights (if requested)
    if(!isTRUE(sorted)) {
#        order <- order(x, na.last=NA)  ## too slow
        order <- order(x)
        x <- x[order]
        weights <- weights[order]  # also works if 'weights' is NULL
    }
    # some preparations
    if(is.null(weights)) rw <- (1:n)/n
    else rw <- cumsum(weights)/sum(weights)
    # obtain quantiles
    q <- sapply(probs, 
        function(p) {
            if (p == 0) return(x[1])
            else if (p == 1) return(x[n])
            select <- min(which(rw >= p))
            if(rw[select] == p) mean(x[select:(select+1)])
            else x[select]
        })
    return(unname(q))
}
