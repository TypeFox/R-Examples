##' Extracts fitted values for training set samples from logistic
##' regression models fitted to each group of samples that describe
##' the probability two samples are analogues (from the same group) as
##' a function of dissimilarity between the paired samples.
##'
##' @title Fitted values for the training set from logistic regression
##' models
##' @param object an object of class \code{"logitreg"} resulting from
##' a call to \code{\link{logitreg}}.
##' @param combined logical; should the fitted values for the overall
##' combined analysis be returned.
##' @param ... arguments passed to other methods
##' @return if \code{combined == FALSE} (the default) then a matrix of
##' fitted probabilities, where the rows are the training set samples
##' and the columns the groupings, is returned. If
##' \code{combined == TRUE}, then a list with components \code{"group"}
##' and \code{"combined"}. \code{"group"} is a matrix of fitted
##' probabilities as above. \code{"combined"} is a vector of fitted
##' values for the entire set of pairwise comparisons considered.
##' @author Gavin Simpson
##' @method fitted logitreg
##' @S3method fitted logitreg
`fitted.logitreg` <- function(object, combined = FALSE, ...) {
    n <- length(object$models)
    groups <- sapply(object$models[seq_len(n - 1)], fitted)
    if(combined) {
        combined <- fitted(object$models[[n]])
        out <- list(groups = groups, combined = combined)
    } else {
        out <- groups
    }
    out
}
