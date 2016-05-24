varExpl <- function(object, ...)
    UseMethod("varExpl")

varExpl.default <- function(object, ...) {
    stop("No default method for 'varExpl()'")
}

varExpl.cca <- function(object, axes = 1L, cumulative = FALSE,
                        pcent = FALSE, ...) {
    if(is.null(object$CCA))
        res <- object$CA$eig[axes]
    else
        res <- object$CCA$eig[axes]
    res <- res / object$tot.chi
    if(cumulative)
        res <- cumsum(res)
    if(pcent)
        res <- 100 * res
    res
}

varExpl.prcurve <- function(object, pcent = FALSE, ...) {
    res <- 1 - object$dist / object$totalDist
    if(pcent)
        res <- 100 * res
    names(res) <- "PrC"
    res
}
