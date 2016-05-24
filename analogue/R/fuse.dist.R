`fuse.dist` <- function(..., weights = NULL) {
    dots <- list(...)
    N <- length(dots)
    if(is.null(weights))
        weights <- rep(1/N, N)
    if(length(weights) != N)
        stop(paste(sQuote("weights"),
                   "must be the same length as number of objects to fuse"))
    if(!isTRUE(all.equal(sum(weights), 1)))
        stop(paste(sQuote("weights"), "must sum to 1"))
    ## reset the storage back to a minimal set
    dots <- lapply(dots, as.dist, diag = FALSE, upper = TRUE)
    ## sanity check to make sure all objects are dist objects
    if(any(!sapply(dots, inherits, c("dist"), USE.NAMES = FALSE)))
        stop(paste("All dissimilarities must be of class",
                   dQuote("dist")))
    ## bind dist vectors to a matrix
    D <- do.call("cbind", dots)
    ## scale each dissimilarity max(d_ij) = 1
    maxs <- apply(D, 2, max, na.rm = TRUE)
    D <- sweep(D, MARGIN = 2, maxs, "/")
    ## weight the N dissimilarities then combines by summation
    retval <- rowSums(D * weights)
    ## cast the return object as a "dist" object
    class(retval) <- "dist"
    attr(retval, "Labels") <- attr(dots[[1]], "Labels")
    attr(retval, "Size") <- attr(dots[[1]], "Size")
    attr(retval, "Diag") <- FALSE
    attr(retval, "Upper") <- FALSE
    attr(retval, "method") <- "fuse"
    attr(retval, "weights") <- weights
    attr(retval, "call") <- match.call()
    return(retval)
}
