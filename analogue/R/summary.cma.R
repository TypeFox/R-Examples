###########################################################################
##                                                                       ##
## summary.cma   - extracts and formats close modern analogues           ##
##                                                                       ##
## Author        : Gavin L. Simpson                                      ##
##                                                                       ##
###########################################################################
summary.cma <- function(object, ...) {
    close <- lapply(object$close, function(x) {
        if(length(x) == 0) {
            res <- NA
            names(res) <- "None"
            return(res)
        } else {
            return(x)}
    })
    each.analogs <- sapply(close, length, USE.NAMES = FALSE)
    max.analogs <- max(each.analogs)
    samples <- distances <- matrix(NA, nrow = max.analogs,
                                   ncol = length(close))
    for (i in seq(along = close)) {
        len <- each.analogs[i]
        distances[1:len,i] <- close[[i]]
        samples[1:len,i] <- names(close[[i]])
    }
    rownames(distances) <- rownames(samples) <- 1:max.analogs
    colnames(distances) <- colnames(samples) <- names(close)
    object <- list(close = object$close,
                   call = object$call, cutoff = object$cutoff,
                   quant = object$quant,
                   prob = object$prob,
                   method = object$method,
                   n.analogs = object$n.analogs, distances = distances,
                   samples = samples)
    class(object) <- "summary.cma"
    return(object)
}

print.summary.cma <- function(x,
                              digits = min(3, getOption("digits") - 4),
                              ...) {
    class(x) <- "cma"
    print(x)
    cat("\nDistances:\n\n")
    print(x$distances, digits = digits, na.print = "")
    cat("\nSamples:\n\n")
    print(x$samples, quote = FALSE, right = TRUE, na.print = "")
    invisible(x)
}
