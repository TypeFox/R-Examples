concov <- function (taxa, clustering, digits = 1, width = 5, typical = TRUE, 
    thresh = 10) 
{
    clustering <- clustify(clustering)
    levels <- levels(clustering)
    clustering <- as.integer(clustering)

    x <- const(taxa, clustering)
    y <- importance(taxa, clustering, typical = typical)
    tmp <- NULL
    keep <- apply(as.matrix(x), 1, max) >= thresh/100
    for (i in 1:length(table(clustering))) {
        a <- formatC(as.numeric(x[, i]) * 100, width = 2, format = "d")
        b <- formatC(as.numeric(y[, i]), width = width, digits = digits, 
            format = "f")
        tmp <- cbind(tmp, paste(a, "(", b, ")", sep = ""))
    }
    tmp <- tmp[keep, ]
    tmp <- data.frame(tmp)
    row.names(tmp) <- names(taxa)[keep]
    names(tmp) <- levels
    attr(tmp,'call') <- match.call()
    attr(tmp,'taxa') <- deparse(substitute(taxa))
    attr(tmp,'clustering') <- clustering
    tmp
}

