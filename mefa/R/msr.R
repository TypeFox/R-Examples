## samples/raw
msr <-
function(x, name, make.unique = FALSE, ...) {
    if (is.null(x$samp))
        stop("'x$samp' must nut be 'NULL'")
    if (missing(name))
        name <- dimnames(x)$taxa
    rval <- data.frame(x$xtab[, name], x$samp, ...)
    if (make.unique && name %in% colnames(x$samp)) {
        colnames(rval) <- c(paste("t.", name, sep = ""),
            paste("s.", colnames(x$samp), sep = ""))
    } else colnames(rval) <- c(name, colnames(x$samp))
    rval
}
