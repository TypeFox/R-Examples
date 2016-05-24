## taxa/raw
mtr <-
function(x, name, make.unique = FALSE, ...) {
    if (is.null(x$taxa))
        stop("'x$taxa' must not be 'NULL'")
    if (missing(name))
        name <- dimnames(x)$samp
    xtab <- if (length(name) == 1)
        data.frame(x$xtab[name,])
        else data.frame(t(x$xtab[name,]))
    rval <- data.frame(xtab, x$taxa, ...)
    if (make.unique && name %in% colnames(x$taxa)) {
        colnames(rval) <- c(paste("s.", name, sep = ""),
            paste("t.", colnames(x$taxa), sep = ""))
    } else colnames(rval) <- c(name, colnames(x$taxa))
    rval
}
