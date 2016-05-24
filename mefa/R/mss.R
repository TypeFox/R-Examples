## samples/summary
mss <-
function(x, name, make.unique = FALSE, ...) {
    if (is.null(x$samp))
        stop("'x$samp' must not be 'NULL'")
    if (missing(name))
        name <- c("s.rich", "s.abu")
    summ <- summary(x)
    rval <- data.frame(data.frame(
        "s.rich" = summ$s.rich,
        "s.abu" = summ$s.abu)[, name],
        x$samp, ...)
    if (make.unique && name %in% colnames(x$samp)) {
        colnames(rval) <- c(name,
            paste("s.", colnames(x$samp), sep = ""))
    } else colnames(rval) <- c(name, colnames(x$samp))
    rval
}
