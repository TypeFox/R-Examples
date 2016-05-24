## taxa/summary
mts <-
function(x, name, make.unique = FALSE, ...) {
    if (is.null(x$taxa))
        stop("'x$taxa' must not be 'NULL'")
    if (missing(name))
        name <- c("t.occ", "t.abu")
    summ <- summary(x)
    rval <- data.frame(data.frame(
        "t.occ" = summ$t.occ,
        "t.abu" = summ$t.abu)[, name],
        x$taxa, ...)
    if (make.unique && name %in% colnames(x$taxa)) {
        colnames(rval) <- c(name,
            paste("t.", colnames(x$taxa), sep = ""))
    } else colnames(rval) <- c(name, colnames(x$taxa))
    rval
}
