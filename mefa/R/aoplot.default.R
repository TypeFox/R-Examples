aoplot.default <- 
function(x, log=TRUE, xlab, ylab, ...) {
    x <- as.matrix(x)
    A <- AA <- colSums(x)
    O <-  colSums(x > 0)
    lx <- ly <- 0:max(O)
    if (missing(xlab))
        xlab <- "Occurrence"
    if (missing(ylab))
        ylab <- ifelse(log, "log10 Abundance", "Abundance")
    if (log) {
        AA <- log10(AA)
        ly <- log10(ly)
    }
    plot(O, AA, ylog=TRUE, xlab=xlab, ylab=ylab, ...)
    lines(lx, ly, lty=2)
    invisible(cbind(A, O))
}
