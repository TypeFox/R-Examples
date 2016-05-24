#' @export
probplot <-
function (x, qdist = qnorm, probs = NULL, line = TRUE, xlab = NULL, 
    ylab = "Percentile", ...) 
{
    DOTARGS <- as.list(substitute(list(...)))[-1]
    DOTARGS <- paste(names(DOTARGS), DOTARGS, sep = "=", collapse = ", ")
    xlab = deparse(substitute(x))
    x <- sort(x)
    QNAME <- deparse(substitute(qdist))
    DOTS <- list(...)
    qdist <- match.fun(qdist)
    QFUN <- function(p) {
        args = DOTS
        args$p = p
        do.call("qdist", args)
    }
    y <- QFUN(ppoints(length(x)))
    if (is.null(probs)) {
        probs <- c(0.01, 0.05, seq(0.1, 0.9, by = 0.1), 0.95, 
            0.99)
        if (length(x) >= 1000) 
            probs <- c(0.001, probs, 0.999)
    }
    qprobs <- QFUN(probs)
    plot(x, y, axes = FALSE, type = "n", ylim = range(c(y, qprobs)), 
        xlab = xlab, ylab = ylab)
    box()
    abline(h = qprobs, col = "grey")
    axis(1)
    axis(2, at = qprobs, labels = 100 * probs)
    points(x, y)
    QTEXT <- paste("Quantile: ", QNAME, sep = "")
    if (nchar(DOTARGS)) 
        QTEXT <- paste(QTEXT, DOTARGS, sep = ", ")
    mtext(QTEXT, side = 1, line = 3, adj = 1)
    xl <- quantile(x, c(0.25, 0.75))
    yl <- qdist(c(0.25, 0.75), ...)
    slope <- diff(yl)/diff(xl)
    int <- yl[1] - slope * xl[1]
    if (line) {
        abline(int, slope, col = "lightskyblue3")
    }
    z <- list(qdist = QFUN, int = int, slope = slope)
    class(z) <- "probplot"
    invisible(z)
}
