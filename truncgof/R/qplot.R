"qplot" <-
function(x, distn, parm, H = NA, plot.it = TRUE, main = "QQ-Plot", 
    xlab = "empirical quantiles", ylab = "theoretical quantiles", ...) 
{
    if (!is.function(try(get(distn), silent = TRUE)))
       stop("'distn' must be a character of a distribution function")
    if (is.na(H)) H <- -Inf
    pdens <- distn
    qdens <- paste("q", substring(distn, 2), sep = "")
    FH <- do.call(pdens, c(list(H), parm))
    n  <- length(x)
    p  <- (n * FH/(1-FH) + c(1:n) - 0.5)/((n * FH/(1-FH) + n))
    sx <- sort(x)
    sy <- do.call(qdens, c(list(p), parm))
    if (plot.it) 
       plot(sx, sy, xlab = xlab, ylab = ylab, main = main, ...)
    abline(0,1)
    invisible(list(x = sx, y = p))
}
