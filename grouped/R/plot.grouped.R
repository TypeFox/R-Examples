"plot.grouped" <-
function(x, B = 100, sub.caption = deparse(formula(x)), ...){
    if(!inherits(x, "grouped"))
        stop("Use only with 'grouped' objects.\n")
    sub. <- x$call
    r <- residuals(x, B = B)$residuals
    fits <- fitted(x)
    plot(fits, r, ylab = "Residuals", xlab = "Fitted values",
                main = "Residuals vs Fitted", ...)
    abline(h = 0.0, lty = 2)
    title(sub = sub.caption, ...)
    invisible()
}

