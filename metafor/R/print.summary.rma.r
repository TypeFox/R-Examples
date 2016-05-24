print.summary.rma <-
function (x, digits, showfit = TRUE, signif.stars = getOption("show.signif.stars"), 
    signif.legend = signif.stars, ...) 
{
    if (!is.element("summary.rma", class(x))) 
        stop("Argument 'x' must be an object of class \"summary.rma\".")
    if (missing(digits)) 
        digits <- x$digits
    class(x) <- class(x)[-1]
    print(x, digits = digits, showfit = showfit, signif.stars = signif.stars, 
        signif.legend = signif.legend, ...)
    invisible()
}
