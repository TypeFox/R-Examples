"plot.fipati" <- function(x, scale, warn = TRUE, ...)
{
    ## Verifications
    if (!inherits(x, "fipati"))
        stop("x should be of class 'fipati'")

    ## graphical settings
    opar <- par(mfrow=n2mfrow(length(x)))

    ## radius
    att <- attr(x, "radii")

    ## Does the specified scale exists?
    ind <- which.min(abs(att-scale))
    if (warn)
        if (abs(att[ind] - scale) > 1e-7)
            warning(paste("No radius equal to ",scale,
                          ", displayed radius is ", att[ind], sep =""))

    ## Draws the plot
    lapply(1:length(x), function(i) {
        u <- x[[i]]
        plot(attr(u,"date"), u[,ind], main = attr(u, "burst"), ylab="FPT",...)
        lines(attr(u,"date"), u[,ind])})
    par(opar)
}

