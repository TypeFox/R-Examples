"qqchi.ltraj" <- function (y, xlab = "Theoretical Quantiles",
                           ylab = "Sample Quantiles (Distances)",
                           ...)
{
    if (!inherits(y,"ltraj"))
        stop("y should be of class ltraj")
    if (!attr(y, "typeII"))
        stop("y should be of type II (time recorded)")
    opar <- par(mfrow=n2mfrow(length(y)))
    on.exit(par(opar))
    toto <- lapply(y, function(i) {
        ys <- i$dist/sqrt(i$dt)
        tutu <- qqchi(ys,xlab=xlab,ylab=ylab,main=attr(i,"burst"),...)
        return(tutu)
    })
    invisible(toto)
}

