"qqnorm.ltraj" <-
function(y, which=c("dx","dy"), ...)
{
    if (!inherits(y,"ltraj"))
        stop("y should be of class ltraj")
    which <- match.arg(which)
    if (!attr(y, "typeII"))
        stop("y should be of type II (time recorded)")
    opar <- par(mfrow=n2mfrow(length(y)))
    on.exit(par(opar))
    toto <- lapply(y, function(i) {
        tutu <- qqnorm(i[[which]]/sqrt(i$dt), main=attr(i,"burst"),
                       xlab=parse(text=paste(which, "sqrt(dt)", sep="/")),...)
        return(tutu)
    })
    invisible(toto)
}

